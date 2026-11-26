import os
import sys
import io
import uuid
import json
import sqlite3
import datetime
import pandas as pd
from flask import Flask, render_template, request, jsonify, send_file, redirect, url_for, flash, Response
from werkzeug.utils import secure_filename
from config import Config

# Modules
from modules.preprocessing import validate_and_preprocess
from modules.data_upload import load_geo_dataset, get_curated_geo_list
from modules.differential_expression import run_differential_expression, benjamini_hochberg
from modules.pathway_enrichment import run_pathway_enrichment
from modules.network_analysis import build_and_analyze_network
from modules.biomarker_validation import run_biomarker_classification
from modules.ai_assistant import ask_ai_assistant, get_gemini_api_key
from modules.report_generator import generate_html_report
from database.db import (
    init_db, save_analysis_run, get_dashboard_stats, get_all_analyses,
    get_analysis_by_id, delete_analysis, save_chat_message, get_chat_history,
    get_db_connection, save_session_snapshot, load_session_snapshot
)

app = Flask(__name__)
app.config.from_object(Config)

# WSGI middleware: gracefully handle any Vercel internal prefix rewrites
class VercelPathFix:
    def __init__(self, wsgi_app):
        self.wsgi_app = wsgi_app

    def __call__(self, environ, start_response):
        path = environ.get('PATH_INFO', '')
        for prefix in ['/api/index.py', '/api/index']:
            if path.startswith(prefix):
                new_path = path[len(prefix):]
                if not new_path.startswith('/'):
                    new_path = '/' + new_path
                environ['PATH_INFO'] = new_path
                break
        return self.wsgi_app(environ, start_response)

app.wsgi_app = VercelPathFix(app.wsgi_app)

# Ensure required directories exist
for folder in [Config.UPLOAD_FOLDER, Config.RESULTS_FOLDER, Config.SAMPLE_DATA_FOLDER, Config.R_SCRIPTS_FOLDER]:
    os.makedirs(folder, exist_ok=True)

# Initialize database schema
init_db()

# In-memory store for active session states (keyed by session_id)
SESSIONS = {}

def _read_text_limited(filepath, max_bytes=1500000):
    """Read small CSV/text file if under max_bytes for snapshot serialization."""
    if not filepath or not os.path.exists(filepath):
        return None
    try:
        if os.path.getsize(filepath) > max_bytes:
            return None
        with open(filepath, 'r', encoding='utf-8') as f:
            return f.read()
    except Exception:
        return None

def build_session_snapshot(sess, include_files=True, include_plots=True):
    """Create a JSON-serializable snapshot of a session for SQLite and client persistence."""
    de_summary = sess.get('de_summary')
    if isinstance(de_summary, dict):
        de_sum_copy = dict(de_summary)
        if not include_plots:
            de_sum_copy.pop('heatmap_plot', None)
            de_sum_copy.pop('pca_plot', None)
    else:
        de_sum_copy = None

    ml_metrics = sess.get('ml_metrics')
    if isinstance(ml_metrics, dict):
        ml_copy = {
            'success': ml_metrics.get('success', True),
            'cv_method': ml_metrics.get('cv_method'),
            'metrics': ml_metrics.get('metrics'),
            'independent_validation': ml_metrics.get('independent_validation'),
        }
        if include_plots:
            ml_copy['roc_plot'] = ml_metrics.get('roc_plot')
            ml_copy['confusion_matrix_plot'] = ml_metrics.get('confusion_matrix_plot')
            ml_copy['expression_comparison_plot'] = ml_metrics.get('expression_comparison_plot')
    else:
        ml_copy = None

    snap = {
        'session_id': sess.get('session_id'),
        'dataset_name': sess.get('dataset_name'),
        'preproc_info': sess.get('preproc_info'),
        'expr_preview': sess.get('expr_preview'),
        'meta_preview': sess.get('meta_preview'),
        'pca_preview': sess.get('pca_preview') if include_plots else None,
        'de_summary': de_sum_copy,
        'volcano_plot': sess.get('volcano_plot') if include_plots else None,
        'pathway_results': sess.get('pathway_results'),
        'pathway_bar_plot': sess.get('pathway_bar_plot') if include_plots else None,
        'network_summary': sess.get('network_summary'),
        'hub_genes': sess.get('hub_genes'),
        'network_edges': sess.get('network_edges'),
        'network_json': sess.get('network_json') if include_plots else None,
        'ml_metrics': ml_copy,
    }

    if include_files:
        snap['processed_expr_csv'] = _read_text_limited(sess.get('processed_expr_path'))
        snap['processed_meta_csv'] = _read_text_limited(sess.get('processed_meta_path'))
        de_file = (de_summary or {}).get('de_file') if isinstance(de_summary, dict) else None
        if de_file and os.path.exists(de_file):
            snap['de_csv_text'] = _read_text_limited(de_file)
        elif sess.get('de_df') is not None and hasattr(sess['de_df'], 'to_csv'):
            snap['de_csv_text'] = sess['de_df'].to_csv(index=False)
        snap['hubs_csv_text'] = _read_text_limited(os.path.join(sess['dir'], "ranked_hub_genes.csv"))
        snap['edges_csv_text'] = _read_text_limited(os.path.join(sess['dir'], "network_edges.csv"))
        snap['metrics_csv_text'] = _read_text_limited(os.path.join(sess['dir'], "validation_metrics.csv"))

    return snap

def restore_session_from_snapshot(session_id, snapshot):
    """Rehydrate in-memory SESSIONS entry and /tmp working files from a snapshot dict."""
    if not session_id or not isinstance(snapshot, dict):
        return None
    sess_dir = os.path.join(Config.RESULTS_FOLDER, session_id)
    os.makedirs(sess_dir, exist_ok=True)

    existing = SESSIONS.get(session_id) or {
        'session_id': session_id,
        'dir': sess_dir,
        'dataset_name': 'Custom Upload',
        'preproc_info': None,
        'expr_preview': None,
        'meta_preview': None,
        'pca_preview': None,
        'de_summary': None,
        'de_df': None,
        'volcano_plot': None,
        'pathway_results': None,
        'pathway_bar_plot': None,
        'network_summary': None,
        'hub_genes': None,
        'network_edges': None,
        'network_json': None,
        'ml_metrics': None,
        'processed_expr_path': None,
        'processed_meta_path': None
    }

    for key in [
        'dataset_name', 'preproc_info', 'expr_preview', 'meta_preview',
        'pca_preview', 'de_summary', 'volcano_plot', 'pathway_results',
        'pathway_bar_plot', 'network_summary', 'hub_genes', 'network_edges',
        'network_json', 'ml_metrics'
    ]:
        if snapshot.get(key) is not None:
            existing[key] = snapshot[key]

    # Reconstruct processed CSVs in sess_dir if needed
    if snapshot.get('processed_expr_csv'):
        p_expr = os.path.join(sess_dir, "processed_expression.csv")
        if not os.path.exists(p_expr):
            with open(p_expr, 'w', encoding='utf-8') as f:
                f.write(snapshot['processed_expr_csv'])
        existing['processed_expr_path'] = p_expr

    if snapshot.get('processed_meta_csv'):
        p_meta = os.path.join(sess_dir, "processed_metadata.csv")
        if not os.path.exists(p_meta):
            with open(p_meta, 'w', encoding='utf-8') as f:
                f.write(snapshot['processed_meta_csv'])
        existing['processed_meta_path'] = p_meta

    if snapshot.get('de_csv_text'):
        de_sum = existing.get('de_summary') or {}
        case_g = de_sum.get('case_group', 'Case')
        ctrl_g = de_sum.get('control_group', 'Control')
        de_file = os.path.join(sess_dir, f"de_results_{case_g}_vs_{ctrl_g}.csv")
        if not os.path.exists(de_file):
            with open(de_file, 'w', encoding='utf-8') as f:
                f.write(snapshot['de_csv_text'])
        if isinstance(existing.get('de_summary'), dict):
            existing['de_summary']['de_file'] = de_file
        try:
            existing['de_df'] = pd.read_csv(io.StringIO(snapshot['de_csv_text']))
        except Exception:
            pass
    elif existing.get('de_df') is None and isinstance(existing.get('de_summary'), dict):
        top_genes = existing['de_summary'].get('top_de_genes')
        if top_genes:
            existing['de_df'] = pd.DataFrame(top_genes)

    for fname, snap_key in [
        ("ranked_hub_genes.csv", "hubs_csv_text"),
        ("network_edges.csv", "edges_csv_text"),
        ("validation_metrics.csv", "metrics_csv_text")
    ]:
        if snapshot.get(snap_key):
            fpath = os.path.join(sess_dir, fname)
            if not os.path.exists(fpath):
                with open(fpath, 'w', encoding='utf-8') as f:
                    f.write(snapshot[snap_key])

    SESSIONS[session_id] = existing
    try:
        save_session_snapshot(session_id, snapshot)
    except Exception:
        pass
    return existing

def sync_analysis_history(sess):
    """Save/update analysis run metadata in SQLite and return the clean record dict."""
    de_sum = sess.get('de_summary') or {}
    preproc = sess.get('preproc_info') or {}
    ml_res = sess.get('ml_metrics') or {}
    ml_m = ml_res.get('metrics') or {}
    hubs = sess.get('hub_genes') or []
    paths = sess.get('pathway_results') or []

    top_hubs_str = ", ".join([h.get('Gene', '') for h in hubs[:4] if h.get('Gene')])
    top_paths_str = ", ".join([p.get('term', '') for p in paths[:3] if p.get('term')])

    record = {
        'id': sess['session_id'],
        'timestamp': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'dataset_name': sess.get('dataset_name', 'Unknown'),
        'input_type': preproc.get('inferred_type', 'N/A'),
        'total_genes': preproc.get('num_genes', 0),
        'total_samples': preproc.get('num_samples', 0),
        'case_group': de_sum.get('case_group', 'Case'),
        'control_group': de_sum.get('control_group', 'Control'),
        'de_method': de_sum.get('engine', 'Statistical Engine'),
        'fdr_threshold': de_sum.get('fdr_threshold', 0.05),
        'log2fc_threshold': de_sum.get('log2fc_threshold', 1.0),
        'sig_genes_count': de_sum.get('significant_total', 0),
        'n_up': de_sum.get('upregulated', 0),
        'n_down': de_sum.get('downregulated', 0),
        'top_hub_genes': top_hubs_str,
        'top_pathways': top_paths_str,
        'classifier_type': ml_m.get('classifier', ''),
        'roc_auc': ml_m.get('roc_auc', 0.0),
        'accuracy': ml_m.get('accuracy', 0.0),
        'precision_score': ml_m.get('precision', 0.0),
        'recall_score': ml_m.get('recall', 0.0),
        'f1_score': ml_m.get('f1_score', 0.0),
        'details': {
            'confusion_matrix': ml_m.get('confusion_matrix', {})
        }
    }

    try:
        save_analysis_run(record)
    except Exception as e:
        print(f"[DB Warning] Could not save analysis run: {e}")

    try:
        save_session_snapshot(sess['session_id'], build_session_snapshot(sess))
    except Exception as e:
        print(f"[DB Warning] Could not save session snapshot: {e}")

    return record

def get_or_create_session(session_id=None):
    """Retrieve existing session (or restore from SQLite snapshot) or instantiate a new one."""
    if session_id and session_id not in SESSIONS:
        try:
            snap = load_session_snapshot(session_id)
            if snap:
                return restore_session_from_snapshot(session_id, snap)
        except Exception:
            pass

    if not session_id or session_id not in SESSIONS:
        new_id = (session_id if session_id else str(uuid.uuid4())[:8])
        sess_dir = os.path.join(Config.RESULTS_FOLDER, new_id)
        os.makedirs(sess_dir, exist_ok=True)
        SESSIONS[new_id] = {
            'session_id': new_id,
            'dir': sess_dir,
            'dataset_name': 'Custom Upload',
            'preproc_info': None,
            'expr_preview': None,
            'meta_preview': None,
            'pca_preview': None,
            'de_summary': None,
            'de_df': None,
            'volcano_plot': None,
            'pathway_results': None,
            'pathway_bar_plot': None,
            'network_summary': None,
            'hub_genes': None,
            'network_edges': None,
            'network_json': None,
            'ml_metrics': None,
            'processed_expr_path': None,
            'processed_meta_path': None
        }
        return SESSIONS[new_id]
    return SESSIONS[session_id]

# ----------------- PAGE NAVIGATION ROUTES -----------------

@app.route('/')
@app.route('/dashboard')
def dashboard():
    """Main Dashboard with live calculated statistics from database."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    stats = get_dashboard_stats()
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('dashboard.html', stats=stats, r_status=r_status, session=sess)

@app.route('/upload')
def upload_page():
    """Module 1: Data Upload & Public GEO Dataset Access."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    geo_list = get_curated_geo_list()
    return render_template('upload.html', session=sess, r_status=r_status, geo_list=geo_list)

@app.route('/preprocessing')
def preprocessing_page():
    """Module 1: Quality Control & Preprocessing Inspection."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('preprocessing.html', session=sess, r_status=r_status)

@app.route('/differential-expression')
def differential_expression_page():
    """Module 2: Differential Expression (Limma / DESeq2)."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('differential_expression.html', session=sess, r_status=r_status)

@app.route('/pathway-enrichment')
def pathway_enrichment_page():
    """Module 2: Biological Pathway Enrichment (GSEApy Enrichr)."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('pathway_enrichment.html', session=sess, r_status=r_status)

@app.route('/network-analysis')
def network_analysis_page():
    """Module 3: Protein-Protein Interaction Network & Hub Genes."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('network_analysis.html', session=sess, r_status=r_status)

@app.route('/biomarker-validation')
def biomarker_validation_page():
    """Module 4: Machine Learning Biomarker Validation & Classification."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('validation.html', session=sess, r_status=r_status)

@app.route('/ai-assistant')
def ai_assistant_page():
    """Module 5: BioPrism AI-Assisted Result Interpretation."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    gemini_key = get_gemini_api_key()
    return render_template('ai_assistant.html', session=sess, r_status=r_status, gemini_active=bool(gemini_key))

@app.route('/reports')
def reports_page():
    """Module 6: Visualization, Capstone Reporting & PDF Exports."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('reports.html', session=sess, r_status=r_status)

@app.route('/history')
def history_page():
    """Analysis History & Provenance Tracking."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    analyses = get_all_analyses()
    return render_template('history.html', session=sess, r_status=r_status, analyses=analyses)

@app.route('/settings')
def settings_page():
    """Application Settings, Runtime Diagnostics & Database Management."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    sys_info = {
        'python_version': sys.version.split()[0],
        'rscript_path': Config.RSCRIPT_PATH,
        'r_ready': r_status,
        'gemini_configured': bool(get_gemini_api_key()),
        'db_path': Config.DATABASE_PATH
    }
    return render_template('settings.html', session=sess, r_status=r_status, sys_info=sys_info)

@app.route('/pipeline')
def pipeline():
    """Interactive Continuous Pipeline Workspace."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('pipeline.html', session=sess, r_status=r_status)

@app.route('/help')
def help_page():
    """Methodology, Capstone Defense Guide & FAQs."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    r_status = os.path.exists(Config.RSCRIPT_PATH)
    return render_template('help.html', session=sess, r_status=r_status)

# ----------------- MODULE 1: UPLOAD & PREPROCESSING APIS -----------------

def _reset_downstream_session(sess):
    """Clear downstream analytical state when a new dataset is loaded into a session."""
    for key in [
        'de_summary', 'de_df', 'volcano_plot',
        'pathway_results', 'pathway_bar_plot',
        'network_summary', 'hub_genes', 'network_edges', 'network_json',
        'ml_metrics'
    ]:
        sess[key] = None

@app.route('/api/session/restore', methods=['POST'])
def api_session_restore():
    """Restore session state from client snapshot into active server container & SQLite."""
    data = request.get_json() or {}
    session_id = data.get('session_id')
    snapshot = data.get('session_snapshot')
    if not session_id or not isinstance(snapshot, dict):
        return jsonify({'success': False, 'message': 'Missing session_id or snapshot.'})
    sess = restore_session_from_snapshot(session_id, snapshot)
    return jsonify({'success': True, 'session_id': sess['session_id']})

@app.route('/api/load-demo', methods=['POST'])
def api_load_demo():
    """Load built-in verified demo datasets."""
    data = request.get_json() or {}
    demo_type = data.get('type', 'microarray')
    session_id = data.get('session_id')
    sess = get_or_create_session(session_id)

    if demo_type == 'rnaseq':
        expr_path = os.path.join(Config.SAMPLE_DATA_FOLDER, 'demo_rnaseq_counts.csv')
        meta_path = os.path.join(Config.SAMPLE_DATA_FOLDER, 'demo_rnaseq_metadata.csv')
        dataset_name = "Synthetic RNA-seq Tumor vs Normal (Raw Counts)"
        inferred = "rnaseq_counts"
    else:
        expr_path = os.path.join(Config.SAMPLE_DATA_FOLDER, 'demo_microarray_expression.csv')
        meta_path = os.path.join(Config.SAMPLE_DATA_FOLDER, 'demo_microarray_metadata.csv')
        dataset_name = "Oncogenomics Microarray Demo (Normalized Log2)"
        inferred = "microarray_log2"

    result = validate_and_preprocess(expr_path, meta_path, input_type=inferred, output_dir=sess['dir'])
    if result['success']:
        _reset_downstream_session(sess)
        sess['dataset_name'] = dataset_name
        sess['preproc_info'] = result['info']
        sess['expr_preview'] = result.get('expr_preview')
        sess['meta_preview'] = result.get('meta_preview')
        sess['pca_preview'] = result.get('pca_preview')
        sess['processed_expr_path'] = result['processed_expr_path']
        sess['processed_meta_path'] = result['processed_meta_path']
        result['session_id'] = sess['session_id']
        result['dataset_name'] = dataset_name
        snap = build_session_snapshot(sess)
        try:
            save_session_snapshot(sess['session_id'], snap)
        except Exception:
            pass
        result['session_snapshot'] = snap
        result['reset_snapshot'] = True

    return jsonify(result)

@app.route('/api/load-geo', methods=['POST'])
def api_load_geo():
    """Load public benchmark GEO transcriptomics cohort."""
    data = request.get_json() or {}
    geo_id = data.get('geo_id', 'GSE42568')
    session_id = data.get('session_id')
    sess = get_or_create_session(session_id)

    geo_res = load_geo_dataset(geo_id, sess['dir'])
    if not geo_res['success']:
        return jsonify(geo_res)

    result = validate_and_preprocess(geo_res['expr_path'], geo_res['meta_path'], input_type=geo_res['type'], output_dir=sess['dir'])
    if result['success']:
        _reset_downstream_session(sess)
        sess['dataset_name'] = geo_res['title']
        sess['preproc_info'] = result['info']
        sess['expr_preview'] = result.get('expr_preview')
        sess['meta_preview'] = result.get('meta_preview')
        sess['pca_preview'] = result.get('pca_preview')
        sess['processed_expr_path'] = result['processed_expr_path']
        sess['processed_meta_path'] = result['processed_meta_path']
        result['session_id'] = sess['session_id']
        result['dataset_name'] = geo_res['title']
        snap = build_session_snapshot(sess)
        try:
            save_session_snapshot(sess['session_id'], snap)
        except Exception:
            pass
        result['session_snapshot'] = snap
        result['reset_snapshot'] = True

    return jsonify(result)

@app.route('/api/upload', methods=['POST'])
def api_upload():
    """Handle custom expression and metadata file uploads."""
    session_id = request.form.get('session_id')
    sess = get_or_create_session(session_id)

    if 'expression_file' not in request.files or 'metadata_file' not in request.files:
        return jsonify({'success': False, 'errors': ['Both Expression Matrix and Metadata files are required.']})

    expr_file = request.files['expression_file']
    meta_file = request.files['metadata_file']
    input_type = request.form.get('input_type', 'auto')
    missing_strategy = request.form.get('missing_strategy', 'impute_median')

    if expr_file.filename == '' or meta_file.filename == '':
        return jsonify({'success': False, 'errors': ['Please select both files before uploading.']})

    expr_fname = secure_filename(f"expr_{sess['session_id']}_{expr_file.filename}")
    meta_fname = secure_filename(f"meta_{sess['session_id']}_{meta_file.filename}")

    expr_path = os.path.join(sess['dir'], expr_fname)
    meta_path = os.path.join(sess['dir'], meta_fname)

    expr_file.save(expr_path)
    meta_file.save(meta_path)

    result = validate_and_preprocess(expr_path, meta_path, input_type=input_type, handle_missing=missing_strategy, output_dir=sess['dir'])
    if result['success']:
        _reset_downstream_session(sess)
        sess['dataset_name'] = expr_file.filename
        sess['preproc_info'] = result['info']
        sess['expr_preview'] = result.get('expr_preview')
        sess['meta_preview'] = result.get('meta_preview')
        sess['pca_preview'] = result.get('pca_preview')
        sess['processed_expr_path'] = result['processed_expr_path']
        sess['processed_meta_path'] = result['processed_meta_path']
        result['session_id'] = sess['session_id']
        result['dataset_name'] = expr_file.filename
        snap = build_session_snapshot(sess)
        try:
            save_session_snapshot(sess['session_id'], snap)
        except Exception:
            pass
        result['session_snapshot'] = snap
        result['reset_snapshot'] = True

    return jsonify(result)

# ----------------- MODULE 2: DE & PATHWAY ENRICHMENT APIS -----------------

@app.route('/api/run-de', methods=['POST'])
def api_run_de():
    """Execute Differential Expression (Limma / DESeq2 / Python Welch's t-test)."""
    data = request.get_json() or {}
    session_id = data.get('session_id')
    if data.get('session_snapshot') and session_id:
        sess = restore_session_from_snapshot(session_id, data['session_snapshot'])
    else:
        sess = get_or_create_session(session_id)

    if not sess.get('processed_expr_path') or not os.path.exists(sess['processed_expr_path']):
        return jsonify({'success': False, 'needs_restore': True, 'message': 'No validated expression data found. Please upload or load a dataset first.'})

    case_grp = data.get('case_group')
    ctrl_grp = data.get('control_group')
    fdr_thresh = float(data.get('fdr_threshold', 0.05))
    fc_thresh = float(data.get('log2fc_threshold', 1.0))
    method = data.get('method', 'limma')

    inferred_type = sess.get('preproc_info', {}).get('inferred_type', '')
    if inferred_type == 'rnaseq_counts' and method == 'limma':
        method = 'deseq2'

    try:
        results_df, volcano_json, summary = run_differential_expression(
            expr_csv=sess['processed_expr_path'],
            meta_csv=sess['processed_meta_path'],
            case_group=case_grp,
            control_group=ctrl_grp,
            fdr_threshold=fdr_thresh,
            log2fc_threshold=fc_thresh,
            method=method,
            output_dir=sess['dir']
        )

        summary['case_group'] = case_grp
        summary['control_group'] = ctrl_grp

        sess['de_summary'] = summary
        sess['de_df'] = results_df
        sess['volcano_plot'] = volcano_json

        hist_rec = sync_analysis_history(sess)
        snap = build_session_snapshot(sess)

        return jsonify({
            'success': True,
            'summary': summary,
            'volcano_plot': volcano_json,
            'top_genes': summary['top_de_genes'],
            'history_record': hist_rec,
            'session_snapshot': snap
        })
    except Exception as e:
        return jsonify({'success': False, 'message': f"Differential expression error: {str(e)}"})

@app.route('/api/run-enrichment', methods=['POST'])
def api_run_enrichment():
    """Run GSEApy biological pathway enrichment."""
    data = request.get_json() or {}
    session_id = data.get('session_id')
    if data.get('session_snapshot') and session_id:
        sess = restore_session_from_snapshot(session_id, data['session_snapshot'])
    else:
        sess = get_or_create_session(session_id)

    gene_set = data.get('gene_set', 'KEGG_2021_Human')
    sig_genes = (sess.get('de_summary') or {}).get('sig_genes', [])

    if not sig_genes:
        if sess.get('de_df') is not None and not sess['de_df'].empty:
            sig_genes = sess['de_df'].head(25)['Gene'].tolist()
        else:
            return jsonify({'success': False, 'needs_restore': True, 'message': 'No differentially expressed genes available for enrichment.'})

    res = run_pathway_enrichment(sig_genes, gene_set=gene_set)
    if res.get('success'):
        sess['pathway_results'] = res.get('results', [])
        sess['pathway_bar_plot'] = res.get('bar_plot')
        res['history_record'] = sync_analysis_history(sess)
        res['session_snapshot'] = build_session_snapshot(sess)
    return jsonify(res)

# ----------------- MODULE 3: NETWORK & HUB GENES APIS -----------------

@app.route('/api/run-network', methods=['POST'])
def api_run_network():
    """Construct PPI network via STRING DB API & NetworkX, rank hub genes."""
    data = request.get_json() or {}
    session_id = data.get('session_id')
    if data.get('session_snapshot') and session_id:
        sess = restore_session_from_snapshot(session_id, data['session_snapshot'])
    else:
        sess = get_or_create_session(session_id)

    species = int(data.get('species', 9606))
    conf_score = int(data.get('confidence_score', 400))

    sig_genes = (sess.get('de_summary') or {}).get('sig_genes', [])
    if not sig_genes and sess.get('de_df') is not None and not sess['de_df'].empty:
        sig_genes = sess['de_df'].head(30)['Gene'].tolist()

    if not sig_genes:
        return jsonify({'success': False, 'needs_restore': True, 'message': 'No candidate genes found from differential expression.'})

    res = build_and_analyze_network(
        sig_genes=sig_genes,
        de_df=sess.get('de_df'),
        species=species,
        confidence_score=conf_score,
        output_dir=sess['dir']
    )

    if res.get('success'):
        sess['network_summary'] = res['summary']
        sess['hub_genes'] = res['hub_genes']
        sess['network_edges'] = res.get('edges', [])
        sess['network_json'] = res['network_json']
        res['history_record'] = sync_analysis_history(sess)
        res['session_snapshot'] = build_session_snapshot(sess)

    return jsonify(res)

# ----------------- MODULE 4: VALIDATION & REPORT APIS -----------------

@app.route('/api/run-validation', methods=['POST'])
def api_run_validation():
    """Run scikit-learn cross-validation with strict data-leakage prevention."""
    data = request.get_json() or {}
    session_id = data.get('session_id')
    if data.get('session_snapshot') and session_id:
        sess = restore_session_from_snapshot(session_id, data['session_snapshot'])
    else:
        sess = get_or_create_session(session_id)

    if not sess.get('processed_expr_path') or not os.path.exists(sess['processed_expr_path']):
        return jsonify({'success': False, 'needs_restore': True, 'message': 'Dataset not preprocessed.'})

    selected_genes = data.get('candidate_genes') or data.get('selected_genes') or []
    if not selected_genes:
        if sess.get('hub_genes'):
            selected_genes = [h['Gene'] for h in sess['hub_genes'][:10]]
        elif sess.get('de_summary', {}).get('sig_genes'):
            selected_genes = sess['de_summary']['sig_genes'][:10]

    clf_type = data.get('classifier', 'logistic_regression')
    n_splits = int(data.get('cv_splits', 5))

    indep_file_path = os.path.join(Config.SAMPLE_DATA_FOLDER, 'independent_validation_data.csv')
    use_indep = data.get('use_independent_cohort', False)
    indep_csv = indep_file_path if use_indep and os.path.exists(indep_file_path) else None

    res = run_biomarker_classification(
        expr_csv=sess['processed_expr_path'],
        meta_csv=sess['processed_meta_path'],
        candidate_genes=selected_genes,
        classifier_type=clf_type,
        n_splits=n_splits,
        independent_csv=indep_csv,
        output_dir=sess['dir']
    )

    if res.get('success'):
        sess['ml_metrics'] = res
        res['history_record'] = sync_analysis_history(sess)
        res['session_snapshot'] = build_session_snapshot(sess)

    return jsonify(res)

# ----------------- MODULE 5: AI-ASSISTED INTERPRETATION API -----------------

@app.route('/api/ai/chat', methods=['POST'])
def api_ai_chat():
    """Grounded AI-assisted scientific interpretation chat endpoint."""
    data = request.get_json() or {}
    session_id = data.get('session_id')
    user_message = data.get('message', '').strip()
    if data.get('session_snapshot') and session_id:
        sess = restore_session_from_snapshot(session_id, data['session_snapshot'])
    else:
        sess = get_or_create_session(session_id)

    if not user_message:
        return jsonify({'success': False, 'message': 'Empty message received.'})

    res = ask_ai_assistant(user_message, sess)
    
    # Save conversation turn in database
    try:
        save_chat_message(sess['session_id'], 'user', user_message)
        save_chat_message(sess['session_id'], 'ai', res.get('response', ''))
    except Exception:
        pass

    return jsonify(res)

# ----------------- MODULE 6: REPORTS & EXPORTS APIS -----------------

def _ensure_report_file(sess):
    """Build BioPrism_Capstone_Report.html for the given session and return path & HTML text."""
    report_path = os.path.join(sess['dir'], "BioPrism_Capstone_Report.html")
    generate_html_report(
        run_metadata={'run_id': sess['session_id'], 'dataset': sess.get('dataset_name')},
        preproc_info=sess.get('preproc_info') or {},
        de_summary=sess.get('de_summary') or {},
        pathway_results=sess.get('pathway_results') or [],
        network_summary=sess.get('network_summary') or {},
        hub_genes=sess.get('hub_genes') or [],
        ml_metrics=sess.get('ml_metrics') or {},
        output_filepath=report_path
    )
    with open(report_path, 'r', encoding='utf-8') as f:
        return report_path, f.read()

@app.route('/api/generate-report', methods=['POST'])
def api_generate_report():
    """Generate standalone publication-ready HTML capstone report and return HTML + URLs."""
    data = request.get_json() or {}
    session_id = data.get('session_id')
    if data.get('session_snapshot') and session_id:
        sess = restore_session_from_snapshot(session_id, data['session_snapshot'])
    else:
        sess = get_or_create_session(session_id)

    _, html_text = _ensure_report_file(sess)
    return jsonify({
        'success': True,
        'report_html': html_text,
        'preview_url': f"/api/preview/html_report?session_id={sess['session_id']}",
        'report_url': f"/api/download/html_report?session_id={sess['session_id']}"
    })

@app.route('/api/preview/html_report')
def api_preview_html_report():
    """Serve standalone HTML report inline for browser preview/printing (never auto-downloads)."""
    session_id = request.args.get('session_id')
    sess = get_or_create_session(session_id)
    _, html_text = _ensure_report_file(sess)
    return Response(html_text, mimetype='text/html; charset=utf-8')

@app.route('/api/download/<file_type>', methods=['GET', 'POST'])
def api_download(file_type):
    """Serve downloadable files (CSVs, HTML reports, templates). Supports snapshot restore on POST."""
    if request.method == 'POST':
        data = request.get_json(silent=True) or {}
        session_id = data.get('session_id') or request.args.get('session_id')
        if data.get('session_snapshot') and session_id:
            sess = restore_session_from_snapshot(session_id, data['session_snapshot'])
        else:
            sess = get_or_create_session(session_id)
    else:
        session_id = request.args.get('session_id')
        sess = get_or_create_session(session_id) if session_id else {}

    sess_dir = sess.get('dir', Config.RESULTS_FOLDER)
    inline_mode = request.args.get('inline') == '1'

    if file_type == 'template_expression':
        return send_file(os.path.join(Config.SAMPLE_DATA_FOLDER, 'template_expression.csv'), as_attachment=True)
    elif file_type == 'template_metadata':
        return send_file(os.path.join(Config.SAMPLE_DATA_FOLDER, 'template_metadata.csv'), as_attachment=True)
    elif file_type == 'demo_expression':
        return send_file(os.path.join(Config.SAMPLE_DATA_FOLDER, 'demo_microarray_expression.csv'), as_attachment=True)
    elif file_type == 'demo_metadata':
        return send_file(os.path.join(Config.SAMPLE_DATA_FOLDER, 'demo_microarray_metadata.csv'), as_attachment=True)
    elif file_type == 'de_results':
        de_file = (sess.get('de_summary') or {}).get('de_file')
        if de_file and os.path.exists(de_file):
            return send_file(de_file, as_attachment=True, download_name="bioprism_differential_expression.csv")
        de_df = sess.get('de_df')
        if de_df is None and (sess.get('de_summary') or {}).get('top_de_genes'):
            de_df = pd.DataFrame(sess['de_summary']['top_de_genes'])
        if de_df is not None and not de_df.empty:
            csv_bytes = io.BytesIO(de_df.to_csv(index=False).encode('utf-8'))
            return send_file(csv_bytes, mimetype='text/csv', as_attachment=True, download_name="bioprism_differential_expression.csv")
    elif file_type == 'hub_genes':
        hub_file = os.path.join(sess_dir, "ranked_hub_genes.csv")
        if os.path.exists(hub_file):
            return send_file(hub_file, as_attachment=True, download_name="bioprism_ranked_hub_genes.csv")
        if sess.get('hub_genes'):
            csv_bytes = io.BytesIO(pd.DataFrame(sess['hub_genes']).to_csv(index=False).encode('utf-8'))
            return send_file(csv_bytes, mimetype='text/csv', as_attachment=True, download_name="bioprism_ranked_hub_genes.csv")
    elif file_type == 'network_edges':
        edges_file = os.path.join(sess_dir, "network_edges.csv")
        if os.path.exists(edges_file):
            return send_file(edges_file, as_attachment=True, download_name="bioprism_network_edges.csv")
        if sess.get('network_edges'):
            csv_bytes = io.BytesIO(pd.DataFrame(sess['network_edges']).to_csv(index=False).encode('utf-8'))
            return send_file(csv_bytes, mimetype='text/csv', as_attachment=True, download_name="bioprism_network_edges.csv")
    elif file_type == 'validation_metrics':
        metrics_file = os.path.join(sess_dir, "validation_metrics.csv")
        if os.path.exists(metrics_file):
            return send_file(metrics_file, as_attachment=True, download_name="bioprism_validation_metrics.csv")
        ml_m = (sess.get('ml_metrics') or {}).get('metrics')
        if ml_m:
            row = {k: v for k, v in ml_m.items() if k in ['classifier', 'roc_auc', 'accuracy', 'precision', 'recall', 'f1_score', 'num_candidate_genes']}
            csv_bytes = io.BytesIO(pd.DataFrame([row]).to_csv(index=False).encode('utf-8'))
            return send_file(csv_bytes, mimetype='text/csv', as_attachment=True, download_name="bioprism_validation_metrics.csv")
    elif file_type == 'html_report':
        if not sess:
            sess = get_or_create_session(session_id)
        report_path, html_text = _ensure_report_file(sess)
        if inline_mode:
            return Response(html_text, mimetype='text/html; charset=utf-8')
        return send_file(report_path, as_attachment=True, download_name="BioPrism_Analysis_Report.html")

    return jsonify({'error': 'Requested file not found. Please run the corresponding analysis step first.'}), 404

# ----------------- HISTORY APIS -----------------

@app.route('/api/history/list', methods=['GET'])
def api_history_list():
    """Return all recorded analysis runs as JSON."""
    return jsonify({'success': True, 'analyses': get_all_analyses()})

@app.route('/api/history/delete', methods=['POST'])
def api_history_delete():
    """Delete a recorded analysis run."""
    data = request.get_json() or {}
    run_id = data.get('id')
    if run_id:
        delete_analysis(run_id)
        return jsonify({'success': True})
    return jsonify({'success': False, 'message': 'Run ID missing.'})

@app.route('/api/history/clear-all', methods=['POST'])
def api_history_clear_all():
    """Clear all analysis records from SQLite."""
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        cur.execute('DELETE FROM analysis_runs')
        cur.execute('DELETE FROM analysis_history')
        cur.execute('DELETE FROM chat_messages')
        cur.execute('DELETE FROM session_snapshots')
        conn.commit()
        conn.close()
        return jsonify({'success': True})
    except Exception as e:
        return jsonify({'success': False, 'message': str(e)})

if __name__ == '__main__':
    print("=" * 65)
    print("🧬 [BioPrism] Integrated Bioinformatics Web Platform")
    print("Differential Expression | STRING PPI | NetworkX Hubs | ML Validation | AI Assistant")
    print(f"Rscript Path: {Config.RSCRIPT_PATH}")
    print(f"Gemini API: {'CONFIGURED' if get_gemini_api_key() else 'OFFLINE DETERMINISTIC INTERPRETER'}")
    print("Serving on http://127.0.0.1:5000")
    print("=" * 65)
    app.run(host='127.0.0.1', port=5000, debug=False)
