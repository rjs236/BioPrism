import os
import sqlite3
import datetime
import json
from config import Config

def get_db_connection():
    """Get a SQLite connection with row factory enabled."""
    conn = sqlite3.connect(Config.DATABASE_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    """Initialize database tables for BioPrism."""
    conn = get_db_connection()
    cur = conn.cursor()
    
    # 1. Datasets table
    cur.execute('''
        CREATE TABLE IF NOT EXISTS datasets (
            id TEXT PRIMARY KEY,
            name TEXT NOT NULL,
            data_type TEXT,
            num_genes INTEGER DEFAULT 0,
            num_samples INTEGER DEFAULT 0,
            groups_json TEXT,
            expr_path TEXT,
            meta_path TEXT,
            uploaded_at TEXT
        )
    ''')

    # 2. Analysis Runs table
    cur.execute('''
        CREATE TABLE IF NOT EXISTS analysis_runs (
            id TEXT PRIMARY KEY,
            timestamp TEXT NOT NULL,
            dataset_name TEXT,
            input_type TEXT,
            total_genes INTEGER DEFAULT 0,
            total_samples INTEGER DEFAULT 0,
            case_group TEXT,
            control_group TEXT,
            de_method TEXT,
            fdr_threshold REAL,
            log2fc_threshold REAL,
            sig_genes_count INTEGER DEFAULT 0,
            n_up INTEGER DEFAULT 0,
            n_down INTEGER DEFAULT 0,
            top_hub_genes TEXT,
            top_pathways TEXT,
            classifier_type TEXT,
            roc_auc REAL DEFAULT 0.0,
            accuracy REAL DEFAULT 0.0,
            precision_score REAL DEFAULT 0.0,
            recall_score REAL DEFAULT 0.0,
            f1_score REAL DEFAULT 0.0,
            details_json TEXT
        )
    ''')

    # 3. Chat Messages table for AI assistant
    cur.execute('''
        CREATE TABLE IF NOT EXISTS chat_messages (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            session_id TEXT NOT NULL,
            role TEXT NOT NULL,
            message TEXT NOT NULL,
            timestamp TEXT NOT NULL
        )
    ''')

    # 4. Settings table
    cur.execute('''
        CREATE TABLE IF NOT EXISTS settings (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL,
            updated_at TEXT
        )
    ''')

    # Legacy table synchronization for existing tests that expect analysis_history
    cur.execute('''
        CREATE TABLE IF NOT EXISTS analysis_history (
            id TEXT PRIMARY KEY,
            timestamp TEXT,
            dataset_name TEXT,
            input_type TEXT,
            total_genes INTEGER,
            total_samples INTEGER,
            case_group TEXT,
            control_group TEXT,
            fdr_threshold REAL,
            log2fc_threshold REAL,
            sig_genes_count INTEGER,
            top_hub_genes TEXT,
            roc_auc REAL,
            accuracy REAL
        )
    ''')

    conn.commit()
    conn.close()

def save_analysis_run(run_data):
    """
    Save complete analysis run to database.
    run_data: dict with run metadata and metrics
    """
    conn = get_db_connection()
    cur = conn.cursor()
    run_id = run_data.get('id', datetime.datetime.now().strftime("%Y%m%d%H%M%S"))
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    top_hubs_str = run_data.get('top_hub_genes', '')
    if isinstance(top_hubs_str, list):
        top_hubs_str = ", ".join([str(h) for h in top_hubs_str])

    top_paths_str = run_data.get('top_pathways', '')
    if isinstance(top_paths_str, list):
        top_paths_str = ", ".join([str(p) for p in top_paths_str])

    details_str = json.dumps(run_data.get('details', {}))

    cur.execute('''
        INSERT OR REPLACE INTO analysis_runs
        (id, timestamp, dataset_name, input_type, total_genes, total_samples,
         case_group, control_group, de_method, fdr_threshold, log2fc_threshold,
         sig_genes_count, n_up, n_down, top_hub_genes, top_pathways,
         classifier_type, roc_auc, accuracy, precision_score, recall_score, f1_score, details_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        run_id, now_str,
        run_data.get('dataset_name', 'Custom Dataset'),
        run_data.get('input_type', 'Microarray Log2'),
        int(run_data.get('total_genes', 0)),
        int(run_data.get('total_samples', 0)),
        run_data.get('case_group', 'Case'),
        run_data.get('control_group', 'Control'),
        run_data.get('de_method', 'limma'),
        float(run_data.get('fdr_threshold', 0.05)),
        float(run_data.get('log2fc_threshold', 1.0)),
        int(run_data.get('sig_genes_count', 0)),
        int(run_data.get('n_up', 0)),
        int(run_data.get('n_down', 0)),
        top_hubs_str,
        top_paths_str,
        run_data.get('classifier_type', 'Logistic Regression'),
        float(run_data.get('roc_auc', 0.0)),
        float(run_data.get('accuracy', 0.0)),
        float(run_data.get('precision_score', 0.0)),
        float(run_data.get('recall_score', 0.0)),
        float(run_data.get('f1_score', 0.0)),
        details_str
    ))

    # Also mirror into legacy analysis_history for backwards compatibility
    cur.execute('''
        INSERT OR REPLACE INTO analysis_history
        (id, timestamp, dataset_name, input_type, total_genes, total_samples,
         case_group, control_group, fdr_threshold, log2fc_threshold,
         sig_genes_count, top_hub_genes, roc_auc, accuracy)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        run_id, now_str,
        run_data.get('dataset_name', 'Custom Dataset'),
        run_data.get('input_type', 'Microarray Log2'),
        int(run_data.get('total_genes', 0)),
        int(run_data.get('total_samples', 0)),
        run_data.get('case_group', 'Case'),
        run_data.get('control_group', 'Control'),
        float(run_data.get('fdr_threshold', 0.05)),
        float(run_data.get('log2fc_threshold', 1.0)),
        int(run_data.get('sig_genes_count', 0)),
        top_hubs_str,
        float(run_data.get('roc_auc', 0.0)),
        float(run_data.get('accuracy', 0.0))
    ))

    conn.commit()
    conn.close()
    return run_id

def get_dashboard_stats():
    """
    Compute actual dynamic dashboard statistics from database records.
    Never returns hardcoded placeholder values.
    """
    conn = get_db_connection()
    cur = conn.cursor()

    # Total runs
    cur.execute('SELECT COUNT(*) FROM analysis_runs')
    total_runs = cur.fetchone()[0]

    # Total distinct datasets analyzed
    cur.execute('SELECT COUNT(DISTINCT dataset_name) FROM analysis_runs')
    total_datasets = cur.fetchone()[0]

    # Total DE genes detected across runs
    cur.execute('SELECT COALESCE(SUM(sig_genes_count), 0) FROM analysis_runs')
    total_de_genes = cur.fetchone()[0]

    # Total unique candidate hub genes identified
    cur.execute('SELECT top_hub_genes FROM analysis_runs WHERE top_hub_genes IS NOT NULL AND top_hub_genes != ""')
    rows = cur.fetchall()
    all_hubs = set()
    for r in rows:
        genes = [g.strip().upper() for g in r[0].split(',') if g.strip()]
        all_hubs.update(genes)
    total_biomarkers = len(all_hubs)

    # Total unique pathways enriched
    cur.execute('SELECT top_pathways FROM analysis_runs WHERE top_pathways IS NOT NULL AND top_pathways != ""')
    p_rows = cur.fetchall()
    all_pathways = set()
    for pr in p_rows:
        paths = [p.strip() for p in pr[0].split(',') if p.strip()]
        all_pathways.update(paths)
    total_pathways = len(all_pathways)

    # Recent analysis history
    cur.execute('''
        SELECT id, timestamp, dataset_name, input_type, total_genes, total_samples,
               case_group, control_group, sig_genes_count, top_hub_genes, classifier_type, roc_auc, accuracy
        FROM analysis_runs
        ORDER BY timestamp DESC
        LIMIT 6
    ''')
    recent_runs = [dict(row) for row in cur.fetchall()]

    # Input types breakdown
    cur.execute('SELECT input_type, COUNT(*) FROM analysis_runs GROUP BY input_type')
    type_counts = dict(cur.fetchall())

    # Accuracy / AUC distribution
    cur.execute('SELECT id, dataset_name, roc_auc, accuracy, timestamp FROM analysis_runs ORDER BY timestamp ASC LIMIT 10')
    trend_data = [dict(row) for row in cur.fetchall()]

    conn.close()

    return {
        'total_runs': total_runs,
        'total_datasets': total_datasets,
        'total_de_genes': total_de_genes,
        'total_biomarkers': total_biomarkers,
        'total_pathways': total_pathways,
        'recent_runs': recent_runs,
        'type_counts': type_counts,
        'trend_data': trend_data
    }

def get_all_analyses():
    """Retrieve full history of analysis runs."""
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute('''
        SELECT id, timestamp, dataset_name, input_type, total_genes, total_samples,
               case_group, control_group, de_method, fdr_threshold, log2fc_threshold,
               sig_genes_count, n_up, n_down, top_hub_genes, top_pathways,
               classifier_type, roc_auc, accuracy, precision_score, recall_score, f1_score
        FROM analysis_runs
        ORDER BY timestamp DESC
    ''')
    analyses = [dict(row) for row in cur.fetchall()]
    conn.close()
    return analyses

def get_analysis_by_id(run_id):
    """Retrieve single analysis run details by ID."""
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute('SELECT * FROM analysis_runs WHERE id = ?', (run_id,))
    row = cur.fetchone()
    conn.close()
    if row:
        data = dict(row)
        if data.get('details_json'):
            try:
                data['details'] = json.loads(data['details_json'])
            except Exception:
                data['details'] = {}
        return data
    return None

def delete_analysis(run_id):
    """Delete an analysis run from database."""
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute('DELETE FROM analysis_runs WHERE id = ?', (run_id,))
    cur.execute('DELETE FROM analysis_history WHERE id = ?', (run_id,))
    conn.commit()
    conn.close()
    return True

def save_chat_message(session_id, role, message):
    """Save an AI chat turn to database."""
    conn = get_db_connection()
    cur = conn.cursor()
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    cur.execute('''
        INSERT INTO chat_messages (session_id, role, message, timestamp)
        VALUES (?, ?, ?, ?)
    ''', (session_id, role, message, now_str))
    conn.commit()
    conn.close()

def get_chat_history(session_id, limit=20):
    """Get chat conversation history for a given session."""
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute('''
        SELECT role, message, timestamp FROM chat_messages
        WHERE session_id = ?
        ORDER BY id ASC
        LIMIT ?
    ''', (session_id, limit))
    msgs = [dict(row) for row in cur.fetchall()]
    conn.close()
    return msgs

def get_setting(key, default=None):
    """Get application setting."""
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute('SELECT value FROM settings WHERE key = ?', (key,))
    row = cur.fetchone()
    conn.close()
    return row[0] if row else default

def set_setting(key, value):
    """Set application setting."""
    conn = get_db_connection()
    cur = conn.cursor()
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    cur.execute('''
        INSERT OR REPLACE INTO settings (key, value, updated_at)
        VALUES (?, ?, ?)
    ''', (key, str(value), now_str))
    conn.commit()
    conn.close()
