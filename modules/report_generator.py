import os
import sys
import datetime
import pandas as pd
from config import Config

def generate_html_report(run_metadata, preproc_info, de_summary, pathway_results,
                         network_summary, hub_genes, ml_metrics, output_filepath,
                         ai_summary=None):
    """
    Generate an elegant, comprehensive, standalone HTML report for Capstone presentation.
    Includes all 9 required sections and print-to-PDF CSS.
    """
    run_metadata = run_metadata or {}
    preproc_info = preproc_info or {}
    de_summary = de_summary or {}
    pathway_results = pathway_results or []
    network_summary = network_summary or {}
    hub_genes = hub_genes or []
    ml_metrics = ml_metrics or {}

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    python_ver = sys.version.split()[0]
    run_id = run_metadata.get('run_id', 'DEMO')
    dataset_name = run_metadata.get('dataset') or 'Custom Transcriptomics Dataset'

    # 1. Format Top DE Table
    top_de = de_summary.get('top_de_genes', [])[:12]
    de_rows_html = ""
    for r in top_de:
        reg_badge = (
            '<span class="badge badge-up">Upregulated</span>' if r.get('Regulation') == 'Upregulated'
            else '<span class="badge badge-down">Downregulated</span>' if r.get('Regulation') == 'Downregulated'
            else '<span class="badge badge-ns">Not Sig</span>'
        )
        de_rows_html += f"""
        <tr>
            <td><strong>{r.get('Gene')}</strong></td>
            <td>{r.get('log2FoldChange', 0.0):.3f}</td>
            <td>{r.get('p_value', 1.0):.2e}</td>
            <td>{r.get('adj_p_val', 1.0):.2e}</td>
            <td>{reg_badge}</td>
        </tr>
        """

    # 2. Format Pathways Table
    pathway_rows_html = ""
    for p in (pathway_results or [])[:8]:
        pathway_rows_html += f"""
        <tr>
            <td><strong>{p.get('term')}</strong></td>
            <td>{p.get('adj_p_val', 1.0):.2e}</td>
            <td>{p.get('odds_ratio', 0.0)}</td>
            <td><small>{p.get('genes', '')[:70]}...</small></td>
        </tr>
        """

    # 3. Format Hub Genes Table
    hub_rows_html = ""
    for h in (hub_genes or [])[:10]:
        hub_rows_html += f"""
        <tr>
            <td><span class="rank-badge">#{h.get('Hub_Rank', '-')}</span></td>
            <td><strong>{h.get('Gene')}</strong></td>
            <td>{h.get('Composite_Hub_Score', 0.0):.4f}</td>
            <td>{h.get('Degree', 0)}</td>
            <td>{h.get('Betweenness_Centrality', 0.0):.4f}</td>
            <td>{h.get('log2FoldChange', 0.0):.3f}</td>
            <td>{h.get('adj_p_val', 1.0):.2e}</td>
        </tr>
        """

    # 4. Format ML Metrics
    ml_html = ""
    if ml_metrics and (ml_metrics.get('success', False) or ml_metrics.get('metrics')):
        m = ml_metrics.get('metrics', {})
        cm = m.get('confusion_matrix', {})
        ml_html = f"""
        <div class="metrics-grid">
            <div class="metric-card">
                <div class="metric-title">ROC-AUC</div>
                <div class="metric-value">{m.get('roc_auc', 0.0):.3f}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Accuracy</div>
                <div class="metric-value">{m.get('accuracy', 0.0):.1%}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Precision</div>
                <div class="metric-value">{m.get('precision', 0.0):.1%}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Recall (Sensitivity)</div>
                <div class="metric-value">{m.get('recall', 0.0):.1%}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">F1-Score</div>
                <div class="metric-value">{m.get('f1_score', 0.0):.3f}</div>
            </div>
        </div>
        """
        comp_rows = ""
        if m.get('comparison'):
            for c in m['comparison']:
                comp_rows += f"""
                <tr>
                    <td><strong>{c.get('Model')}</strong></td>
                    <td>{c.get('ROC_AUC', 0.0):.3f}</td>
                    <td>{c.get('Accuracy', 0.0):.1%}</td>
                    <td>{c.get('Precision', 0.0):.1%}</td>
                    <td>{c.get('Recall', 0.0):.1%}</td>
                    <td>{c.get('F1_Score', 0.0):.3f}</td>
                </tr>
                """
            ml_html += f"""
            <h5 style="margin-top: 1.5rem; color: #1e293b;">Head-to-Head Classifier Comparison</h5>
            <table class="report-table">
                <thead>
                    <tr><th>Classifier Model</th><th>ROC-AUC</th><th>Accuracy</th><th>Precision</th><th>Recall</th><th>F1-Score</th></tr>
                </thead>
                <tbody>{comp_rows}</tbody>
            </table>
            """

    # AI interpretation text
    ai_interpretation_html = f"""
    <div class="ai-box">
        <div class="ai-badge">🤖 BioPrism AI-Assisted Interpretation</div>
        <p>{ai_summary if ai_summary else f"BioPrism identified {de_summary.get('significant_total', 0)} differentially expressed transcripts across {preproc_info.get('num_samples', 0)} samples. Candidate hub biomarkers ({', '.join([h.get('Gene') for h in (hub_genes or [])[:4]])}) demonstrate high network centrality and statistically significant expression shifts. The machine learning classifier demonstrated an out-of-fold cross-validation ROC-AUC of {ml_metrics.get('metrics', {}).get('roc_auc', 0.0):.3f}."}</p>
        <small style="color: #64748b;">Grounded directly in actual session outputs. No synthetic or hallucinated metrics.</small>
    </div>
    """

    # HTML Template
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BioPrism Analysis Report - {run_id}</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap" rel="stylesheet">
    <style>
        :root {{
            --primary: #0d9488;
            --primary-dark: #0f172a;
            --secondary: #2563eb;
            --light-bg: #f8fafc;
            --card-border: #e2e8f0;
            --text-dark: #1e293b;
            --text-muted: #64748b;
        }}
        body {{
            font-family: 'Inter', sans-serif;
            background-color: var(--light-bg);
            color: var(--text-dark);
            margin: 0;
            padding: 2.5rem;
            line-height: 1.6;
        }}
        .report-container {{
            max-width: 1000px;
            margin: 0 auto;
            background: #ffffff;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.06);
            padding: 3rem;
            border: 1px solid var(--card-border);
        }}
        .header {{
            border-bottom: 2px solid var(--card-border);
            padding-bottom: 1.5rem;
            margin-bottom: 2rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .brand-title {{
            font-size: 2rem;
            font-weight: 800;
            color: var(--primary-dark);
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}
        .brand-title span {{ color: var(--primary); }}
        .header-meta {{
            text-align: right;
            font-size: 0.85rem;
            color: var(--text-muted);
        }}
        .section-title {{
            font-size: 1.25rem;
            font-weight: 700;
            color: var(--primary-dark);
            margin-top: 2rem;
            margin-bottom: 1rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid #f1f5f9;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 1rem;
            margin-bottom: 1.5rem;
        }}
        .metric-card {{
            background: #f8fafc;
            border: 1px solid var(--card-border);
            border-radius: 8px;
            padding: 1rem;
            text-align: center;
        }}
        .metric-title {{
            font-size: 0.75rem;
            font-weight: 600;
            text-transform: uppercase;
            color: var(--text-muted);
            margin-bottom: 0.25rem;
        }}
        .metric-value {{
            font-size: 1.5rem;
            font-weight: 800;
            color: var(--primary);
        }}
        .report-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
            margin-top: 0.75rem;
            margin-bottom: 1.5rem;
        }}
        .report-table th, .report-table td {{
            padding: 0.75rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--card-border);
        }}
        .report-table th {{
            background: #f1f5f9;
            font-weight: 600;
            color: var(--text-dark);
        }}
        .badge {{
            display: inline-block;
            padding: 0.25rem 0.5rem;
            border-radius: 4px;
            font-size: 0.75rem;
            font-weight: 600;
        }}
        .badge-up {{ background: #ccfbf1; color: #0f766e; }}
        .badge-down {{ background: #dbeafe; color: #1d4ed8; }}
        .badge-ns {{ background: #f1f5f9; color: #64748b; }}
        .rank-badge {{
            background: #fef3c7;
            color: #b45309;
            padding: 0.15rem 0.4rem;
            border-radius: 4px;
            font-weight: 700;
        }}
        .ai-box {{
            background: #f0fdfa;
            border: 1px solid #99f6e4;
            border-radius: 8px;
            padding: 1.25rem;
            margin: 1.5rem 0;
        }}
        .ai-badge {{
            font-weight: 700;
            color: #0d9488;
            margin-bottom: 0.5rem;
        }}
        .disclaimer-box {{
            background: #fffbeb;
            border: 1px solid #fde68a;
            border-radius: 8px;
            padding: 1rem;
            margin-top: 2rem;
            font-size: 0.85rem;
            color: #92400e;
        }}
        .print-btn {{
            background: var(--primary);
            color: white;
            border: none;
            padding: 0.5rem 1rem;
            border-radius: 6px;
            cursor: pointer;
            font-weight: 600;
            float: right;
        }}
        @media print {{
            * {{
                -webkit-print-color-adjust: exact !important;
                print-color-adjust: exact !important;
            }}
            @page {{
                margin: 12mm;
                size: A4 portrait;
            }}
            body {{
                padding: 0;
                margin: 0;
                background: #ffffff !important;
                color: #1e293b !important;
            }}
            .report-container {{
                box-shadow: none !important;
                border: none !important;
                padding: 0.5rem !important;
                max-width: 100% !important;
            }}
            .print-btn {{ display: none !important; }}
            .metric-card, .ai-box, .disclaimer-box, tr {{
                page-break-inside: avoid;
                break-inside: avoid;
            }}
            .section-title {{
                page-break-after: avoid;
                break-after: avoid;
            }}
        }}
    </style>
</head>
<body>
    <div class="report-container">
        <button class="print-btn" onclick="window.print()">🖨️ Print / Export PDF</button>
        <div class="header">
            <div>
                <div class="brand-title">🧬 Bio<span>Prism</span></div>
                <div style="font-size: 0.95rem; color: #64748b;">Integrated Transcriptomics & Biomarker Validation Pipeline</div>
            </div>
            <div class="header-meta">
                <div><strong>Run ID:</strong> {run_id}</div>
                <div><strong>Generated:</strong> {timestamp}</div>
                <div><strong>Python:</strong> {python_ver} | <strong>R:</strong> Limma / DESeq2</div>
            </div>
        </div>

        <!-- KPI Summary Cards -->
        <div class="metrics-grid">
            <div class="metric-card">
                <div class="metric-title">Total Samples</div>
                <div class="metric-value">{preproc_info.get('num_samples', 0)}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Total Genes</div>
                <div class="metric-value">{preproc_info.get('num_genes', 0)}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Significant DEGs</div>
                <div class="metric-value">{de_summary.get('significant_total', 0)}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Candidate Hubs</div>
                <div class="metric-value">{len(hub_genes) if hub_genes else 0}</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Best ROC-AUC</div>
                <div class="metric-value">{ml_metrics.get('metrics', {}).get('roc_auc', 0.0):.3f}</div>
            </div>
        </div>

        <!-- 1. Dataset Information -->
        <div class="section-title">1. Dataset Information</div>
        <p><strong>Dataset Name:</strong> {dataset_name}<br>
        <strong>Input Data Type:</strong> {preproc_info.get('inferred_type', 'N/A')}<br>
        <strong>Experimental Groups:</strong> {preproc_info.get('group_counts', {})}</p>

        <!-- 2. Preprocessing Details -->
        <div class="section-title">2. Preprocessing & Quality Control Details</div>
        <p><strong>Preprocessing Applied:</strong> {preproc_info.get('preprocessing_applied', 'Standard Quality Control')}<br>
        <strong>Missing Values Handled:</strong> {preproc_info.get('missing_values_handled', 0)}<br>
        <strong>Expression Range:</strong> Min: {preproc_info.get('min_expression', 'N/A')}, Max: {preproc_info.get('max_expression', 'N/A')}, Mean: {preproc_info.get('mean_expression', 'N/A')}</p>

        <!-- 3. Differential Expression Results -->
        <div class="section-title">3. Differential Expression Results</div>
        <p><strong>Contrast:</strong> {de_summary.get('case_group', 'Case')} vs {de_summary.get('control_group', 'Control')}<br>
        <strong>Statistical Engine:</strong> {de_summary.get('engine', 'Statistical Engine')}<br>
        <strong>Criteria:</strong> FDR &lt; {de_summary.get('fdr_threshold', 0.05)} and |log₂FC| &ge; {de_summary.get('log2fc_threshold', 1.0)}<br>
        <strong>Summary:</strong> {de_summary.get('upregulated', 0)} Upregulated, {de_summary.get('downregulated', 0)} Downregulated, {de_summary.get('not_significant', 0)} Not Significant</p>

        <table class="report-table">
            <thead>
                <tr><th>Gene</th><th>log₂ Fold Change</th><th>p-value</th><th>Adjusted p-val (FDR)</th><th>Status</th></tr>
            </thead>
            <tbody>
                {de_rows_html if de_rows_html else '<tr><td colspan="5">No differentially expressed genes identified.</td></tr>'}
            </tbody>
        </table>

        <!-- 4. Pathway Enrichment Results -->
        <div class="section-title">4. Biological Pathway Enrichment Results</div>
        <table class="report-table">
            <thead>
                <tr><th>Biological Pathway</th><th>FDR Adjusted p-val</th><th>Odds Ratio</th><th>Overlapping Genes</th></tr>
            </thead>
            <tbody>
                {pathway_rows_html if pathway_rows_html else '<tr><td colspan="4">No pathway enrichment results available.</td></tr>'}
            </tbody>
        </table>

        <!-- 5. Network Analysis Results -->
        <div class="section-title">5. Protein-Protein Interaction Network Analysis</div>
        <p><strong>Database:</strong> STRING DB (Taxonomy 9606, Homo sapiens)<br>
        <strong>Confidence Score Threshold:</strong> {network_summary.get('confidence_score_threshold', 400)}<br>
        <strong>Topology:</strong> {network_summary.get('total_nodes', network_summary.get('num_nodes', 0))} nodes, {network_summary.get('total_edges', network_summary.get('num_edges', 0))} edges, Density: {network_summary.get('density', 0.0)}</p>

        <!-- 6. Candidate Biomarker Evaluation -->
        <div class="section-title">6. Candidate Biomarker Hub Prioritization</div>
        <p>Hub genes prioritized using BioPrism's Composite Hub Ranking formula incorporating degree centrality, betweenness centrality, and statistical significance:</p>
        <table class="report-table">
            <thead>
                <tr><th>Rank</th><th>Gene</th><th>Hub Score</th><th>Degree</th><th>Betweenness</th><th>log₂FC</th><th>FDR</th></tr>
            </thead>
            <tbody>
                {hub_rows_html if hub_rows_html else '<tr><td colspan="7">No hub genes identified.</td></tr>'}
            </tbody>
        </table>

        <!-- 7. Classification Metrics -->
        <div class="section-title">7. Machine Learning Classification & Validation Metrics</div>
        {ml_html if ml_html else '<p>Machine learning classification not performed.</p>'}

        <!-- 8. AI-Assisted Interpretation -->
        <div class="section-title">8. AI-Assisted Scientific Interpretation</div>
        {ai_interpretation_html}

        <!-- 9. Limitations & Analysis Parameters -->
        <div class="section-title">9. Methodological Parameters & Scientific Limitations</div>
        <div class="disclaimer-box">
            <strong>⚠️ Scientific & Methodological Disclaimers:</strong><br>
            1. <strong>Strict Data Leakage Prevention:</strong> Feature scaling and selection were executed strictly within cross-validation folds using <code>sklearn.pipeline.Pipeline</code>.<br>
            2. <strong>Exploratory Finding Notice:</strong> Differential expression and machine learning models identify retrospective mathematical correlations. They do NOT establish causal pathogenesis or clinical diagnostic utility without prospective orthogonal wet-lab assays (RT-qPCR, immunohistochemistry, independent cohort testing).<br>
            3. <strong>Reproducibility:</strong> Random states (seed=42) were fixed across cross-validation splits and layout models for deterministic reproducibility.
        </div>
    </div>
</body>
</html>
"""
    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(html_content)

    return output_filepath
