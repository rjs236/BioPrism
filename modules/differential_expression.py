import os
import subprocess
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
import plotly.graph_objects as go
from config import Config

def benjamini_hochberg(p_values):
    """Compute Benjamini-Hochberg adjusted p-values (FDR)."""
    p_values = np.asarray(p_values)
    n = len(p_values)
    if n == 0:
        return np.array([])
    
    order = np.argsort(p_values)
    sorted_p = p_values[order]
    
    ranks = np.arange(1, n + 1)
    adj_p = sorted_p * n / ranks
    
    # Enforce monotonicity: cummin from right to left
    adj_p = np.minimum.accumulate(adj_p[::-1])[::-1]
    adj_p = np.clip(adj_p, 0.0, 1.0)
    
    result = np.empty(n)
    result[order] = adj_p
    return result

def run_differential_expression(expr_csv, meta_csv, case_group, control_group,
                                fdr_threshold=0.05, log2fc_threshold=1.0,
                                method="limma", output_dir=None):
    """
    Module 2: Differential Expression using Limma (R) or DESeq2 (R) with Python statistical fallback.
    Returns:
        results_df: pd.DataFrame with Gene, log2FoldChange, AveExpr, p_value, adj_p_val, Regulation
        volcano_json: dict of Plotly figure
        summary: dict of summary counts, file paths, top DE genes, heatmap and PCA plots
    """
    if output_dir is None:
        output_dir = os.path.dirname(expr_csv)
    os.makedirs(output_dir, exist_ok=True)
    
    de_out_csv = os.path.join(output_dir, f"de_results_{case_group}_vs_{control_group}.csv")
    rscript_bin = Config.RSCRIPT_PATH
    r_used = False
    engine_name = "Limma (R Linear Model & Empirical Bayes)"
    notes = []

    if method == "deseq2":
        r_script = os.path.join(Config.R_SCRIPTS_FOLDER, "run_deseq2.R")
        engine_name = "DESeq2 (R Negative Binomial GLM for Counts)"
    else:
        r_script = os.path.join(Config.R_SCRIPTS_FOLDER, "run_limma.R")
        engine_name = "Limma (R Linear Model & Empirical Bayes)"

    # Try executing R script if Rscript is available
    if rscript_bin and os.path.exists(rscript_bin) and os.path.exists(r_script):
        cmd = [rscript_bin, r_script, expr_csv, meta_csv, de_out_csv, case_group, control_group]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if res.returncode == 0 and os.path.exists(de_out_csv):
                r_used = True
                notes.append(f"Successfully executed native {engine_name} via Rscript.")
            else:
                err_msg = res.stderr or res.stdout
                notes.append(f"R script returned code {res.returncode}. Python statistical fallback engaged.")
        except Exception as e:
            notes.append(f"R execution error: {str(e)[:100]}. Fallback engaged.")

    # Fallback to Python Welch's t-test + BH-FDR if R is not available or failed
    expr_df = pd.read_csv(expr_csv, index_col=0)
    meta_df = pd.read_csv(meta_csv)
    meta_df.columns = [c.strip() for c in meta_df.columns]
    
    sample_col = meta_df.columns[0]
    group_col = meta_df.columns[1]
    
    case_samples = meta_df[meta_df[group_col] == case_group][sample_col].astype(str).tolist()
    ctrl_samples = meta_df[meta_df[group_col] == control_group][sample_col].astype(str).tolist()
    
    case_cols = [c for c in case_samples if c in expr_df.columns]
    ctrl_cols = [c for c in ctrl_samples if c in expr_df.columns]
    
    if len(case_cols) < 2 or len(ctrl_cols) < 2:
        raise ValueError(f"Insufficient samples for comparison: Case ({case_group}) has {len(case_cols)}, Control ({control_group}) has {len(ctrl_cols)}. Minimum 2 replicates required.")

    if not r_used:
        engine_name = "Python Welch's t-test (Two-sided) + Benjamini-Hochberg FDR"
        case_data = expr_df[case_cols].values
        ctrl_data = expr_df[ctrl_cols].values
        
        mean_case = np.nanmean(case_data, axis=1)
        mean_ctrl = np.nanmean(ctrl_data, axis=1)
        ave_expr = np.nanmean(np.hstack([case_data, ctrl_data]), axis=1)
        log2fc = mean_case - mean_ctrl
        
        t_stat, p_vals = stats.ttest_ind(case_data, ctrl_data, axis=1, equal_var=False, nan_policy='omit')
        p_vals = np.nan_to_num(p_vals, nan=1.0)
        p_vals = np.clip(p_vals, 1e-300, 1.0)
        
        adj_p = benjamini_hochberg(p_vals)
        
        de_df = pd.DataFrame({
            'Gene': expr_df.index,
            'log2FoldChange': log2fc,
            'AveExpr': ave_expr,
            't_statistic': t_stat,
            'p_value': p_vals,
            'adj_p_val': adj_p
        })
        de_df = de_df.sort_values(by='p_value').reset_index(drop=True)
        de_df.to_csv(de_out_csv, index=False)
        notes.append("Computed differential expression using Python statistical pipeline.")

    # Read back results to classify and summarize
    results_df = pd.read_csv(de_out_csv)
    results_df['p_value'] = results_df['p_value'].fillna(1.0)
    results_df['adj_p_val'] = results_df['adj_p_val'].fillna(1.0)
    results_df['log2FoldChange'] = results_df['log2FoldChange'].fillna(0.0)

    fdr_thresh = float(fdr_threshold)
    fc_thresh = float(log2fc_threshold)

    conditions = [
        (results_df['adj_p_val'] < fdr_thresh) & (results_df['log2FoldChange'] >= fc_thresh),
        (results_df['adj_p_val'] < fdr_thresh) & (results_df['log2FoldChange'] <= -fc_thresh)
    ]
    choices = ['Upregulated', 'Downregulated']
    results_df['Regulation'] = np.select(conditions, choices, default='Not Significant')

    results_df.to_csv(de_out_csv, index=False)

    counts = results_df['Regulation'].value_counts().to_dict()
    n_up = counts.get('Upregulated', 0)
    n_down = counts.get('Downregulated', 0)
    n_ns = counts.get('Not Significant', 0)
    n_total = len(results_df)

    up_genes = results_df[results_df['Regulation'] == 'Upregulated']['Gene'].tolist()
    down_genes = results_df[results_df['Regulation'] == 'Downregulated']['Gene'].tolist()
    sig_genes = up_genes + down_genes

    # 1. Generate Volcano Plot
    fig_volcano = go.Figure()
    plot_df = results_df.copy()
    plot_df['neg_log10_fdr'] = -np.log10(np.clip(plot_df['adj_p_val'], 1e-50, 1.0))

    # Not Significant
    ns_subset = plot_df[plot_df['Regulation'] == 'Not Significant']
    fig_volcano.add_trace(go.Scatter(
        x=ns_subset['log2FoldChange'],
        y=ns_subset['neg_log10_fdr'],
        mode='markers',
        name=f"Not Sig ({n_ns})",
        marker=dict(color='#94a3b8', size=6, opacity=0.5),
        text=ns_subset['Gene'],
        hovertemplate='<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(FDR): %{y:.2f}<extra></extra>'
    ))

    # Upregulated
    up_subset = plot_df[plot_df['Regulation'] == 'Upregulated']
    fig_volcano.add_trace(go.Scatter(
        x=up_subset['log2FoldChange'],
        y=up_subset['neg_log10_fdr'],
        mode='markers',
        name=f"Upregulated ({n_up})",
        marker=dict(color='#0d9488', size=8, opacity=0.85),
        text=up_subset['Gene'],
        hovertemplate='<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(FDR): %{y:.2f}<br><span style="color:#0d9488">Upregulated</span><extra></extra>'
    ))

    # Downregulated
    down_subset = plot_df[plot_df['Regulation'] == 'Downregulated']
    fig_volcano.add_trace(go.Scatter(
        x=down_subset['log2FoldChange'],
        y=down_subset['neg_log10_fdr'],
        mode='markers',
        name=f"Downregulated ({n_down})",
        marker=dict(color='#2563eb', size=8, opacity=0.85),
        text=down_subset['Gene'],
        hovertemplate='<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(FDR): %{y:.2f}<br><span style="color:#2563eb">Downregulated</span><extra></extra>'
    ))

    # Add threshold cutoffs
    neg_log10_thresh = -np.log10(fdr_thresh)
    fig_volcano.add_hline(y=neg_log10_thresh, line_dash="dash", line_color="#ef4444", annotation_text=f"FDR = {fdr_thresh}", annotation_position="bottom right")
    fig_volcano.add_vline(x=fc_thresh, line_dash="dot", line_color="#64748b", annotation_text=f"+{fc_thresh}")
    fig_volcano.add_vline(x=-fc_thresh, line_dash="dot", line_color="#64748b", annotation_text=f"-{fc_thresh}")

    fig_volcano.update_layout(
        title=f"Volcano Plot: {case_group} vs. {control_group}",
        xaxis_title="log₂ Fold Change",
        yaxis_title="-log₁₀ (Adjusted p-value / FDR)",
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=50, b=40)
    )
    volcano_json = json.loads(fig_volcano.to_json())

    # 2. Generate Expression Heatmap (Top DE genes across all samples)
    heatmap_json = None
    try:
        top_heat_genes = results_df.sort_values(by='adj_p_val').head(25)['Gene'].tolist()
        if top_heat_genes:
            ordered_samples = case_cols + ctrl_cols
            heat_mat = expr_df.loc[top_heat_genes, ordered_samples]
            
            # Row-wise Z-score standard normalization for visual contrast
            z_mat = heat_mat.apply(lambda r: (r - r.mean()) / (r.std() if r.std() > 0 else 1.0), axis=1)

            fig_heat = go.Figure(data=go.Heatmap(
                z=z_mat.values,
                x=ordered_samples,
                y=top_heat_genes,
                colorscale=[[0, '#2563eb'], [0.5, '#f8fafc'], [1, '#0d9488']],
                colorbar=dict(title="Z-Score"),
                hoverongaps=False,
                hovertemplate="<b>Gene:</b> %{y}<br><b>Sample:</b> %{x}<br><b>Z-Score:</b> %{z:.2f}<extra></extra>"
            ))
            fig_heat.update_layout(
                title="Top Differentially Expressed Genes (Standardized Expression Heatmap)",
                xaxis_title="Samples (Sorted by Group)",
                yaxis_title="Top Candidate Genes",
                template="plotly_white",
                margin=dict(l=80, r=40, t=50, b=70)
            )
            heatmap_json = json.loads(fig_heat.to_json())
    except Exception:
        heatmap_json = None

    # 3. Generate PCA Plot (PCA on Samples using DE genes or top variable genes)
    pca_json = None
    try:
        pca_genes = sig_genes[:100] if len(sig_genes) >= 3 else results_df.head(100)['Gene'].tolist()
        sub_expr = expr_df.loc[pca_genes, case_cols + ctrl_cols].T
        if sub_expr.shape[0] >= 3 and sub_expr.shape[1] >= 2:
            pca_model = PCA(n_components=2)
            coords = pca_model.fit_transform(sub_expr)
            var1, var2 = pca_model.explained_variance_ratio_

            fig_pca = go.Figure()
            # Case samples
            case_mask = [s in case_cols for s in sub_expr.index]
            ctrl_mask = [s in ctrl_cols for s in sub_expr.index]

            fig_pca.add_trace(go.Scatter(
                x=coords[case_mask, 0],
                y=coords[case_mask, 1],
                mode='markers+text',
                name=f"{case_group} (n={len(case_cols)})",
                text=[sub_expr.index[i] for i, m in enumerate(case_mask) if m],
                textposition="top center",
                marker=dict(size=12, color='#0d9488', symbol='circle', line=dict(width=1, color='#0f172a')),
                hovertemplate="<b>Sample:</b> %{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>"
            ))

            fig_pca.add_trace(go.Scatter(
                x=coords[ctrl_mask, 0],
                y=coords[ctrl_mask, 1],
                mode='markers+text',
                name=f"{control_group} (n={len(ctrl_cols)})",
                text=[sub_expr.index[i] for i, m in enumerate(ctrl_mask) if m],
                textposition="top center",
                marker=dict(size=12, color='#2563eb', symbol='diamond', line=dict(width=1, color='#0f172a')),
                hovertemplate="<b>Sample:</b> %{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>"
            ))

            fig_pca.update_layout(
                title=f"Sample Separation PCA (PC1: {var1:.1%}, PC2: {var2:.1%})",
                xaxis_title=f"Principal Component 1 ({var1:.1%} Variance)",
                yaxis_title=f"Principal Component 2 ({var2:.1%} Variance)",
                template="plotly_white",
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                margin=dict(l=40, r=40, t=50, b=40)
            )
            pca_json = json.loads(fig_pca.to_json())
    except Exception:
        pca_json = None

    # Summary Output
    summary = {
        'total_genes': n_total,
        'upregulated': n_up,
        'downregulated': n_down,
        'not_significant': n_ns,
        'significant_total': n_up + n_down,
        'fdr_threshold': fdr_thresh,
        'log2fc_threshold': fc_thresh,
        'de_file': de_out_csv,
        'sig_genes': sig_genes,
        'up_genes': up_genes,
        'down_genes': down_genes,
        'engine': engine_name,
        'notes': notes,
        'top_de_genes': results_df.head(20).to_dict(orient='records'),
        'heatmap_plot': heatmap_json,
        'pca_plot': pca_json
    }

    return results_df, volcano_json, summary
