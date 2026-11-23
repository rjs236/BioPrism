import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import json

def validate_and_preprocess(expr_filepath, meta_filepath, input_type="auto", handle_missing="impute_median", log2_transform_if_needed=True):
    """
    Module 1: Comprehensive Input Validation and Preprocessing.
    Returns:
        dict: {
            'success': bool,
            'errors': list of str,
            'warnings': list of str,
            'info': dict (summary stats, group distribution, inferred format),
            'expr_preview': list of dicts,
            'meta_preview': list of dicts,
            'pca_preview': dict (Plotly JSON spec for PCA preview),
            'processed_expr_path': str,
            'processed_meta_path': str
        }
    """
    errors = []
    warnings = []
    
    # 1. Check files existence
    if not os.path.exists(expr_filepath):
        return {'success': False, 'errors': ['Expression matrix file not found.'], 'warnings': []}
    if not os.path.exists(meta_filepath):
        return {'success': False, 'errors': ['Metadata file not found.'], 'warnings': []}

    # 2. Read Metadata
    try:
        sep = '\t' if meta_filepath.endswith('.tsv') or meta_filepath.endswith('.txt') else ','
        meta_df = pd.read_csv(meta_filepath, sep=sep)
    except Exception as e:
        return {'success': False, 'errors': [f'Failed to parse metadata file: {str(e)}'], 'warnings': []}
        
    if meta_df.empty:
        return {'success': False, 'errors': ['Metadata file is empty.'], 'warnings': []}
        
    # Check minimum columns
    if meta_df.shape[1] < 2:
        return {'success': False, 'errors': ['Metadata must contain at least 2 columns: Sample Identifier and Group Label.'], 'warnings': []}
        
    # Standardize metadata column names
    sample_col = meta_df.columns[0]
    group_col = meta_df.columns[1]
    meta_df = meta_df.rename(columns={sample_col: 'Sample', group_col: 'Group'})
    
    # Clean strings
    meta_df['Sample'] = meta_df['Sample'].astype(str).str.strip()
    meta_df['Group'] = meta_df['Group'].astype(str).str.strip()
    
    # Check duplicated samples in metadata
    dup_samples = meta_df['Sample'][meta_df['Sample'].duplicated()].unique().tolist()
    if dup_samples:
        errors.append(f"Duplicated sample identifiers found in metadata: {', '.join(dup_samples[:5])}")

    # Check group labels
    groups = meta_df['Group'].value_counts().to_dict()
    if len(groups) < 2:
        errors.append(f"Metadata must contain at least 2 distinct groups for differential comparison. Found only: {list(groups.keys())}")
    for grp, count in groups.items():
        if count < 2:
            warnings.append(f"Group '{grp}' has only {count} sample(s). Statistical tests require at least 2 replicates per condition.")

    # 3. Read Expression Matrix
    try:
        sep = '\t' if expr_filepath.endswith('.tsv') or expr_filepath.endswith('.txt') else ','
        expr_df = pd.read_csv(expr_filepath, sep=sep)
    except Exception as e:
        return {'success': False, 'errors': [f'Failed to parse expression matrix: {str(e)}'], 'warnings': warnings}

    if expr_df.empty:
        return {'success': False, 'errors': ['Expression matrix is empty.'], 'warnings': warnings}

    # First column must be Gene ID / Symbol
    gene_col = expr_df.columns[0]
    expr_df = expr_df.rename(columns={gene_col: 'Gene'})
    expr_df['Gene'] = expr_df['Gene'].astype(str).str.strip()

    # Check for empty gene names
    if expr_df['Gene'].str.len().eq(0).any() or expr_df['Gene'].isin(['nan', 'None', 'NA']).any():
        warnings.append("Some rows had missing gene identifiers; invalid entries were dropped.")
        expr_df = expr_df[~expr_df['Gene'].isin(['nan', 'None', 'NA', ''])].copy()

    # Check for duplicated genes
    dup_genes = expr_df['Gene'][expr_df['Gene'].duplicated()].unique().tolist()
    if dup_genes:
        warnings.append(f"Detected {len(dup_genes)} duplicated gene symbols (e.g. {', '.join([str(g) for g in dup_genes[:3]])}). Collapsed duplicates by mean expression.")
        expr_df = expr_df.groupby('Gene', as_index=False).mean(numeric_only=True)

    expr_df = expr_df.set_index('Gene')
    
    # 4. Check sample alignment between expression matrix and metadata
    expr_samples = [str(c).strip() for c in expr_df.columns]
    expr_df.columns = expr_samples
    meta_samples = meta_df['Sample'].tolist()

    common_samples = [s for s in meta_samples if s in expr_samples]
    missing_in_expr = [s for s in meta_samples if s not in expr_samples]
    missing_in_meta = [s for s in expr_samples if s not in meta_samples]

    if len(common_samples) < 2:
        errors.append(f"Expression matrix columns and metadata sample IDs have insufficient overlap. Overlapping samples: {len(common_samples)}. Samples in metadata: {len(meta_samples)}, Columns in matrix: {len(expr_samples)}")
        return {'success': False, 'errors': errors, 'warnings': warnings}

    if missing_in_expr:
        warnings.append(f"{len(missing_in_expr)} sample(s) listed in metadata are missing from the expression matrix columns (e.g. {', '.join(missing_in_expr[:3])}).")
    if missing_in_meta:
        warnings.append(f"{len(missing_in_meta)} expression matrix column(s) are not present in metadata and will be excluded (e.g. {', '.join(missing_in_meta[:3])}).")

    # Filter and align columns
    expr_df = expr_df[common_samples]
    meta_df = meta_df[meta_df['Sample'].isin(common_samples)].reset_index(drop=True)

    # 5. Check numeric conversion and missing values
    for col in expr_df.columns:
        if not pd.api.types.is_numeric_dtype(expr_df[col]):
            converted = pd.to_numeric(expr_df[col], errors='coerce')
            n_nan = converted.isna().sum()
            if n_nan > 0:
                warnings.append(f"Column '{col}' had {n_nan} non-numeric value(s) converted to NaN.")
            expr_df[col] = converted

    total_missing = int(expr_df.isna().sum().sum())
    if total_missing > 0:
        if handle_missing == "drop":
            before_len = len(expr_df)
            expr_df = expr_df.dropna()
            warnings.append(f"Dropped {before_len - len(expr_df)} genes with missing values.")
        elif handle_missing == "impute_mean":
            warnings.append(f"Imputed {total_missing} missing values using row-wise mean.")
            expr_df = expr_df.apply(lambda row: row.fillna(row.mean() if not np.isnan(row.mean()) else 0), axis=1)
        else: # impute_median (default)
            warnings.append(f"Imputed {total_missing} missing values using row-wise median.")
            expr_df = expr_df.apply(lambda row: row.fillna(row.median() if not np.isnan(row.median()) else 0), axis=1)

    # Drop zero variance genes
    variances = expr_df.var(axis=1)
    zero_var_genes = variances[variances == 0].index.tolist()
    if zero_var_genes:
        warnings.append(f"Removed {len(zero_var_genes)} gene(s) with zero variance across all samples.")
        expr_df = expr_df[variances > 0]

    if expr_df.empty:
        errors.append("No valid genes remain after quality filtering.")
        return {'success': False, 'errors': errors, 'warnings': warnings}

    # 6. Infer data type (Microarray log2, RNA-seq counts, or RNA-seq normalized)
    min_val = float(expr_df.values.min())
    max_val = float(expr_df.values.max())
    
    # Check if integers (likely RNA-seq counts)
    sample_values = expr_df.values.flatten()
    # sample 1000 values to check
    sample_sub = sample_values[:min(len(sample_values), 1000)]
    is_all_non_negative = (min_val >= 0)
    is_mostly_integer = np.allclose(sample_sub, np.round(sample_sub), atol=1e-3)
    
    inferred_type = "microarray_log2"
    preprocessing_applied = "Standard QC Filtering"

    if input_type == "auto":
        if is_all_non_negative and is_mostly_integer and max_val > 50:
            inferred_type = "rnaseq_counts"
            preprocessing_applied = "Raw Integer Count Matrix Detected (Compatible with DESeq2)"
        elif min_val < 0 or (max_val <= 25 and min_val >= 0):
            inferred_type = "microarray_log2"
            preprocessing_applied = "Log2-transformed expression matrix detected"
        else:
            inferred_type = "normalized_linear"
            if log2_transform_if_needed and min_val >= 0:
                warnings.append("Data appears to be linear non-log expression (Max > 25, non-integer). Applied log2(x + 1) transformation for linear modeling.")
                expr_df = np.log2(expr_df + 1.0)
                inferred_type = "microarray_log2"
                preprocessing_applied = "Log2(x + 1) transformation applied"
    else:
        inferred_type = input_type
        if input_type == "rnaseq_counts" and (min_val < 0 or not is_mostly_integer):
            warnings.append("Selected input type is 'RNA-seq Counts', but data contains negative or non-integer values.")

    # 7. Generate Quality Control PCA Preview
    pca_json = None
    try:
        # Transpose to Samples x Genes
        X_pca = expr_df.T
        if X_pca.shape[0] >= 3 and X_pca.shape[1] >= 2:
            pca = PCA(n_components=2)
            coords = pca.fit_transform(X_pca)
            var_exp = pca.explained_variance_ratio_

            # Map sample groups
            meta_dict = dict(zip(meta_df['Sample'], meta_df['Group']))
            sample_groups = [meta_dict.get(s, 'Unknown') for s in expr_df.columns]

            fig = go.Figure()
            unique_grps = list(dict.fromkeys(sample_groups))
            colors = ['#0d9488', '#2563eb', '#e11d48', '#d97706', '#7c3aed']

            for idx, grp in enumerate(unique_grps):
                mask = [g == grp for g in sample_groups]
                fig.add_trace(go.Scatter(
                    x=coords[mask, 0],
                    y=coords[mask, 1],
                    mode='markers',
                    name=grp,
                    marker=dict(size=11, color=colors[idx % len(colors)], line=dict(width=1, color='#ffffff')),
                    text=[expr_df.columns[i] for i, m in enumerate(mask) if m],
                    hovertemplate='<b>Sample:</b> %{text}<br><b>PC1:</b> %{x:.2f}<br><b>PC2:</b> %{y:.2f}<extra></extra>'
                ))

            fig.update_layout(
                title=f"Sample Clustering PCA Preview (PC1: {var_exp[0]:.1%}, PC2: {var_exp[1]:.1%})",
                xaxis_title=f"PC1 ({var_exp[0]:.1%} variance)",
                yaxis_title=f"PC2 ({var_exp[1]:.1%} variance)",
                template="plotly_white",
                margin=dict(l=40, r=40, t=50, b=40),
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
            )
            pca_json = json.loads(fig.to_json())
    except Exception:
        pca_json = None

    # 8. Save Processed Files
    proc_dir = os.path.dirname(expr_filepath)
    processed_expr_path = os.path.join(proc_dir, "processed_expression.csv")
    processed_meta_path = os.path.join(proc_dir, "processed_metadata.csv")

    expr_df.to_csv(processed_expr_path)
    meta_df.to_csv(processed_meta_path, index=False)

    # 9. Previews
    preview_rows = 8
    preview_cols = min(6, len(expr_df.columns))
    sub_expr = expr_df.iloc[:preview_rows, :preview_cols].reset_index()
    expr_preview = sub_expr.to_dict(orient='records')
    meta_preview = meta_df.head(preview_rows).to_dict(orient='records')

    # Summary Info
    info = {
        'num_genes': int(len(expr_df)),
        'num_samples': int(len(common_samples)),
        'available_groups': list(groups.keys()),
        'group_counts': groups,
        'inferred_type': inferred_type,
        'preprocessing_applied': preprocessing_applied,
        'min_expression': round(float(expr_df.values.min()), 3),
        'max_expression': round(float(expr_df.values.max()), 3),
        'mean_expression': round(float(expr_df.values.mean()), 3),
        'missing_values_handled': total_missing
    }

    if errors:
        return {'success': False, 'errors': errors, 'warnings': warnings}

    return {
        'success': True,
        'errors': [],
        'warnings': warnings,
        'info': info,
        'expr_preview': expr_preview,
        'meta_preview': meta_preview,
        'pca_preview': pca_json,
        'processed_expr_path': processed_expr_path,
        'processed_meta_path': processed_meta_path
    }
