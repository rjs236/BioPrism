import os
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
import plotly.graph_objects as go
from modules.classification import get_classifier, evaluate_classifier_cv

def run_biomarker_classification(expr_csv, meta_csv, candidate_genes,
                                 classifier_type="logistic_regression",
                                 n_splits=5, independent_csv=None,
                                 output_dir=None):
    """
    Module 4: Machine Learning Biomarker Validation.
    Performs rigorous cross-validation with complete data-leakage prevention
    by strictly encapsulating scaling inside CV folds.
    Supports Logistic Regression, Random Forest, SVM, and Comparison mode.
    Also tests against independent cohort if provided.
    """
    if output_dir is None:
        output_dir = "."
    os.makedirs(output_dir, exist_ok=True)

    # 1. Load Data
    expr_df = pd.read_csv(expr_csv, index_col=0)
    meta_df = pd.read_csv(meta_csv)
    meta_df.columns = [c.strip() for c in meta_df.columns]
    sample_col = meta_df.columns[0]
    group_col = meta_df.columns[1]

    common_samples = [s for s in meta_df[sample_col].astype(str) if s in expr_df.columns]
    if len(common_samples) < 4:
        return {
            'success': False,
            'message': f"Insufficient samples ({len(common_samples)}) for cross-validation (minimum 4 samples required)."
        }

    meta_sub = meta_df[meta_df[sample_col].astype(str).isin(common_samples)].copy()
    groups = meta_sub[group_col].unique().tolist()
    if len(groups) != 2:
        return {
            'success': False,
            'message': f"Binary classification requires exactly 2 groups. Found: {groups}"
        }

    case_grp = groups[0]
    ctrl_grp = groups[1]

    valid_genes = [g for g in candidate_genes if g in expr_df.index]
    if len(valid_genes) == 0:
        return {
            'success': False,
            'message': "None of the selected candidate genes were found in the expression matrix."
        }

    # Transpose to Samples x Genes format
    X_df = expr_df.loc[valid_genes, common_samples].T
    y = np.array([1 if g == case_grp else 0 for g in meta_sub[group_col]])
    X = X_df.values

    n_samples, n_features = X_df.shape
    class_counts = pd.Series(y).value_counts().to_dict()

    min_class_count = min(class_counts.values())
    if min_class_count < 2:
        return {
            'success': False,
            'message': f"Cannot perform cross-validation: smallest class has only {min_class_count} sample(s)."
        }

    actual_splits = min(n_splits, min_class_count)
    if actual_splits >= 3:
        cv = StratifiedKFold(n_splits=actual_splits, shuffle=True, random_state=42)
        cv_method = f"{actual_splits}-Fold Stratified Cross-Validation"
    else:
        cv = LeaveOneOut()
        cv_method = "Leave-One-Out Cross-Validation (LOOCV due to small cohort)"

    models_to_run = []
    if classifier_type == "compare_both":
        # Specific 2-model comparison (LR vs RF)
        est_lr, name_lr, col_lr = get_classifier("logistic_regression")
        est_rf, name_rf, col_rf = get_classifier("random_forest")
        models_to_run = [(est_lr, name_lr, col_lr), (est_rf, name_rf, col_rf)]
    elif classifier_type in ["compare_all", "compare_three"]:
        # 3-model comparison (LR vs RF vs SVM)
        est_lr, name_lr, col_lr = get_classifier("logistic_regression")
        est_rf, name_rf, col_rf = get_classifier("random_forest")
        est_svm, name_svm, col_svm = get_classifier("svm")
        models_to_run = [(est_lr, name_lr, col_lr), (est_rf, name_rf, col_rf), (est_svm, name_svm, col_svm)]
    else:
        est, name, col = get_classifier(classifier_type)
        models_to_run = [(est, name, col)]

    results_list = []
    for estimator, name, col in models_to_run:
        eval_res = evaluate_classifier_cv(X, y, estimator, name, col, cv, valid_genes)
        results_list.append(eval_res)

    # Primary model is the first one or the one with highest AUC
    primary_res = sorted(results_list, key=lambda x: x['roc_auc'], reverse=True)[0] if len(results_list) > 1 else results_list[0]

    # 1. Interactive ROC Plot
    fig_roc = go.Figure()
    for m in results_list:
        fig_roc.add_trace(go.Scatter(
            x=m['fpr'], y=m['tpr'],
            mode='lines',
            name=f"{m['model_name']} (AUC = {m['roc_auc']:.3f})",
            line=dict(color=m['color'], width=2.5)
        ))

    fig_roc.add_trace(go.Scatter(
        x=[0, 1], y=[0, 1],
        mode='lines',
        name='Random Chance (AUC = 0.50)',
        line=dict(color='#94a3b8', dash='dash', width=1.5)
    ))

    fig_roc.update_layout(
        title=f"Receiver Operating Characteristic (ROC) — {cv_method}",
        xaxis_title="False Positive Rate (1 - Specificity)",
        yaxis_title="True Positive Rate (Sensitivity / Recall)",
        template="plotly_white",
        legend=dict(x=0.45, y=0.15, bgcolor='rgba(255,255,255,0.85)'),
        margin=dict(l=40, r=40, t=50, b=40)
    )
    roc_json = json.loads(fig_roc.to_json())

    # 2. Interactive Confusion Matrix Plot
    cm_dict = primary_res['confusion_matrix']
    cm_matrix = [[cm_dict['tn'], cm_dict['fp']], [cm_dict['fn'], cm_dict['tp']]]
    fig_cm = go.Figure(data=go.Heatmap(
        z=cm_matrix,
        x=[f"Pred {ctrl_grp}", f"Pred {case_grp}"],
        y=[f"True {ctrl_grp}", f"True {case_grp}"],
        colorscale=[[0, '#f8fafc'], [1, '#0d9488']],
        text=[[f"TN: {cm_dict['tn']}", f"FP: {cm_dict['fp']}"],
              [f"FN: {cm_dict['fn']}", f"TP: {cm_dict['tp']}"]],
        texttemplate="%{text}",
        textfont={"size": 16, "family": "Inter"},
        showscale=False
    ))
    fig_cm.update_layout(
        title=f"Out-of-Fold Confusion Matrix ({primary_res['model_name']})",
        template="plotly_white",
        margin=dict(l=60, r=40, t=50, b=40)
    )
    cm_json = json.loads(fig_cm.to_json())

    # 3. Independent Validation Cohort Evaluation & Expression Comparison
    indep_results = None
    expr_comparison_plot = None
    if independent_csv and os.path.exists(independent_csv):
        try:
            indep_df = pd.read_csv(independent_csv)
            indep_cols = [c.strip() for c in indep_df.columns]
            indep_df.columns = indep_cols

            indep_genes = [g for g in valid_genes if g in indep_df.columns]
            if len(indep_genes) >= 1 and 'Group' in indep_df.columns:
                # Align labels
                indep_case_mask = indep_df['Group'] == case_grp
                y_indep = np.array([1 if g == case_grp else 0 for g in indep_df['Group']])

                # Impute or fill any missing candidate genes with zeros
                X_indep_df = pd.DataFrame(index=indep_df.index)
                for g in valid_genes:
                    if g in indep_df.columns:
                        X_indep_df[g] = indep_df[g].values
                    else:
                        X_indep_df[g] = 0.0

                X_indep = X_indep_df.values
                trained_pipe = primary_res['pipeline']
                indep_preds = trained_pipe.predict(X_indep)

                if hasattr(trained_pipe.named_steps['classifier'], 'predict_proba'):
                    indep_probs = trained_pipe.predict_proba(X_indep)[:, 1]
                else:
                    indep_probs = trained_pipe.decision_function(X_indep)

                indep_acc = float(accuracy_score(y_indep, indep_preds))
                try:
                    indep_auc = float(roc_auc_score(y_indep, indep_probs))
                except Exception:
                    indep_auc = 0.5

                indep_results = {
                    'dataset_name': os.path.basename(independent_csv),
                    'total_samples': len(indep_df),
                    'matching_genes': len(indep_genes),
                    'accuracy': indep_acc,
                    'roc_auc': indep_auc,
                    'message': f"Evaluated on {len(indep_df)} independent samples ({len(indep_genes)}/{len(valid_genes)} candidate genes matched)."
                }

                # Cross-cohort candidate gene expression comparison plot
                fig_comp = go.Figure()
                top_comp_genes = valid_genes[:4]
                for g in top_comp_genes:
                    if g in expr_df.index:
                        fig_comp.add_trace(go.Box(
                            y=expr_df.loc[g, common_samples].values,
                            name=f"Discovery: {g}",
                            boxpoints='all',
                            jitter=0.3,
                            marker=dict(color='#0d9488', size=4)
                        ))
                    if g in indep_df.columns:
                        fig_comp.add_trace(go.Box(
                            y=indep_df[g].values,
                            name=f"Independent: {g}",
                            boxpoints='all',
                            jitter=0.3,
                            marker=dict(color='#2563eb', size=4)
                        ))

                fig_comp.update_layout(
                    title="Cross-Cohort Candidate Biomarker Expression Distribution",
                    yaxis_title="Expression Level",
                    template="plotly_white",
                    margin=dict(l=40, r=40, t=50, b=50)
                )
                expr_comparison_plot = json.loads(fig_comp.to_json())
        except Exception as e:
            indep_results = {
                'message': f"Failed to evaluate on independent cohort: {str(e)[:80]}"
            }

    # Save metrics table
    metrics_summary = {
        'classifier': primary_res['model_name'],
        'model_name': primary_res['model_name'],
        'cv_method': cv_method,
        'candidate_genes': valid_genes,
        'num_candidate_genes': len(valid_genes),
        'accuracy': primary_res['accuracy'],
        'precision': primary_res['precision'],
        'recall': primary_res['recall'],
        'f1_score': primary_res['f1_score'],
        'roc_auc': primary_res['roc_auc'],
        'confusion_matrix': cm_dict,
        'gene_importances': primary_res['gene_importances'],
        'comparison': [
            {
                'Model': m['model_name'],
                'ROC_AUC': m['roc_auc'],
                'Accuracy': m['accuracy'],
                'Precision': m['precision'],
                'Recall': m['recall'],
                'F1_Score': m['f1_score']
            } for m in results_list
        ]
    }

    # Write metrics to CSV
    metrics_df = pd.DataFrame([metrics_summary['comparison'][0]])
    metrics_csv_path = os.path.join(output_dir, "validation_metrics.csv")
    pd.DataFrame(metrics_summary['comparison']).to_csv(metrics_csv_path, index=False)

    return {
        'success': True,
        'classifier_type': classifier_type,
        'cv_method': cv_method,
        'metrics': metrics_summary,
        'roc_plot': roc_json,
        'confusion_matrix_plot': cm_json,
        'independent_validation': indep_results,
        'expression_comparison_plot': expr_comparison_plot
    }
