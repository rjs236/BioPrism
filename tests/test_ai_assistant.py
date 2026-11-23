import os
import pytest
from modules.ai_assistant import (
    build_session_context_summary,
    generate_deterministic_scientific_response,
    ask_ai_assistant
)

def test_session_context_building():
    """Verify that build_session_context_summary compiles actual session metrics."""
    mock_session = {
        'session_id': 'TEST_SESS',
        'dataset_name': 'TCGA-BRCA Benchmark',
        'preproc_info': {
            'inferred_type': 'microarray_log2',
            'num_genes': 500,
            'num_samples': 40,
            'group_counts': {'Tumor': 20, 'Normal': 20}
        },
        'de_summary': {
            'case_group': 'Tumor',
            'control_group': 'Normal',
            'fdr_threshold': 0.05,
            'log2fc_threshold': 1.0,
            'significant_total': 15,
            'upregulated': 10,
            'downregulated': 5,
            'top_de_genes': [
                {'Gene': 'TP53', 'log2FoldChange': 3.2, 'adj_p_val': 1.2e-6, 'Regulation': 'Upregulated'},
                {'Gene': 'PTEN', 'log2FoldChange': -2.5, 'adj_p_val': 3.4e-5, 'Regulation': 'Downregulated'}
            ]
        },
        'pathway_results': [
            {'term': 'p53 signaling pathway', 'adj_p_val': 1.5e-5, 'odds_ratio': 14.2, 'genes': 'TP53, MDM2'}
        ],
        'network_summary': {
            'total_nodes': 15,
            'total_edges': 24,
            'density': 0.22,
            'confidence_score_threshold': 400
        },
        'hub_genes': [
            {'Hub_Rank': 1, 'Gene': 'TP53', 'Composite_Hub_Score': 0.895, 'Degree': 8, 'Betweenness_Centrality': 0.42, 'log2FoldChange': 3.2}
        ],
        'ml_metrics': {
            'success': True,
            'metrics': {
                'model_name': 'Logistic Regression (L2)',
                'roc_auc': 0.945,
                'accuracy': 0.92,
                'precision': 0.90,
                'recall': 0.95,
                'f1_score': 0.925,
                'candidate_genes': ['TP53', 'PTEN'],
                'confusion_matrix': {'tp': 19, 'tn': 18, 'fp': 2, 'fn': 1}
            },
            'cv_method': '5-Fold Stratified Cross-Validation'
        }
    }

    ctx = build_session_context_summary(mock_session)
    assert "TCGA-BRCA Benchmark" in ctx
    assert "TP53" in ctx
    assert "0.945" in ctx
    assert "p53 signaling pathway" in ctx

def test_ai_deterministic_grounded_response():
    """Verify deterministic scientific interpreter answers without hallucination."""
    mock_session = {
        'session_id': 'TEST_SESS',
        'dataset_name': 'Test Cohort',
        'de_summary': {
            'case_group': 'Disease',
            'control_group': 'Control',
            'fdr_threshold': 0.05,
            'log2fc_threshold': 1.0,
            'significant_total': 8,
            'upregulated': 5,
            'downregulated': 3
        },
        'hub_genes': [
            {'Hub_Rank': 1, 'Gene': 'EGFR', 'Composite_Hub_Score': 0.912, 'Degree': 6, 'Betweenness_Centrality': 0.35, 'log2FoldChange': 2.8, 'adj_p_val': 1e-4, 'Regulation': 'Upregulated'}
        ],
        'pathway_results': [
            {'term': 'ErbB signaling pathway', 'adj_p_val': 2.1e-4, 'odds_ratio': 12.0, 'genes': 'EGFR, AKT1'}
        ],
        'ml_metrics': {
            'success': True,
            'metrics': {
                'model_name': 'Logistic Regression',
                'roc_auc': 0.92,
                'accuracy': 0.88,
                'precision': 0.85,
                'recall': 0.90,
                'f1_score': 0.87
            },
            'cv_method': '5-Fold Stratified Cross-Validation'
        }
    }

    # Query about specific gene EGFR
    resp = generate_deterministic_scientific_response("Why was EGFR selected as a biomarker?", mock_session)
    assert "EGFR" in resp
    assert "0.912" in resp or "Rank #1" in resp

    # Query about FDR / log2FC
    resp_stats = generate_deterministic_scientific_response("Explain log2 fold change and FDR threshold", mock_session)
    assert "log₂ Fold Change" in resp_stats
    assert "False Discovery Rate" in resp_stats or "Benjamini–Hochberg" in resp_stats

    # Query about pathways
    resp_path = generate_deterministic_scientific_response("Summarize enriched pathways", mock_session)
    assert "ErbB signaling pathway" in resp_path

    # Query via main ask_ai_assistant endpoint
    full_res = ask_ai_assistant("Evaluate machine learning classification performance", mock_session)
    assert full_res["success"] is True
    assert "0.92" in full_res["response"] or "ROC-AUC" in full_res["response"]
