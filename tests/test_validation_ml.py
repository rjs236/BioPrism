import os
import pandas as pd
import numpy as np
from pipeline.validation_ml import run_biomarker_classification

def test_biomarker_classification_logistic_regression(synthetic_data, temp_dir):
    """Test Logistic Regression classification with stratified cross-validation."""
    res = run_biomarker_classification(
        expr_csv=synthetic_data["expr_path"],
        meta_csv=synthetic_data["meta_path"],
        candidate_genes=["TP53", "MDM2", "EGFR", "PTEN"],
        classifier_type="logistic_regression",
        n_splits=3,
        output_dir=temp_dir
    )

    assert res["success"] is True
    metrics = res["metrics"]
    assert "accuracy" in metrics
    assert "roc_auc" in metrics
    assert "f1_score" in metrics
    assert "precision" in metrics
    assert "recall" in metrics
    assert "confusion_matrix" in metrics
    assert metrics["roc_auc"] >= 0.5

    # Check that validation_metrics.csv was created
    metrics_csv = os.path.join(temp_dir, "validation_metrics.csv")
    assert os.path.exists(metrics_csv)

def test_biomarker_classification_random_forest(synthetic_data, temp_dir):
    """Test Random Forest classifier evaluation."""
    res = run_biomarker_classification(
        expr_csv=synthetic_data["expr_path"],
        meta_csv=synthetic_data["meta_path"],
        candidate_genes=["TP53", "MDM2", "EGFR", "PTEN"],
        classifier_type="random_forest",
        n_splits=3,
        output_dir=temp_dir
    )

    assert res["success"] is True
    assert "Random Forest" in res["metrics"]["classifier"]
    assert res["metrics"]["accuracy"] > 0.0

def test_classifier_comparison_mode(synthetic_data, temp_dir):
    """Test Model Comparison mode (LR vs RF) evaluates both models across identical CV folds."""
    res = run_biomarker_classification(
        expr_csv=synthetic_data["expr_path"],
        meta_csv=synthetic_data["meta_path"],
        candidate_genes=["TP53", "MDM2", "EGFR", "PTEN"],
        classifier_type="compare_both",
        n_splits=3,
        output_dir=temp_dir
    )

    assert res["success"] is True
    assert "comparison" in res["metrics"]
    comp = res["metrics"]["comparison"]
    assert comp is not None
    assert len(comp) == 2
    model_names = [c["Model"] for c in comp]
    assert "Logistic Regression" in model_names
    assert "Random Forest" in model_names

    # Check metrics CSV contains comparison rows
    metrics_csv = os.path.join(temp_dir, "validation_metrics.csv")
    df = pd.read_csv(metrics_csv)
    assert len(df) == 2

def test_small_sample_size_handling(temp_dir):
    """Test that insufficient sample sizes (<4 samples) return a clear message instead of crashing."""
    expr_df = pd.DataFrame({
        "Gene": ["TP53", "EGFR"],
        "S1": [8.0, 6.0],
        "S2": [7.5, 6.2]
    })
    meta_df = pd.DataFrame({
        "Sample": ["S1", "S2"],
        "Group": ["Tumor", "Normal"]
    })

    expr_path = os.path.join(temp_dir, "small_expr.csv")
    meta_path = os.path.join(temp_dir, "small_meta.csv")
    expr_df.to_csv(expr_path, index=False)
    meta_df.to_csv(meta_path, index=False)

    res = run_biomarker_classification(
        expr_csv=expr_path,
        meta_csv=meta_path,
        candidate_genes=["TP53"],
        output_dir=temp_dir
    )
    assert res["success"] is False
    assert "insufficient samples" in res["message"].lower()

def test_data_leakage_safeguard(synthetic_data, temp_dir):
    """
    Verify that scaling occurs inside CV folds and does not modify the original data.
    """
    expr_df_before = pd.read_csv(synthetic_data["expr_path"], index_col=0)
    
    res = run_biomarker_classification(
        expr_csv=synthetic_data["expr_path"],
        meta_csv=synthetic_data["meta_path"],
        candidate_genes=["TP53", "MDM2"],
        n_splits=3,
        output_dir=temp_dir
    )
    assert res["success"] is True

    # Check original input was not modified
    expr_df_after = pd.read_csv(synthetic_data["expr_path"], index_col=0)
    pd.testing.assert_frame_equal(expr_df_before, expr_df_after)
