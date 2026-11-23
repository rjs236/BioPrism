import os
import pandas as pd
import numpy as np
from pipeline.preprocessor import validate_and_preprocess

def test_valid_preprocessing(synthetic_data):
    """Test standard preprocessing with valid synthetic data."""
    res = validate_and_preprocess(synthetic_data["expr_path"], synthetic_data["meta_path"])
    assert res["success"] is True
    assert len(res["errors"]) == 0
    assert res["info"]["num_genes"] == 25
    assert res["info"]["num_samples"] == 10
    assert os.path.exists(res["processed_expr_path"])
    assert os.path.exists(res["processed_meta_path"])
    assert "Tumor" in res["info"]["available_groups"]
    assert "Normal" in res["info"]["available_groups"]

def test_missing_and_duplicate_genes(temp_dir):
    """Test handling of duplicate gene symbols and empty/missing gene names."""
    expr_data = {
        "Gene": ["TP53", "TP53", "EGFR", "", "NA", "PTEN"],
        "S1": [8.5, 9.0, 7.2, 5.0, 6.1, 4.3],
        "S2": [8.1, 8.9, 7.5, 5.2, 6.0, 4.5],
        "S3": [4.2, 4.1, 8.0, 5.1, 6.3, 7.8],
        "S4": [4.0, 4.3, 8.2, 5.0, 6.2, 7.9]
    }
    expr_df = pd.DataFrame(expr_data)
    meta_df = pd.DataFrame({
        "Sample": ["S1", "S2", "S3", "S4"],
        "Group": ["Tumor", "Tumor", "Normal", "Normal"]
    })

    expr_path = os.path.join(temp_dir, "dup_expr.csv")
    meta_path = os.path.join(temp_dir, "dup_meta.csv")
    expr_df.to_csv(expr_path, index=False)
    meta_df.to_csv(meta_path, index=False)

    res = validate_and_preprocess(expr_path, meta_path)
    assert res["success"] is True
    # Duplicate TP53 collapsed into 1 row, invalid gene names dropped
    assert res["info"]["num_genes"] == 3  # TP53, EGFR, PTEN
    
    # Check processed expression file has mean TP53: (8.5+9.0)/2 = 8.75
    proc_df = pd.read_csv(res["processed_expr_path"], index_col=0)
    assert "TP53" in proc_df.index
    assert abs(proc_df.loc["TP53", "S1"] - 8.75) < 1e-4

def test_invalid_numerical_values_and_imputation(temp_dir):
    """Test coercion of non-numeric values and median imputation."""
    expr_data = {
        "Gene": ["TP53", "EGFR", "PTEN", "MYC"],
        "S1": [8.5, "INVALID_STR", 7.2, 5.0],
        "S2": [8.1, 6.5, np.nan, 5.2],
        "S3": [4.2, 7.0, 8.0, 4.8],
        "S4": [4.0, 7.2, 7.9, 5.1]
    }
    expr_df = pd.DataFrame(expr_data)
    meta_df = pd.DataFrame({
        "Sample": ["S1", "S2", "S3", "S4"],
        "Group": ["Tumor", "Tumor", "Normal", "Normal"]
    })

    expr_path = os.path.join(temp_dir, "nan_expr.csv")
    meta_path = os.path.join(temp_dir, "nan_meta.csv")
    expr_df.to_csv(expr_path, index=False)
    meta_df.to_csv(meta_path, index=False)

    res = validate_and_preprocess(expr_path, meta_path, handle_missing="impute_median")
    assert res["success"] is True
    # Check warnings exist
    assert any("non-numeric" in w.lower() for w in res["warnings"])
    
    # Check processed expression file has no NaNs
    proc_df = pd.read_csv(res["processed_expr_path"], index_col=0)
    assert not proc_df.isna().any().any()

def test_mismatched_samples(temp_dir):
    """Test detection of mismatched metadata and expression sample identifiers."""
    expr_df = pd.DataFrame({
        "Gene": ["TP53", "EGFR"],
        "Sample_A": [8.0, 6.0],
        "Sample_B": [7.5, 6.2]
    })
    meta_df = pd.DataFrame({
        "Sample": ["Sample_X", "Sample_Y"],
        "Group": ["Tumor", "Normal"]
    })

    expr_path = os.path.join(temp_dir, "mismatch_expr.csv")
    meta_path = os.path.join(temp_dir, "mismatch_meta.csv")
    expr_df.to_csv(expr_path, index=False)
    meta_df.to_csv(meta_path, index=False)

    res = validate_and_preprocess(expr_path, meta_path)
    # Should report error due to insufficient overlapping samples
    assert res["success"] is False
    assert any("insufficient overlap" in e.lower() for e in res["errors"])

def test_rnaseq_count_auto_inference(temp_dir):
    """Test auto-inference of integer raw RNA-seq count data."""
    expr_df = pd.DataFrame({
        "Gene": ["TP53", "EGFR", "PTEN", "MYC"],
        "S1": [1200, 450, 890, 3200],
        "S2": [1350, 420, 840, 3100],
        "S3": [310, 850, 1600, 1200],
        "S4": [290, 890, 1550, 1150]
    })
    meta_df = pd.DataFrame({
        "Sample": ["S1", "S2", "S3", "S4"],
        "Group": ["Tumor", "Tumor", "Normal", "Normal"]
    })

    expr_path = os.path.join(temp_dir, "counts_expr.csv")
    meta_path = os.path.join(temp_dir, "counts_meta.csv")
    expr_df.to_csv(expr_path, index=False)
    meta_df.to_csv(meta_path, index=False)

    res = validate_and_preprocess(expr_path, meta_path, input_type="auto")
    assert res["success"] is True
    assert res["info"]["inferred_type"] == "rnaseq_counts"
