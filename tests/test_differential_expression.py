import os
import numpy as np
import pandas as pd
from pipeline.de_analysis import (
    benjamini_hochberg,
    run_differential_expression,
    run_pathway_enrichment
)

def test_benjamini_hochberg_properties():
    """Verify mathematical properties of the Benjamini-Hochberg FDR procedure."""
    # Monotonicity test with known ordered p-values
    p_vals = [0.001, 0.01, 0.03, 0.04, 0.5, 0.8]
    adj = benjamini_hochberg(p_vals)
    
    assert len(adj) == len(p_vals)
    assert np.all(adj >= 0.0)
    assert np.all(adj <= 1.0)
    # Adjusted p-values should be >= raw p-values
    assert np.all(adj >= np.array(p_vals))
    # Order monotonicity: if p_i <= p_j, then adj_i <= adj_j for sorted inputs
    for i in range(len(adj) - 1):
        assert adj[i] <= adj[i + 1]

def test_differential_expression_execution(synthetic_data, temp_dir):
    """Test differential expression execution, fold change and gene categorization."""
    df_results, volcano_json, summary = run_differential_expression(
        expr_csv=synthetic_data["expr_path"],
        meta_csv=synthetic_data["meta_path"],
        case_group=synthetic_data["case_group"],
        control_group=synthetic_data["control_group"],
        fdr_threshold=0.05,
        log2fc_threshold=1.0,
        output_dir=temp_dir
    )

    assert not df_results.empty
    assert "log2FoldChange" in df_results.columns
    assert "p_value" in df_results.columns
    assert "adj_p_val" in df_results.columns

    # Verify induced differentially expressed genes were detected
    # TP53 was induced with +3.5 in tumor
    tp53_row = df_results[df_results["Gene"] == "TP53"]
    assert not tp53_row.empty
    assert tp53_row.iloc[0]["log2FoldChange"] > 1.5

    # Check upregulated vs downregulated separation
    assert summary["upregulated"] >= 1
    assert summary["downregulated"] >= 1
    assert summary["significant_total"] == summary["upregulated"] + summary["downregulated"]

    # Verify interactive volcano plot payload is valid
    assert "data" in volcano_json
    assert "layout" in volcano_json

def test_pathway_enrichment_handling():
    """Test GSEApy pathway enrichment and offline fallback handling."""
    # Test with known set of cancer-related genes
    genes = ["TP53", "EGFR", "PTEN", "MYC", "MDM2", "CDKN1A", "AKT1", "BRCA1"]
    res = run_pathway_enrichment(genes, gene_set="KEGG_2021_Human")

    assert res["success"] is True
    assert "results" in res
    assert isinstance(res["results"], list)

    # Test with empty gene list
    res_empty = run_pathway_enrichment([], gene_set="KEGG_2021_Human")
    assert res_empty["success"] is False
