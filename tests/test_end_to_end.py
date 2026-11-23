import os
import pandas as pd
from pipeline.preprocessor import validate_and_preprocess
from pipeline.de_analysis import run_differential_expression, run_pathway_enrichment
from pipeline.network_analysis import build_and_analyze_network
from pipeline.validation_ml import run_biomarker_classification
from pipeline.report_generator import generate_html_report

def test_full_pipeline_end_to_end(synthetic_data, temp_dir):
    """
    End-to-end integration test for the entire BioPrism pipeline.
    NOTE: The data utilized in this test is EXPLICITLY SYNTHETIC and is solely
    designed to verify software integration and pipeline flow without claiming
    biological or clinical discoveries.
    """
    expr_file = synthetic_data["expr_path"]
    meta_file = synthetic_data["meta_path"]
    case_grp = synthetic_data["case_group"]
    ctrl_grp = synthetic_data["control_group"]

    # 1. Module 1: Preprocessing
    step1_res = validate_and_preprocess(expr_file, meta_file)
    assert step1_res["success"] is True
    proc_expr = step1_res["processed_expr_path"]
    proc_meta = step1_res["processed_meta_path"]
    assert os.path.exists(proc_expr)
    assert os.path.exists(proc_meta)

    # 2. Module 2: Differential Expression & Pathway Enrichment
    de_df, volcano_json, de_summary = run_differential_expression(
        expr_csv=proc_expr,
        meta_csv=proc_meta,
        case_group=case_grp,
        control_group=ctrl_grp,
        fdr_threshold=0.05,
        log2fc_threshold=1.0,
        output_dir=temp_dir
    )
    assert not de_df.empty
    assert "log2FoldChange" in de_df.columns
    assert "adj_p_val" in de_df.columns

    sig_genes = de_summary.get("sig_genes", [])
    if not sig_genes:
        sig_genes = de_df.head(10)["Gene"].tolist()

    enrich_res = run_pathway_enrichment(sig_genes, gene_set="KEGG_2021_Human")
    assert enrich_res["success"] is True

    # 3. Module 3: Network Analysis & Hub Gene Identification
    net_res = build_and_analyze_network(
        sig_genes=sig_genes,
        de_df=de_df,
        species=9606,
        confidence_score=400,
        output_dir=temp_dir
    )
    assert net_res["success"] is True
    hub_genes = net_res["hub_genes"]
    assert len(hub_genes) > 0

    # 4. Module 4: Machine Learning Validation (Model Comparison)
    candidate_genes = [h["Gene"] for h in hub_genes[:6]]
    ml_res = run_biomarker_classification(
        expr_csv=proc_expr,
        meta_csv=proc_meta,
        candidate_genes=candidate_genes,
        classifier_type="compare_both",
        n_splits=3,
        output_dir=temp_dir
    )
    assert ml_res["success"] is True
    assert ml_res["metrics"]["roc_auc"] >= 0.5
    assert ml_res["metrics"]["comparison"] is not None

    # 5. Final Report Generation
    report_file = os.path.join(temp_dir, "BioPrism_EndToEnd_Report.html")
    report_path = generate_html_report(
        run_metadata={"run_id": "E2E-TEST", "dataset": "Synthetic Verification Dataset"},
        preproc_info=step1_res["info"],
        de_summary=de_summary,
        pathway_results=enrich_res.get("results", []),
        network_summary=net_res["summary"],
        hub_genes=hub_genes,
        ml_metrics=ml_res,
        output_filepath=report_file
    )
    assert os.path.exists(report_path)
    assert os.path.getsize(report_path) > 1000
