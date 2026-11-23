import os
import io
from pipeline.report_generator import generate_html_report

def test_html_report_generation(temp_dir):
    """Test standalone HTML report generation with provenance and tables."""
    report_path = os.path.join(temp_dir, "test_report.html")
    
    run_meta = {"run_id": "TEST-1234", "dataset": "Synthetic Benchmark"}
    preproc = {"inferred_type": "Microarray Log2", "preprocessing_applied": "Log2 Normalized", "num_genes": 25, "num_samples": 10}
    de_summary = {
        "engine": "Welch t-test + BH-FDR",
        "case_group": "Tumor",
        "control_group": "Normal",
        "fdr_threshold": 0.05,
        "log2fc_threshold": 1.0,
        "significant_total": 4,
        "upregulated": 2,
        "downregulated": 2,
        "top_de_genes": [{"Gene": "TP53", "log2FoldChange": 3.5, "p_value": 0.0001, "adj_p_val": 0.001, "Regulation": "Upregulated"}]
    }
    pathway_res = [{"term": "Pathways in Cancer", "adj_p_val": 0.002, "odds_ratio": 4.5, "genes": "TP53, MDM2"}]
    net_sum = {"num_nodes": 6, "num_edges": 8, "density": 0.25}
    hubs = [{"Hub_Rank": 1, "Gene": "TP53", "Composite_Hub_Score": 0.95, "Degree": 4, "Betweenness_Centrality": 0.35, "log2FoldChange": 3.5, "adj_p_val": 0.001}]
    ml_metrics = {
        "success": True,
        "metrics": {
            "roc_auc": 0.95,
            "accuracy": 0.90,
            "precision": 0.85,
            "recall": 0.90,
            "f1_score": 0.875,
            "cv_method": "5-Fold Stratified Cross-Validation",
            "classifier": "Logistic Regression",
            "confusion_matrix": {"tp": 5, "tn": 4, "fp": 1, "fn": 0}
        }
    }

    out = generate_html_report(
        run_metadata=run_meta,
        preproc_info=preproc,
        de_summary=de_summary,
        pathway_results=pathway_res,
        network_summary=net_sum,
        hub_genes=hubs,
        ml_metrics=ml_metrics,
        output_filepath=report_path
    )

    assert os.path.exists(out)
    with open(out, "r", encoding="utf-8") as f:
        content = f.read()
    
    assert "BioPrism Analysis Report" in content
    assert "TEST-1234" in content
    assert "TP53" in content
    assert "Pathways in Cancer" in content
    assert "Print / Export PDF" in content
    assert "Scientific & Methodological Disclaimers" in content

def test_flask_web_routes(flask_client):
    """Test core Flask HTML pages load properly with 200 OK."""
    # Homepage
    resp = flask_client.get("/")
    assert resp.status_code == 200
    assert b"BioPrism" in resp.data

    # Pipeline Workspace
    resp = flask_client.get("/pipeline")
    assert resp.status_code == 200
    assert b"Analysis Workspace" in resp.data

    # Methodology / Help page
    resp = flask_client.get("/help")
    assert resp.status_code == 200
    assert b"Methodology" in resp.data

    # All dedicated module pages
    for route in [
        "/upload", "/preprocessing", "/differential-expression",
        "/pathway-enrichment", "/network-analysis", "/biomarker-validation",
        "/ai-assistant", "/reports", "/history", "/settings"
    ]:
        r = flask_client.get(route)
        assert r.status_code == 200, f"Route {route} failed with status {r.status_code}"

def test_flask_ai_chat_api(flask_client):
    """Test AI assistant chat endpoint."""
    resp = flask_client.post("/api/ai/chat", json={
        "session_id": "test_ai_sess",
        "message": "Explain why TP53 was prioritized"
    })
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["success"] is True
    assert "response" in data
    assert len(data["response"]) > 0

def test_flask_demo_and_download_routes(flask_client):
    """Test loading demo data and downloading templates."""
    # Load demo dataset
    resp = flask_client.post("/api/load-demo", json={"type": "microarray"})
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["success"] is True
    assert "session_id" in data

    # Download template expression
    resp_dl = flask_client.get("/api/download/template_expression")
    assert resp_dl.status_code == 200
    assert b"Gene" in resp_dl.data

    # Download template metadata
    resp_meta = flask_client.get("/api/download/template_metadata")
    assert resp_meta.status_code == 200
    assert b"Sample" in resp_meta.data

    # Download non-existent file
    resp_404 = flask_client.get("/api/download/invalid_file_type")
    assert resp_404.status_code == 404
