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


def test_reports_no_auto_download_and_results_persistence(flask_client):
    """Ensure /reports does not auto-download, inline preview works for PDF print, HTML/CSV exports include results, and history records runs."""
    sid = "reg_test_sess"

    # 1. /reports page must NOT embed attachment download url in iframe src
    r_page = flask_client.get(f"/reports?session_id={sid}")
    assert r_page.status_code == 200
    assert b'iframe id="reportFrame" src="/api/download/html_report' not in r_page.data

    # 2. Run pipeline steps via Flask APIs
    r_demo = flask_client.post("/api/load-demo", json={"type": "microarray", "session_id": sid})
    assert r_demo.status_code == 200
    assert r_demo.get_json()["success"] is True

    r_de = flask_client.post("/api/run-de", json={
        "session_id": sid,
        "case_group": "Tumor",
        "control_group": "Normal",
        "method": "python_ttest",
        "fdr_threshold": 0.05,
        "log2fc_threshold": 1.0
    })
    assert r_de.status_code == 200
    de_json = r_de.get_json()
    assert de_json["success"] is True
    snap = de_json["session_snapshot"]

    r_net = flask_client.post("/api/run-network", json={"session_id": sid, "confidence_score": 400})
    assert r_net.status_code == 200
    net_json = r_net.get_json()
    assert net_json["success"] is True
    snap = net_json["session_snapshot"]

    # 3. Verify Analysis History already contains the run after DE/Network
    r_hist = flask_client.get("/history")
    assert r_hist.status_code == 200
    assert sid.encode() in r_hist.data
    assert b"Oncogenomics Microarray Demo" in r_hist.data

    # 4. Simulate serverless cold container (wipe SESSIONS in-memory dict)
    from app import SESSIONS
    SESSIONS.pop(sid, None)

    # 5. Inline preview must NOT have Content-Disposition: attachment and MUST contain results
    r_prev = flask_client.get(f"/api/preview/html_report?session_id={sid}")
    assert r_prev.status_code == 200
    assert "attachment" not in r_prev.headers.get("Content-Disposition", "").lower()
    assert b"Oncogenomics Microarray Demo" in r_prev.data
    assert b"45" in r_prev.data  # 45 significant DEGs

    # 6. POST /api/generate-report with snapshot works even if both RAM and /tmp are cold
    SESSIONS.pop(sid, None)
    r_gen = flask_client.post("/api/generate-report", json={"session_id": sid, "session_snapshot": snap})
    assert r_gen.status_code == 200
    gen_data = r_gen.get_json()
    assert gen_data["success"] is True
    assert "Oncogenomics Microarray Demo" in gen_data["report_html"]

    # 7. Explicit HTML download & CSV downloads with POST snapshot succeed with real data
    r_dl = flask_client.post(f"/api/download/html_report?session_id={sid}", json={"session_id": sid, "session_snapshot": snap})
    assert r_dl.status_code == 200
    assert "attachment" in r_dl.headers.get("Content-Disposition", "").lower()
    assert b"Oncogenomics Microarray Demo" in r_dl.data

    r_de_csv = flask_client.post(f"/api/download/de_results?session_id={sid}", json={"session_id": sid, "session_snapshot": snap})
    assert r_de_csv.status_code == 200
    assert b"log2FoldChange" in r_de_csv.data

    # 8. Cross-container merge check: ensure running validation preserves network_summary & hub_genes
    SESSIONS.pop(sid, None)
    r_val = flask_client.post("/api/run-validation", json={
        "session_id": sid,
        "candidate_genes": ["TP53", "AKT1"],
        "classifier": "logistic_regression",
        "cv_splits": 3,
        "session_snapshot": snap
    })
    assert r_val.status_code == 200
    val_json = r_val.get_json()
    assert val_json["success"] is True
    merged_snap = val_json["session_snapshot"]
    assert merged_snap.get("network_summary", {}).get("total_nodes", 0) > 0
    assert len(merged_snap.get("hub_genes", [])) > 0

    r_gen2 = flask_client.post("/api/generate-report", json={"session_id": sid, "session_snapshot": merged_snap})
    assert r_gen2.status_code == 200
    html2 = r_gen2.get_json()["report_html"]
    assert "0 nodes, 0 edges" not in html2
    assert "No hub genes identified" not in html2

    # Clean up test run from history
    flask_client.post("/api/history/delete", json={"id": sid})
