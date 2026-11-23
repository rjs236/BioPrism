import os
import sqlite3
import pytest
from database.db import (
    init_db, save_analysis_run, get_dashboard_stats,
    get_all_analyses, get_analysis_by_id, delete_analysis,
    save_chat_message, get_chat_history
)
from modules.data_upload import get_curated_geo_list, load_geo_dataset

def test_database_operations(temp_dir):
    """Test SQLite initialization, record creation, retrieval, and statistics."""
    init_db()

    # Save a run
    run_id = "TEST_RUN_999"
    run_data = {
        'id': run_id,
        'dataset_name': 'Unit Test GSE999',
        'input_type': 'Microarray Log2',
        'total_genes': 120,
        'total_samples': 20,
        'case_group': 'Disease',
        'control_group': 'Healthy',
        'de_method': 'limma',
        'fdr_threshold': 0.05,
        'log2fc_threshold': 1.0,
        'sig_genes_count': 12,
        'n_up': 7,
        'n_down': 5,
        'top_hub_genes': 'TP53, EGFR, MYC',
        'top_pathways': 'Cell Cycle, Pathways in Cancer',
        'classifier_type': 'Logistic Regression',
        'roc_auc': 0.965,
        'accuracy': 0.95,
        'precision_score': 0.92,
        'recall_score': 0.98,
        'f1_score': 0.95
    }

    saved_id = save_analysis_run(run_data)
    assert saved_id == run_id

    # Retrieve by ID
    record = get_analysis_by_id(run_id)
    assert record is not None
    assert record['dataset_name'] == 'Unit Test GSE999'
    assert record['sig_genes_count'] == 12
    assert abs(record['roc_auc'] - 0.965) < 1e-4

    # Check dashboard stats
    stats = get_dashboard_stats()
    assert stats['total_runs'] >= 1
    assert stats['total_datasets'] >= 1
    assert stats['total_de_genes'] >= 12
    assert stats['total_biomarkers'] >= 3  # TP53, EGFR, MYC

    # Chat history operations (using unique session ID)
    import uuid
    sess_key = f"CHAT_TEST_{uuid.uuid4().hex}"
    save_chat_message(sess_key, "user", "What are the top hub genes?")
    save_chat_message(sess_key, "ai", "The top hub genes are TP53 and EGFR.")
    chat_history = get_chat_history(sess_key)
    assert len(chat_history) == 2
    assert chat_history[0]['role'] == 'user'
    assert chat_history[1]['role'] == 'ai'

    # Clean up test run
    del_ok = delete_analysis(run_id)
    assert del_ok is True
    assert get_analysis_by_id(run_id) is None

def test_curated_geo_dataset_loader(temp_dir):
    """Test loading curated GEO benchmarks."""
    geo_list = get_curated_geo_list()
    assert "GSE42568" in geo_list
    assert "GSE19804" in geo_list

    res = load_geo_dataset("GSE42568", temp_dir)
    assert res["success"] is True
    assert os.path.exists(res["expr_path"])
    assert os.path.exists(res["meta_path"])
    assert res["geo_id"] == "GSE42568"
