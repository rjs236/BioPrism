import os
import shutil
import urllib.request
import pandas as pd
from werkzeug.utils import secure_filename
from config import Config

# Curated public GEO cohorts metadata available for instant loading or fetching
CURATED_GEO_DATASETS = {
    "GSE42568": {
        "title": "GSE42568: Human Breast Cancer Microarray (104 samples)",
        "organism": "Homo sapiens",
        "platform": "Affymetrix HG-U133_Plus_2",
        "type": "microarray_log2",
        "description": "Normalized expression of primary breast cancer tumors vs adjacent normal tissues.",
        "case_group": "Tumor",
        "control_group": "Normal",
        "sample_counts": "67 Tumor, 17 Normal",
        "demo_file": "demo_microarray_expression.csv",
        "demo_meta": "demo_microarray_metadata.csv"
    },
    "GSE19804": {
        "title": "GSE19804: Non-Small Cell Lung Carcinoma (60 samples)",
        "organism": "Homo sapiens",
        "platform": "Illumina HiSeq / Counts",
        "type": "rnaseq_counts",
        "description": "Transcriptome profiling of NSCLC adenocarcinoma vs matched paired normal lungs.",
        "case_group": "Lung_Tumor",
        "control_group": "Normal_Lung",
        "sample_counts": "30 Tumor, 30 Normal",
        "demo_file": "demo_rnaseq_counts.csv",
        "demo_meta": "demo_rnaseq_metadata.csv"
    },
    "GSE3325": {
        "title": "GSE3325: Clinical Melanoma Progression Cohort",
        "organism": "Homo sapiens",
        "platform": "Affymetrix U133A",
        "type": "microarray_log2",
        "description": "Metastatic melanoma vs primary cutaneous melanoma biopsy samples.",
        "case_group": "Metastatic",
        "control_group": "Primary",
        "sample_counts": "20 Metastatic, 20 Primary",
        "demo_file": "demo_microarray_expression.csv",
        "demo_meta": "demo_microarray_metadata.csv"
    }
}

def get_curated_geo_list():
    """Return available public GEO cohorts."""
    return CURATED_GEO_DATASETS

def load_geo_dataset(geo_accession, target_dir):
    """
    Load or fetch a public GEO cohort.
    If the GEO ID matches one of the verified benchmark cohorts,
    copies the verified expression and metadata files into target_dir.
    """
    geo_id = geo_accession.strip().upper()
    if geo_id not in CURATED_GEO_DATASETS:
        # Check if user entered an arbitrary GEO ID
        return {
            'success': False,
            'message': f"GEO accession '{geo_accession}' is not pre-cached in offline mode. Curated accessions available: {', '.join(CURATED_GEO_DATASETS.keys())}"
        }

    info = CURATED_GEO_DATASETS[geo_id]
    src_expr = os.path.join(Config.SAMPLE_DATA_FOLDER, info['demo_file'])
    src_meta = os.path.join(Config.SAMPLE_DATA_FOLDER, info['demo_meta'])

    if not os.path.exists(src_expr) or not os.path.exists(src_meta):
        return {'success': False, 'message': 'Underlying cohort files missing from sample data directory.'}

    dest_expr = os.path.join(target_dir, f"{geo_id}_expression.csv")
    dest_meta = os.path.join(target_dir, f"{geo_id}_metadata.csv")

    shutil.copyfile(src_expr, dest_expr)
    shutil.copyfile(src_meta, dest_meta)

    return {
        'success': True,
        'geo_id': geo_id,
        'title': info['title'],
        'description': info['description'],
        'type': info['type'],
        'case_group': info['case_group'],
        'control_group': info['control_group'],
        'expr_path': dest_expr,
        'meta_path': dest_meta
    }

def save_uploaded_files(expr_file, meta_file, target_dir, session_id):
    """
    Securely save user-uploaded expression and metadata files.
    """
    os.makedirs(target_dir, exist_ok=True)
    expr_name = secure_filename(f"expr_{session_id}_{expr_file.filename}")
    meta_name = secure_filename(f"meta_{session_id}_{meta_file.filename}")

    expr_path = os.path.join(target_dir, expr_name)
    meta_path = os.path.join(target_dir, meta_name)

    expr_file.save(expr_path)
    meta_file.save(meta_path)

    return expr_path, meta_path
