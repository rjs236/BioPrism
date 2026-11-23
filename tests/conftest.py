import os
import shutil
import tempfile
import pytest
import pandas as pd
import numpy as np
from app import app

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test artifacts and clean up afterwards."""
    d = tempfile.mkdtemp(prefix="bioprism_test_")
    yield d
    shutil.rmtree(d, ignore_errors=True)

@pytest.fixture
def synthetic_data(temp_dir):
    """
    Generate synthetic toy microarray expression and metadata files.
    Clearly labeled as ARTIFICIAL SYNTHETIC DATA.
    """
    np.random.seed(42)
    genes = [f"GENE_{i:02d}" for i in range(1, 26)]
    # Include some known cancer-associated genes for network tests
    genes[:6] = ["TP53", "MDM2", "EGFR", "PTEN", "MYC", "CDKN1A"]

    samples_case = [f"Tumor_{i:02d}" for i in range(1, 6)]
    samples_ctrl = [f"Normal_{i:02d}" for i in range(1, 6)]
    samples = samples_case + samples_ctrl

    # Create expression values (log2 normalized range 3 to 14)
    expr_vals = np.random.normal(loc=7.0, scale=1.5, size=(len(genes), len(samples)))
    
    # Induce strong differential expression in first 4 genes for testing
    expr_vals[0, :5] += 4.0  # TP53 upregulated in tumor
    expr_vals[1, :5] += 3.5  # MDM2 upregulated in tumor
    expr_vals[2, 5:] += 4.0  # EGFR upregulated in normal (down in tumor)
    expr_vals[3, 5:] += 3.5  # PTEN upregulated in normal

    expr_df = pd.DataFrame(expr_vals, index=genes, columns=samples)
    expr_df.index.name = "Gene"

    meta_df = pd.DataFrame({
        "Sample": samples,
        "Group": ["Tumor"] * 5 + ["Normal"] * 5
    })

    expr_path = os.path.join(temp_dir, "synthetic_expression.csv")
    meta_path = os.path.join(temp_dir, "synthetic_metadata.csv")

    expr_df.to_csv(expr_path)
    meta_df.to_csv(meta_path, index=False)

    return {
        "expr_path": expr_path,
        "meta_path": meta_path,
        "genes": genes,
        "samples": samples,
        "case_group": "Tumor",
        "control_group": "Normal"
    }

@pytest.fixture
def flask_client():
    """Configured Flask test client."""
    app.config["TESTING"] = True
    with app.test_client() as client:
        yield client
