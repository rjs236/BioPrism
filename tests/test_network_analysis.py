import os
import pandas as pd
from unittest.mock import patch, MagicMock
from pipeline.network_analysis import (
    fetch_string_network,
    build_and_analyze_network
)

def test_string_network_mock_success():
    """Test STRING API data parsing when API returns valid protein interactions."""
    mock_response = [
        {
            "preferredName_A": "TP53",
            "preferredName_B": "MDM2",
            "score": 0.99
        },
        {
            "preferredName_A": "TP53",
            "preferredName_B": "CDKN1A",
            "score": 0.98
        }
    ]

    with patch("requests.post") as mock_post:
        mock_obj = MagicMock()
        mock_obj.status_code = 200
        mock_obj.json.return_value = mock_response
        mock_post.return_value = mock_obj

        edges, msg = fetch_string_network(["TP53", "MDM2", "CDKN1A"], species=9606, confidence_score=400)
        assert len(edges) == 2
        assert edges[0]["source"] == "TP53"
        assert edges[0]["target"] == "MDM2"
        assert edges[0]["score"] == 990.0

def test_string_network_api_failure_and_fallback():
    """Test STRING API timeout / connection failure handling."""
    with patch("requests.post", side_effect=Exception("Connection timed out")):
        # Request with genes present in local curated fallback
        edges, msg = fetch_string_network(["TP53", "MDM2", "CDKN1A"], species=9606, confidence_score=400)
        # Should gracefully return edges from curated fallback or informative message
        assert isinstance(edges, list)
        assert isinstance(msg, str)

def test_network_centrality_and_hub_ranking(temp_dir):
    """Test NetworkX graph building, centrality calculations, and hub gene ranking."""
    sig_genes = ["TP53", "MDM2", "EGFR", "PTEN", "MYC", "CDKN1A", "AKT1"]
    de_df = pd.DataFrame({
        "Gene": sig_genes,
        "log2FoldChange": [3.2, 2.5, -2.8, -2.1, 2.9, 1.8, 2.0],
        "adj_p_val": [1e-6, 1e-5, 1e-5, 1e-4, 1e-5, 1e-3, 1e-4],
        "Regulation": ["Upregulated", "Upregulated", "Downregulated", "Downregulated", "Upregulated", "Upregulated", "Upregulated"]
    })

    res = build_and_analyze_network(
        sig_genes=sig_genes,
        de_df=de_df,
        species=9606,
        confidence_score=400,
        output_dir=temp_dir
    )

    assert res["success"] is True
    assert "summary" in res
    assert "hub_genes" in res
    assert "network_json" in res

    # Check network statistics
    summary = res["summary"]
    assert summary["num_nodes"] >= 2
    assert summary["num_edges"] >= 1
    assert "density" in summary

    # Check hub gene candidates
    hub_genes = res["hub_genes"]
    assert len(hub_genes) > 0
    top_hub = hub_genes[0]
    assert "Gene" in top_hub
    assert "Composite_Hub_Score" in top_hub
    assert "Degree" in top_hub
    assert "Betweenness_Centrality" in top_hub
    assert top_hub["Hub_Rank"] == 1

    # Check exported CSV files
    edges_csv = os.path.join(temp_dir, "network_edges.csv")
    hubs_csv = os.path.join(temp_dir, "ranked_hub_genes.csv")
    assert os.path.exists(edges_csv)
    assert os.path.exists(hubs_csv)
