import os
import json
import urllib.parse
import requests
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Built-in curated high-confidence human interactions fallback (STRING / BioGRID reference)
BUILTIN_PPI_INTERACTIONS = [
    ("TP53", "MDM2", 990, "Experimental Assay"),
    ("TP53", "CDKN1A", 985, "Experimental Assay"),
    ("TP53", "ATM", 980, "Experimental Assay"),
    ("TP53", "CHEK2", 975, "Experimental Assay"),
    ("TP53", "BAX", 970, "Experimental Assay"),
    ("TP53", "PTEN", 950, "Curated Database"),
    ("TP53", "BRCA1", 965, "Experimental Assay"),
    ("TP53", "MYC", 940, "Co-expression / Text mining"),
    ("EGFR", "ERBB2", 995, "Experimental Assay"),
    ("EGFR", "GRB2", 990, "Experimental Assay"),
    ("EGFR", "PIK3CA", 970, "Experimental Assay"),
    ("EGFR", "AKT1", 960, "Curated Database"),
    ("EGFR", "STAT3", 940, "Experimental Assay"),
    ("EGFR", "VEGFA", 920, "Predicted / Co-expression"),
    ("ERBB2", "PIK3CA", 980, "Experimental Assay"),
    ("ERBB2", "AKT1", 975, "Experimental Assay"),
    ("ERBB2", "GRB2", 985, "Experimental Assay"),
    ("MYC", "CCND1", 960, "Curated Database"),
    ("MYC", "CDK4", 950, "Experimental Assay"),
    ("MYC", "RB1", 940, "Experimental Assay"),
    ("MYC", "TP53", 940, "Co-expression / Text mining"),
    ("CCND1", "CDK4", 995, "Experimental Assay"),
    ("CCND1", "RB1", 990, "Experimental Assay"),
    ("CDK4", "RB1", 985, "Experimental Assay"),
    ("CDK4", "CDKN1A", 970, "Experimental Assay"),
    ("PTEN", "AKT1", 990, "Experimental Assay"),
    ("PTEN", "PIK3CA", 985, "Experimental Assay"),
    ("AKT1", "MTOR", 990, "Experimental Assay"),
    ("AKT1", "BAD", 970, "Experimental Assay"),
    ("AKT1", "MDM2", 950, "Experimental Assay"),
    ("VEGFA", "KDR", 995, "Experimental Assay"),
    ("VEGFA", "FLT1", 990, "Experimental Assay"),
    ("VEGFA", "HIF1A", 960, "Experimental Assay"),
    ("BRCA1", "BRCA2", 995, "Experimental Assay"),
    ("BRCA1", "ATM", 985, "Experimental Assay"),
    ("BRCA1", "RAD51", 980, "Experimental Assay"),
    ("BRCA1", "CHEK1", 970, "Experimental Assay"),
    ("CASP3", "CASP8", 990, "Experimental Assay"),
    ("CASP3", "CASP9", 985, "Experimental Assay"),
    ("CASP3", "BAX", 960, "Experimental Assay"),
    ("CASP3", "BCL2", 950, "Experimental Assay"),
    ("BCL2", "BAX", 995, "Experimental Assay"),
    ("CD274", "PDCD1", 995, "Experimental Assay"),
    ("CD80", "CTLA4", 995, "Experimental Assay"),
    ("CD86", "CTLA4", 990, "Experimental Assay"),
    ("TNF", "IL6", 960, "Curated Database"),
    ("IL6", "STAT3", 990, "Experimental Assay"),
    ("IFNG", "STAT1", 995, "Experimental Assay")
]

def fetch_string_network(genes, species=9606, confidence_score=400):
    """
    Query the STRING Database API for protein-protein functional associations.
    Distinguishes experimentally established interactions from predicted associations.
    """
    cleaned_genes = [str(g).strip().upper() for g in genes if str(g).strip()]
    if not cleaned_genes:
        return [], "No valid genes provided."

    query_genes = cleaned_genes[:150]
    string_api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(query_genes),
        "species": int(species),
        "required_score": int(confidence_score),
        "caller_identity": "BioPrism_Pipeline"
    }

    try:
        response = requests.post(string_api_url, data=params, timeout=12)
        if response.status_code == 200:
            data = response.json()
            if isinstance(data, list) and len(data) > 0:
                edges = []
                for item in data:
                    pA = item.get('preferredName_A', item.get('stringId_A', ''))
                    pB = item.get('preferredName_B', item.get('stringId_B', ''))
                    score = float(item.get('score', 0.4)) * 1000.0
                    
                    # Inspect sub-scores if present to distinguish experimental vs predicted
                    escore = float(item.get('escore', 0.0)) # experimental score
                    dscore = float(item.get('dscore', 0.0)) # database score
                    if escore > 0.4:
                        evidence_type = "Experimentally Validated"
                    elif dscore > 0.4:
                        evidence_type = "Curated Knowledgebase"
                    else:
                        evidence_type = "Predicted Functional Association"

                    if pA and pB and pA != pB:
                        edges.append({
                            'source': pA,
                            'target': pB,
                            'score': round(score, 1),
                            'evidence': evidence_type
                        })
                return edges, f"Retrieved {len(edges)} interactions from STRING DB API (Species: {species}, Min Score: {confidence_score})."
            else:
                return [], "STRING API returned no interactions above the selected confidence threshold."
        else:
            return [], f"STRING API returned HTTP status {response.status_code}."
    except Exception as e:
        # Fallback to local curated interactions
        query_set = set(query_genes)
        fallback_edges = []
        for pA, pB, score, ev in BUILTIN_PPI_INTERACTIONS:
            if pA in query_set and pB in query_set and score >= confidence_score:
                fallback_edges.append({
                    'source': pA,
                    'target': pB,
                    'score': score,
                    'evidence': ev
                })
        if fallback_edges:
            return fallback_edges, f"Network query used BioPrism Curated Reference PPI (API connection unreachable: {str(e)[:50]})."
        return [], f"Unable to reach STRING API and no local PPI matches found: {str(e)[:80]}"


def build_and_analyze_network(sig_genes, de_df, species=9606, confidence_score=400, output_dir=None):
    """
    Module 3: Construct graph with NetworkX, calculate centralities, rank hub genes,
    and generate interactive Plotly network visualization.
    """
    if output_dir is None:
        output_dir = "."
    os.makedirs(output_dir, exist_ok=True)

    edges, source_msg = fetch_string_network(sig_genes, species=species, confidence_score=confidence_score)

    de_dict = {}
    if de_df is not None and not de_df.empty:
        for _, row in de_df.iterrows():
            g = str(row['Gene']).upper()
            de_dict[g] = {
                'log2FoldChange': float(row.get('log2FoldChange', 0.0)),
                'adj_p_val': float(row.get('adj_p_val', 1.0)),
                'p_value': float(row.get('p_value', 1.0)),
                'Regulation': str(row.get('Regulation', 'Not Significant'))
            }

    G = nx.Graph()

    for g in sig_genes:
        g_clean = str(g).strip().upper()
        info = de_dict.get(g_clean, {
            'log2FoldChange': 0.0,
            'adj_p_val': 1.0,
            'p_value': 1.0,
            'Regulation': 'Significant'
        })
        G.add_node(g_clean, **info)

    for edge in edges:
        s, t = edge['source'].upper(), edge['target'].upper()
        if s not in G:
            G.add_node(s, log2FoldChange=0.0, adj_p_val=1.0, p_value=1.0, Regulation='Interacting Partner')
        if t not in G:
            G.add_node(t, log2FoldChange=0.0, adj_p_val=1.0, p_value=1.0, Regulation='Interacting Partner')
        G.add_edge(s, t, weight=edge['score'] / 1000.0, score=edge['score'], evidence=edge['evidence'])

    if G.number_of_nodes() == 0:
        return {
            'success': False,
            'message': 'No nodes available to build network.',
            'hub_genes': [],
            'network_json': None,
            'edges': []
        }

    # Compute Centralities
    degree_dict = dict(G.degree())
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G, normalized=True)
    closeness_centrality = nx.closeness_centrality(G)

    # Composite Hub Ranking Formula
    max_deg = max(degree_dict.values()) if degree_dict.values() else 1
    hub_candidates = []

    for node in G.nodes():
        deg = degree_dict.get(node, 0)
        deg_cent = degree_centrality.get(node, 0.0)
        bet_cent = betweenness_centrality.get(node, 0.0)
        close_cent = closeness_centrality.get(node, 0.0)
        
        node_data = G.nodes[node]
        fc = abs(node_data.get('log2FoldChange', 0.0))
        fdr = node_data.get('adj_p_val', 1.0)
        stat_significance = min(1.0, (fc * (-np.log10(max(fdr, 1e-15)))) / 10.0)

        # Composite Hub Score: Degree (40%) + Betweenness (35%) + Statistical Significance (25%)
        hub_score = (0.40 * (deg / max_deg)) + (0.35 * bet_cent) + (0.25 * stat_significance)

        hub_candidates.append({
            'Gene': node,
            'Composite_Hub_Score': round(hub_score, 4),
            'Degree': deg,
            'Degree_Centrality': round(deg_cent, 4),
            'Betweenness_Centrality': round(bet_cent, 4),
            'Closeness_Centrality': round(close_cent, 4),
            'log2FoldChange': round(node_data.get('log2FoldChange', 0.0), 3),
            'adj_p_val': node_data.get('adj_p_val', 1.0),
            'Regulation': node_data.get('Regulation', 'Not Significant')
        })

    hub_candidates.sort(key=lambda x: x['Composite_Hub_Score'], reverse=True)
    for idx, h in enumerate(hub_candidates):
        h['Hub_Rank'] = idx + 1

    # Save outputs to CSV
    hubs_df = pd.DataFrame(hub_candidates)
    hubs_csv_path = os.path.join(output_dir, "ranked_hub_genes.csv")
    hubs_df.to_csv(hubs_csv_path, index=False)

    edges_df = pd.DataFrame(edges)
    edges_csv_path = os.path.join(output_dir, "network_edges.csv")
    if not edges_df.empty:
        edges_df.to_csv(edges_csv_path, index=False)

    # Generate 2D Force-Directed Layout via NetworkX spring_layout
    pos = nx.spring_layout(G, k=0.45, iterations=60, seed=42)

    # Edge traces
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.0, color='#cbd5e1'),
        hoverinfo='none',
        mode='lines'
    )

    # Node trace
    node_x = []
    node_y = []
    node_text = []
    node_color = []
    node_size = []

    top_hub_genes_set = set([h['Gene'] for h in hub_candidates[:5]])

    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        
        node_data = G.nodes[node]
        reg = node_data.get('Regulation', 'Not Significant')
        deg = degree_dict.get(node, 0)
        bet = betweenness_centrality.get(node, 0.0)
        fc = node_data.get('log2FoldChange', 0.0)
        fdr = node_data.get('adj_p_val', 1.0)

        # Highlight top hub genes with distinct styling
        is_hub = node in top_hub_genes_set
        
        # Color coding: Teal for Up, Navy/Blue for Down, Slate for other
        if reg == 'Upregulated':
            color = '#0d9488'
        elif reg == 'Downregulated':
            color = '#2563eb'
        else:
            color = '#64748b'
            
        node_color.append(color)
        
        # Size by degree centrality
        size = 14 + (deg / max_deg) * 26
        if is_hub:
            size += 6
        node_size.append(size)

        hub_label = " ⭐ [TOP CANDIDATE HUB]" if is_hub else ""
        node_text.append(
            f"<b>{node}</b>{hub_label}<br>"
            f"Regulation: {reg}<br>"
            f"Degree: {deg}<br>"
            f"Betweenness Centrality: {bet:.4f}<br>"
            f"log2FC: {fc:.2f}<br>"
            f"FDR: {fdr:.2e}"
        )

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=[node if node in top_hub_genes_set or degree_dict.get(node, 0) >= 4 else "" for node in G.nodes()],
        textposition="top center",
        hovertext=node_text,
        marker=dict(
            color=node_color,
            size=node_size,
            line=dict(width=2, color='#ffffff')
        )
    )

    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            title="Protein-Protein Interaction Network (STRING Database & NetworkX Topology)",
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=20, r=20, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            template="plotly_white"
        )
    )
    network_json = json.loads(fig.to_json())

    # Summary dictionary
    summary = {
        'num_nodes': G.number_of_nodes(),
        'total_nodes': G.number_of_nodes(),
        'num_edges': G.number_of_edges(),
        'total_edges': G.number_of_edges(),
        'density': round(nx.density(G), 4),
        'connected_components': nx.number_connected_components(G),
        'source_message': source_msg,
        'confidence_score_threshold': confidence_score,
        'top_hubs_preview': hub_candidates[:8],
        'edges_csv_path': edges_csv_path,
        'hubs_csv_path': hubs_csv_path
    }

    return {
        'success': True,
        'message': f"Network built with {G.number_of_nodes()} proteins and {G.number_of_edges()} interactions.",
        'summary': summary,
        'hub_genes': hub_candidates,
        'network_json': network_json,
        'edges': edges
    }
