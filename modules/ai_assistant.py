import os
import json
import urllib.request
import urllib.error

def get_gemini_api_key():
    """Retrieve Gemini API key from environment variables."""
    return os.environ.get('GEMINI_API_KEY') or os.environ.get('GOOGLE_API_KEY') or ''

def build_session_context_summary(session_data):
    """
    Construct a dense, factual scientific context string from the active session data.
    Ensures the AI assistant is strictly grounded in actual computed values.
    """
    if not session_data:
        return "No active analysis session found. User has not uploaded or run any dataset yet."

    dataset_name = session_data.get('dataset_name', 'Unknown Dataset')
    preproc = session_data.get('preproc_info') or {}
    de = session_data.get('de_summary') or {}
    pathways = session_data.get('pathway_results') or []
    network = session_data.get('network_summary') or {}
    hubs = session_data.get('hub_genes') or []
    ml = session_data.get('ml_metrics') or {}

    lines = []
    lines.append(f"DATASET METADATA:")
    lines.append(f"- Dataset Name: {dataset_name}")
    lines.append(f"- Data Type: {preproc.get('inferred_type', 'N/A')}")
    lines.append(f"- Total Filtered Genes: {preproc.get('num_genes', 0)}")
    lines.append(f"- Total Analyzed Samples: {preproc.get('num_samples', 0)}")
    lines.append(f"- Cohort Groups: {preproc.get('group_counts', {})}")

    lines.append(f"\nDIFFERENTIAL EXPRESSION (Module 2):")
    lines.append(f"- Comparison: {de.get('case_group', 'Case')} vs {de.get('control_group', 'Control')}")
    lines.append(f"- Engine Used: {de.get('engine', 'Statistical Engine')}")
    lines.append(f"- Thresholds: FDR < {de.get('fdr_threshold', 0.05)}, |log2FC| >= {de.get('log2fc_threshold', 1.0)}")
    lines.append(f"- Total Significant Genes: {de.get('significant_total', 0)} (Upregulated: {de.get('upregulated', 0)}, Downregulated: {de.get('downregulated', 0)})")
    
    top_de = de.get('top_de_genes', [])[:8]
    if top_de:
        lines.append("- Top DE Genes:")
        for g in top_de:
            lines.append(f"  * {g.get('Gene')}: log2FC = {g.get('log2FoldChange', 0.0):.2f}, FDR = {g.get('adj_p_val', 1.0):.2e}, Status = {g.get('Regulation')}")

    lines.append(f"\nPATHWAY ENRICHMENT (Module 2):")
    if pathways:
        lines.append(f"- Top Enriched Biological Pathways (Total: {len(pathways)}):")
        for p in pathways[:6]:
            lines.append(f"  * {p.get('term')}: FDR = {p.get('adj_p_val', 1.0):.2e}, Odds Ratio = {p.get('odds_ratio', 0.0)}, Genes = {p.get('genes', '')[:60]}")
    else:
        lines.append("- Pathway enrichment has not been executed or yielded no significant terms.")

    lines.append(f"\nPROTEIN INTERACTION NETWORK (Module 3):")
    lines.append(f"- Total Nodes: {network.get('total_nodes', 0)}, Total Edges: {network.get('total_edges', 0)}")
    lines.append(f"- Network Density: {network.get('density', 0.0)}, Confidence Score Threshold: {network.get('confidence_score_threshold', 400)}")
    if hubs:
        lines.append("- Top Ranked Hub Genes (Composite Hub Score):")
        for h in hubs[:6]:
            lines.append(f"  * Rank #{h.get('Hub_Rank')}: {h.get('Gene')} (Hub Score: {h.get('Composite_Hub_Score'):.3f}, Degree: {h.get('Degree')}, Betweenness: {h.get('Betweenness_Centrality'):.4f}, log2FC: {h.get('log2FoldChange'):.2f})")
    else:
        lines.append("- PPI network analysis has not been executed yet.")

    lines.append(f"\nBIOMARKER VALIDATION & ML CLASSIFICATION (Module 4):")
    if ml and ml.get('success'):
        m = ml.get('metrics', {})
        lines.append(f"- Evaluated Model: {m.get('model_name')}")
        lines.append(f"- Validation Strategy: {ml.get('cv_method')}")
        lines.append(f"- Metrics: ROC-AUC = {m.get('roc_auc', 0.0):.3f}, Accuracy = {m.get('accuracy', 0.0):.1%}, Precision = {m.get('precision', 0.0):.1%}, Recall = {m.get('recall', 0.0):.1%}, F1 = {m.get('f1_score', 0.0):.3f}")
        lines.append(f"- Candidate Genes Evaluated: {', '.join(m.get('candidate_genes', []))}")
        cm = m.get('confusion_matrix', {})
        lines.append(f"- Confusion Matrix: TP={cm.get('tp',0)}, TN={cm.get('tn',0)}, FP={cm.get('fp',0)}, FN={cm.get('fn',0)}")
    else:
        lines.append("- Machine learning classification has not been completed yet.")

    return "\n".join(lines)

def generate_deterministic_scientific_response(user_query, session_data):
    """
    Generate rigorous, domain-specific scientific answers directly from the
    active session's computed data when no external LLM API key is present.
    """
    query_lower = user_query.lower()
    de = session_data.get('de_summary') or {}
    hubs = session_data.get('hub_genes') or []
    pathways = session_data.get('pathway_results') or []
    ml = session_data.get('ml_metrics') or {}
    preproc = session_data.get('preproc_info') or {}
    dataset_name = session_data.get('dataset_name', 'your active dataset')

    # Check for specific gene inquiry
    top_gene_names = [h.get('Gene', '').upper() for h in hubs]
    for g in top_gene_names:
        if g.lower() in query_lower:
            hub_info = next((h for h in hubs if h.get('Gene', '').upper() == g), None)
            if hub_info:
                return (
                    f"### 🧬 Biological & Topological Assessment of **{g}**\n\n"
                    f"In the analysis of **{dataset_name}**, gene **{g}** was prioritized as a **Rank #{hub_info.get('Hub_Rank')} candidate hub biomarker** based on the following converging evidence:\n\n"
                    f"1. **Differential Expression Significance:**\n"
                    f"   - **log₂ Fold Change:** `{hub_info.get('log2FoldChange'):.3f}` ({hub_info.get('Regulation')})\n"
                    f"   - **Adjusted p-value (FDR):** `{hub_info.get('adj_p_val'):.2e}`\n"
                    f"   - This indicates a robust statistical change in expression between `{de.get('case_group', 'Case')}` and `{de.get('control_group', 'Control')}`.\n\n"
                    f"2. **Network Centrality & Topological Indispensability:**\n"
                    f"   - **Network Degree:** `{hub_info.get('Degree')}` functional protein interactions\n"
                    f"   - **Betweenness Centrality:** `{hub_info.get('Betweenness_Centrality'):.4f}` (acting as a critical bottleneck/bridge in signal flow)\n"
                    f"   - **Composite Hub Score:** `{hub_info.get('Composite_Hub_Score'):.4f}`\n\n"
                    f"3. **Clinical / Experimental Context:**\n"
                    f"   - While `{g}` demonstrates strong retrospective discriminative potential, computational prioritization constitutes an initial exploratory hypothesis that requires orthogonal RT-qPCR and targeted proteomic validation."
                )

    # Question about log2 fold change / p-value / FDR
    if "fold change" in query_lower or "p-value" in query_lower or "fdr" in query_lower or "volcano" in query_lower:
        fdr_th = de.get('fdr_threshold', 0.05)
        fc_th = de.get('log2fc_threshold', 1.0)
        n_up = de.get('upregulated', 0)
        n_down = de.get('downregulated', 0)
        return (
            f"### 📊 Explanation of Differential Expression Statistics in BioPrism\n\n"
            f"In your analysis of **{dataset_name}**, genes were evaluated using two primary statistical axes:\n\n"
            f"1. **log₂ Fold Change (Effect Size):**\n"
            f"   - Measures the magnitude and direction of expression shift between `{de.get('case_group')}` and `{de.get('control_group')}`.\n"
            f"   - A $\\log_2\\text{{FC}} = +1.0$ corresponds to a 2-fold upregulation in the case group.\n"
            f"   - A $\\log_2\\text{{FC}} = -1.0$ corresponds to a 2-fold downregulation (50% reduction).\n"
            f"   - Current cutoff applied: `|log2FC| >= {fc_th}`.\n\n"
            f"2. **Adjusted p-value (False Discovery Rate / FDR):**\n"
            f"   - Because testing thousands of transcripts simultaneously leads to massive false positives (*multiple hypothesis testing problem*), BioPrism applies the Benjamini–Hochberg procedure.\n"
            f"   - Current threshold: `FDR < {fdr_th}` (less than {int(fdr_th*100)}% expected false discoveries among flagged genes).\n\n"
            f"3. **Your Results:**\n"
            f"   - **{n_up}** genes met the criteria for significant upregulation.\n"
            f"   - **{n_down}** genes met the criteria for significant downregulation."
        )

    # Question about pathways
    if "pathway" in query_lower or "kegg" in query_lower or "go" in query_lower or "enrich" in query_lower:
        if not pathways:
            return "No enriched biological pathways are currently loaded in your session. Please execute Module 2 (Pathway Enrichment) first."
        path_list_md = "\n".join([f"- **{p.get('term')}** (FDR: `{p.get('adj_p_val'):.2e}`, Odds Ratio: `{p.get('odds_ratio')}`, Overlapping Genes: `{p.get('genes')}`)" for p in pathways[:6]])
        return (
            f"### 🔬 Functional Pathway Enrichment Summary\n\n"
            f"Pathway analysis maps your statistically significant differentially expressed genes to curated biological databases (KEGG / GO / Reactome):\n\n"
            f"{path_list_md}\n\n"
            f"**Biological Significance:** The convergence of multiple candidate genes into these signaling pathways indicates that the observed expression shifts reflect coordinate cellular pathway dysregulation rather than isolated transcriptional noise."
        )

    # Question about hub genes / network
    if "hub" in query_lower or "network" in query_lower or "centrality" in query_lower or "string" in query_lower:
        if not hubs:
            return "Protein interaction network analysis has not been executed yet. Run Module 3 to construct the STRING PPI network and identify hub genes."
        top_hubs_md = "\n".join([f"{h.get('Hub_Rank')}. **{h.get('Gene')}** — Score: `{h.get('Composite_Hub_Score'):.4f}` (Degree: `{h.get('Degree')}`, Betweenness: `{h.get('Betweenness_Centrality'):.4f}`, log2FC: `{h.get('log2FoldChange'):.2f}`)" for h in hubs[:5]])
        return (
            f"### 🕸️ Protein-Protein Interaction Network & Hub Gene Rankings\n\n"
            f"Proteins rarely operate in isolation. In Module 3, significant genes were queried against the **STRING Database** and topologically evaluated via **NetworkX**:\n\n"
            f"**BioPrism Composite Hub Score Formula:**\n"
            f"$$\\text{{Hub Score}} = 0.40 \\times \\frac{{\\text{{Degree}}}}{{\\text{{Max Degree}}}} + 0.35 \\times C_{{\\text{{Betweenness}}}} + 0.25 \\times \\text{{Significance Score}}$$\n\n"
            f"**Top Candidate Hub Biomarkers Identified:**\n"
            f"{top_hubs_md}\n\n"
            f"High betweenness centrality indicates a 'bottleneck' protein that controls communication between different functional modules of the interactome."
        )

    # Question about classification / ML / validation
    if "model" in query_lower or "auc" in query_lower or "roc" in query_lower or "classif" in query_lower or "valid" in query_lower or "accuracy" in query_lower:
        if not ml or not ml.get('success'):
            return "Machine learning validation has not yet been executed in this session. Run Module 4 to train and validate classifiers."
        m = ml.get('metrics', {})
        return (
            f"### 🎯 Machine Learning Biomarker Validation Summary\n\n"
            f"The candidate biomarker panel was evaluated using **{m.get('model_name')}** under **{ml.get('cv_method')}**:\n\n"
            f"- **ROC-AUC:** `{m.get('roc_auc', 0.0):.3f}` (Discriminative power across all classification thresholds)\n"
            f"- **Cross-Validation Accuracy:** `{m.get('accuracy', 0.0):.1%}`\n"
            f"- **Precision:** `{m.get('precision', 0.0):.1%}`\n"
            f"- **Recall (Sensitivity):** `{m.get('recall', 0.0):.1%}`\n"
            f"- **F1-Score:** `{m.get('f1_score', 0.0):.3f}`\n\n"
            f"**Data Leakage Prevention:** Feature standardization (`StandardScaler`) was strictly fit only on the training folds inside `sklearn.pipeline.Pipeline`, guaranteeing that no test sample information leaked into training."
        )

    # General overview summary
    hubs_preview = ", ".join([h.get('Gene') for h in hubs[:4]]) if hubs else "None identified yet"
    return (
        f"### 📋 BioPrism Analysis Synthesis for **{dataset_name}**\n\n"
        f"**1. Preprocessing:** {preproc.get('num_samples', 0)} samples analyzed across {preproc.get('num_genes', 0)} genes (`{preproc.get('inferred_type', 'transcriptomic data')}`).\n\n"
        f"**2. Differential Expression:** Detected **{de.get('significant_total', 0)} differentially expressed genes** ({de.get('upregulated', 0)} upregulated, {de.get('downregulated', 0)} downregulated in `{de.get('case_group', 'Case')}` vs `{de.get('control_group', 'Control')}`).\n\n"
        f"**3. Enriched Pathways:** Top pathways include *{pathways[0].get('term') if pathways else 'N/A'}*.\n\n"
        f"**4. Topological Hub Genes:** Key interactome drivers identified: **{hubs_preview}**.\n\n"
        f"**5. Machine Learning Discriminative Power:** Candidate hub panel achieved a cross-validation **ROC-AUC of {ml.get('metrics', {}).get('roc_auc', 0.0):.3f}**.\n\n"
        f"*Notice: These computational findings serve as prioritizations for prospective bench validation.*"
    )

def ask_ai_assistant(user_message, session_data):
    """
    Main entry point for AI-assisted interpretation.
    If GEMINI_API_KEY is available, queries Google Gemini API with grounded context.
    Otherwise, gracefully employs the BioPrism Deterministic Scientific Interpreter.
    """
    api_key = get_gemini_api_key()
    context_str = build_session_context_summary(session_data)

    if not api_key:
        response_text = generate_deterministic_scientific_response(user_message, session_data)
        return {
            'success': True,
            'response': response_text,
            'engine': 'BioPrism Deterministic Scientific Interpreter',
            'api_configured': False,
            'notice': 'Gemini API key not configured in environment variables. Providing domain-specific scientific interpretation directly from active session data.'
        }

    # Call Gemini API
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={api_key}"
    
    system_instruction = (
        "You are the BioPrism Bioinformatics AI Assistant, an expert in transcriptomics, differential expression (Limma/DESeq2), "
        "pathway enrichment (GSEApy/Enrichr), protein-protein interaction networks (STRING DB/NetworkX), and biomarker machine learning validation.\n\n"
        "RULES:\n"
        "1. Strictly use the provided active BioPrism session results below to answer the user's inquiry.\n"
        "2. Do NOT invent or hallucinate genes, p-values, fold changes, or pathway terms that are not in the context.\n"
        "3. If the user asks about an analysis step that has not been executed yet, clearly state that it is not available in the current session.\n"
        "4. Always maintain scientific rigor: emphasize that computational biomarker validation is hypothesis-generating and does not constitute clinical diagnostic proof without prospective clinical trials.\n"
        "5. Format your answers clearly using Markdown headers, bullet points, and LaTeX notation where appropriate ($log_2FC$)."
    )

    prompt = f"ACTIVE BIOPRISM SESSION CONTEXT:\n{context_str}\n\nUSER QUESTION:\n{user_message}"

    payload = {
        "contents": [
            {
                "parts": [
                    {"text": system_instruction + "\n\n" + prompt}
                ]
            }
        ],
        "generationConfig": {
            "temperature": 0.2,
            "maxOutputTokens": 1024
        }
    }

    try:
        req = urllib.request.Request(
            url,
            data=json.dumps(payload).encode('utf-8'),
            headers={'Content-Type': 'application/json'}
        )
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            candidates = data.get('candidates', [])
            if candidates and 'content' in candidates[0] and 'parts' in candidates[0]['content']:
                reply = candidates[0]['content']['parts'][0].get('text', '')
                return {
                    'success': True,
                    'response': reply,
                    'engine': 'Google Gemini 2.5 Flash',
                    'api_configured': True
                }
            else:
                fallback_reply = generate_deterministic_scientific_response(user_message, session_data)
                return {
                    'success': True,
                    'response': fallback_reply,
                    'engine': 'BioPrism Deterministic Scientific Interpreter (API parse fallback)',
                    'api_configured': True
                }
    except Exception as e:
        # Graceful fallback on network or API error
        fallback_reply = generate_deterministic_scientific_response(user_message, session_data)
        return {
            'success': True,
            'response': fallback_reply,
            'engine': 'BioPrism Deterministic Scientific Interpreter (API connection fallback)',
            'api_configured': True,
            'notice': f"External LLM API call encountered: {str(e)[:60]}. Provided grounded interpretation from local analysis data."
        }
