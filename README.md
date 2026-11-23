# 🧬 BioPrism: An Integrated Web-Based Pipeline for AI-Assisted Biomarker Discovery and Validation

[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)
[![Flask 3.x](https://img.shields.io/badge/Web_Framework-Flask_3.x-green.svg)](https://flask.palletsprojects.com/)
[![R 4.4.2](https://img.shields.io/badge/R-4.4.2-276DC3.svg)](https://www.r-project.org/)
[![Bioconductor Limma & DESeq2](https://img.shields.io/badge/Bioconductor-limma_&_DESeq2-brightgreen.svg)](https://bioconductor.org/)
[![Plotly](https://img.shields.io/badge/Visualization-Plotly.js-orange.svg)](https://plotly.com/)
[![Gemini AI](https://img.shields.io/badge/AI_Engine-Google_Gemini-teal.svg)](https://ai.google.dev/)
[![SQLite](https://img.shields.io/badge/Database-SQLite-blue.svg)](https://www.sqlite.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A complete, functional, modern bioinformatics web platform named **BioPrism**.
BioPrism unifies high-throughput transcriptomics preprocessing, differential expression modeling, biological pathway enrichment, protein-protein interaction networks, hub gene topological ranking, rigorous machine learning classification, and grounded AI-assisted interpretation into a seamless, modern application.

---

## 📌 Table of Contents
1. [Project Overview](#project-overview)
2. [Technology Stack](#technology-stack)
3. [Application Design & Sidebar Navigation](#application-design--sidebar-navigation)
4. [Core Analysis Capabilities](#core-analysis-capabilities)
5. [System Architecture & Directory Structure](#system-architecture--directory-structure)
6. [Installation & Execution on Windows](#installation--execution-on-windows)
7. [API Configuration & External Services](#api-configuration--external-services)
8. [Scientific & Methodological Highlights](#scientific--methodological-highlights)
9. [Automated Test Suite](#automated-test-suite)
10. [Important Scientific Disclaimers](#important-scientific-disclaimers)

---

## 🔬 1. Project Overview

High-throughput transcriptomics (microarrays, bulk RNA-seq) frequently faces the *"curse of dimensionality"* ($p \gg n$), with tens of thousands of genes measured across small patient cohorts. Discovering actionable, robust, and reproducible biomarkers requires integrating multiple levels of biological evidence:

1. **Statistical Effect Size & Significance** (Differential Expression: Limma / DESeq2 / Welch's t-test with Benjamini-Hochberg FDR)
2. **Higher-Order Cellular Context** (Biological Pathway Enrichment: GSEApy / Enrichr for KEGG, GO, Reactome, MSigDB Hallmark)
3. **Topological Indispensability** (Protein-Protein Interaction Networks: STRING DB API & NetworkX Centralities)
4. **Generalizable Predictive Performance** (Machine Learning Validation: Logistic Regression, Random Forest, SVM with strict data leakage safeguards)
5. **Grounded AI Interpretation** (LLM Assistant grounded directly in actual computed statistics without hallucination)
6. **Publication-Ready Reporting** (Interactive Plotly visualizations, downloadable CSV matrices, and printable PDF reports)

BioPrism unifies these capabilities into an elegant, responsive web platform designed to eliminate data leakage and avoid over-interpretation.

---

## 💻 2. Technology Stack

- **Frontend & UI:** HTML5, CSS3, JavaScript (ES6+), Bootstrap 5.3, Bootstrap Icons, Plotly.js 2.35, Marked.js
- **Backend & Web Server:** Python 3.13, Flask 3.1
- **Data Science & Bioinformatics:**
  - `Pandas` (v3.1) and `NumPy` (v2.5) for matrix transformations and QC
  - `SciPy` (v1.18) for statistical tests and Benjamini-Hochberg FDR monotonic adjustment
  - `R 4.4.2` with Bioconductor `limma` (microarrays/normalized log2) and `DESeq2` (RNA-seq integer counts), with automatic Python fallback
  - `GSEApy` (v1.3) with Enrichr API and curated offline pathway database
  - `STRING Database REST API` (Taxonomy 9606 Homo sapiens)
  - `NetworkX` (v3.7) for degree, betweenness, and closeness centrality computation
  - `Scikit-learn` (v1.9) with `StandardScaler` strictly encapsulated inside `Pipeline` with `StratifiedKFold`
- **Database:** SQLite (`bioprism.db`) for run provenance, session metrics, chat history, and settings
- **Artificial Intelligence:**
  - Google Gemini API (`gemini-2.5-flash`) via secure backend environment variables (`GEMINI_API_KEY`)
  - Grounded prompt injection preventing biological hallucinations
  - Deterministic Scientific Interpreter fallback when offline or no API key is set

---

## 🎨 3. Application Design & Sidebar Navigation

BioPrism features a modern bioinformatics aesthetic:
- **Color Palette:** Clean white and light backgrounds (`#f8fafc`, `#ffffff`) with teal (`#0d9488`, `#0f766e`), navy blue (`#0f172a`, `#1e293b`), and cyan accents (`#06b6d4`, `#ecfeff`).
- **Logo:** Custom SVG badge combining a DNA double helix, molecular interactome network nodes, and biomarker prism motif.
- **Responsive Layout:** Desktop sidebar with mobile/tablet collapsible drawer.
- **Sidebar Navigation:**
  - **Dashboard** (`/` or `/dashboard`): Dynamic overview KPI cards, modality breakdown donut chart, performance trend chart, and recent analysis runs.
  - **Data Upload** (`/upload`): Custom file upload (CSV/TSV), downloadable templates, demo loaders, and public GEO cohort access.
  - **Data Preprocessing** (`/preprocessing`): Quality control summary, imputation options, 2D sample clustering PCA preview, and expression matrix table preview.
  - **Differential Expression** (`/differential-expression`): Group selection, Limma/DESeq2 execution, interactive Volcano plot, Clustered Heatmap, and PCA plot.
  - **Pathway Enrichment** (`/pathway-enrichment`): KEGG/GO/Reactome enrichment, horizontal -log10(FDR) bar plot, and downloadable pathway tables.
  - **Network Analysis** (`/network-analysis`): STRING PPI integration, NetworkX centralities, Composite Hub Score ranking, and 2D force-directed Plotly network graph.
  - **Biomarker Validation** (`/biomarker-validation`): Regularized ML classification (Logistic Regression, Random Forest, SVM), ROC curves, Confusion Matrix, and independent cohort testing.
  - **AI Assistant** (`/ai-assistant`): Context-grounded chatbot interface with preset inquiry buttons and zero-hallucination verification.
  - **Reports & Exports** (`/reports`): Standalone HTML report preview, print/save to PDF, and CSV downloads for all pipeline outputs.
  - **Analysis History** (`/history`): Searchable database of past analysis runs with one-click reloading into the workspace.
  - **Settings** (`/settings`): Runtime environment diagnostics, API key status, default thresholds, and SQLite maintenance.
  - **Pipeline Workspace** (`/pipeline`): Unified 4-step continuous wizard for rapid analysis.

---

## 📦 4. Core Analysis Capabilities

### Data Input and Preprocessing
- Accepts gene expression matrices (CSV/TSV) and matched sample metadata.
- Pre-packaged standardized templates (`template_expression.csv`, `template_metadata.csv`).
- Comprehensive quality control:
  - Collapses duplicate gene symbols by arithmetic mean.
  - Removes missing rows or applies median/mean row-wise imputation.
  - Removes uninformative zero-variance genes.
  - Verifies sample-to-column parity and minimum replicates ($n \ge 2$ per group).
  - Auto-infers data type: Microarray Log2 vs. RNA-seq Raw Counts vs. Linear Normalized.
  - Generates 2D quality control PCA clustering plot.
- Supports public benchmark Gene Expression Omnibus (GEO) cohorts (e.g. GSE42568, GSE19804, GSE3325).

### Differential Expression and Pathway Enrichment
- Flexible pairwise group comparison (e.g. Tumor vs. Normal, Case vs. Control).
- Statistical workflow:
  - Native execution of R `limma` (microarrays/normalized log2) via `r_scripts/run_limma.R`
  - Native execution of R `DESeq2` (RNA-seq integer counts) via `r_scripts/run_deseq2.R`
  - Statistical fallback: Python Welch's t-test with monotonic Benjamini-Hochberg FDR calculation.
- Interactive Plotly Visualizations:
  - **Volcano Plot:** Hover tooltips displaying gene symbol, log2FC, -log10(FDR), and status.
  - **Clustered Expression Heatmap:** Standardized Z-score expression of top candidate genes across samples.
  - **PCA Separation Plot:** Principal Component Analysis visualizing cohort separation.
- Biological Pathway Enrichment:
  - GSEApy Enrichr integration for `KEGG_2021_Human`, `GO_Biological_Process_2023`, `Reactome_2022`, `MSigDB_Hallmark_2020`.
  - Curated offline pathway fallback database for air-gapped or network-restricted environments.
  - Interactive horizontal bar chart sorted by $-\log_{10}(\text{FDR})$ with odds ratio color scale.

### Network Analysis and Hub Gene Identification
- Submits significant genes to the **STRING Database REST API** (Species 9606 Homo sapiens).
- Configurable confidence score cutoff (Medium: 400, High: 700, Highest: 900).
- Distinguishes experimentally validated interactions from predicted functional associations.
- Computes graph properties using **NetworkX**: Degree, Degree Centrality, Betweenness Centrality, Closeness Centrality.
- Computes the **Composite Hub Score**:
  $$\text{Hub Score} = 0.40 \times \frac{\text{Degree}}{\text{Max Degree}} + 0.35 \times C_{\text{Betweenness}} + 0.25 \times \min\left(1.0, \frac{|\log_2\text{FC}| \cdot (-\log_{10}\text{FDR})}{10}\right)$$
- Interactive 2D force-directed Plotly network visualization with node sizing by degree and color coding by regulation status (Teal=Up, Blue=Down, Gold=Top Hub).
- Downloadable ranked hub genes table and network edge list in CSV format.

### Biomarker Validation and Classification
- Prioritizes candidate biomarkers based on differential expression and network topology.
- Machine Learning Classifiers:
  - **Logistic Regression** (L2 regularized)
  - **Random Forest Classifier** (100 estimators)
  - **Support Vector Machine (SVM)** (Linear kernel)
  - **Model Comparison Mode** comparing all three algorithms head-to-head.
- **Strict Data Leakage Prevention:** Feature scaling (`StandardScaler`) is strictly instantiated and fit inside training folds within `sklearn.pipeline.Pipeline`.
- Cross-Validation: Stratified $k$-Fold ($k=3, 5, 10$) or Leave-One-Out (for small cohorts).
- Metrics: Accuracy, Precision, Recall, F1-Score, and ROC-AUC.
- Interactive Plotly **ROC Curves** and **Confusion Matrix Heatmap**.
- Independent cohort testing with candidate gene expression distribution boxplots across cohorts.

### AI-Assisted Result Interpretation
- Dedicated conversational chatbot interface.
- Grounded context injection: passes the exact active session results (genes, log2FC, FDR, pathways, hub scores, classifier AUC) into the prompt context.
- Zero-hallucination guarantee: instructed never to invent genes or statistics.
- Pre-configured prompt shortcuts:
  - "Summarize Key Findings"
  - "Explain Why Top Hub Genes Were Selected"
  - "Explain log2 Fold Change & FDR Statistics"
  - "Interpret Biological Pathway Enrichments"
  - "Evaluate Classifier Reliability & Leakage Safeguards"
  - "Draft Report Abstract"
- Google Gemini API integration with automatic fallback to BioPrism's Deterministic Scientific Interpreter when no API key is configured.

### Visualization and Report Generation
- Standalone, publication-ready HTML Capstone Report including all 9 formal sections:
  1. Dataset Information
  2. Preprocessing & QC Details
  3. Differential Expression Results
  4. Biological Pathway Enrichment Results
  5. Protein-Protein Interaction Network Analysis
  6. Candidate Biomarker Hub Prioritization
  7. Machine Learning Classification & Cross-Validation Metrics
  8. AI-Assisted Scientific Interpretation
  9. Methodological Parameters & Scientific Limitations
- Print-ready `@media print` CSS for browser "Print / Save as PDF".
- Direct download of CSV matrices for each stage.

---

## 🏛️ 5. System Architecture & Directory Structure

```
BioPrism/
│
├── app.py                     # Main Flask web application & API routing
├── config.py                  # Environment settings, paths, and configurations
├── requirements.txt           # Python package dependencies
├── .env.example               # Example environment variables (GEMINI_API_KEY)
├── bioprism.db                # SQLite database for history, runs, and chats
├── generate_sample_data.py    # Generator for verified synthetic benchmark cohorts
├── README.md                  # Comprehensive documentation and setup guide
│
├── modules/                   # Core modular bioinformatics & ML engines
│   ├── __init__.py
│   ├── data_upload.py         # File ingestion & GEO cohort loader
│   ├── preprocessing.py       # Quality control, imputation, normalization, PCA
│   ├── differential_expression.py # Limma, DESeq2, Welch's t-test, Volcano, Heatmap
│   ├── pathway_enrichment.py  # GSEApy Enrichr biological pathway enrichment
│   ├── network_analysis.py    # STRING API client, NetworkX graph theory, Hubs
│   ├── biomarker_validation.py# ML validation pipeline & independent cohort testing
│   ├── classification.py      # Logistic Regression, Random Forest, SVM models
│   ├── ai_assistant.py        # Gemini API client & grounded scientific interpreter
│   └── report_generator.py    # Standalone HTML report generator with 9 sections
│
├── database/                  # SQLite database management
│   ├── __init__.py
│   └── db.py                  # Schema creation, run logging, dashboard stats, chat
│
├── pipeline/                  # Backward compatibility re-export layer for tests
│   ├── preprocessor.py
│   ├── de_analysis.py
│   ├── network_analysis.py
│   ├── validation_ml.py
│   └── report_generator.py
│
├── templates/                 # Modern responsive HTML5 templates
│   ├── base.html              # Core layout with sidebar navigation & topbar
│   ├── dashboard.html         # Main dashboard with live SQLite statistics & charts
│   ├── upload.html            # File upload & GEO dataset importer
│   ├── preprocessing.html     # QC inspection & PCA preview
│   ├── differential_expression.html # DEGs, Volcano, Heatmap, PCA
│   ├── pathway_enrichment.html# GSEApy pathways & bar chart
│   ├── network_analysis.html  # STRING interactome & Hub gene table
│   ├── validation.html        # ML validation, ROC, Confusion Matrix
│   ├── ai_assistant.html      # Grounded conversational AI assistant
│   ├── reports.html           # Live report preview & PDF/CSV exports
│   ├── history.html           # SQLite analysis provenance & session reloader
│   ├── settings.html          # Runtime diagnostics & default thresholds
│   ├── pipeline.html          # Continuous 4-step wizard workspace
│   └── help.html              # Capstone defense questions & methodology guide
│
├── static/
│   ├── css/
│   │   └── style.css          # Bioinformatics theme (teal, navy, cyan, white)
│   ├── js/
│   │   ├── main.js            # Core frontend utilities & session sync
│   │   └── app.js             # Interactive pipeline stepper logic
│   ├── images/
│   │   └── logo.svg           # Custom DNA helix + molecular network SVG logo
│   └── sample_data/           # Standardized templates & demo benchmark cohorts
│       ├── template_expression.csv
│       ├── template_metadata.csv
│       ├── demo_microarray_expression.csv
│       ├── demo_microarray_metadata.csv
│       ├── demo_rnaseq_counts.csv
│       ├── demo_rnaseq_metadata.csv
│       └── independent_validation_data.csv
│
├── r_scripts/                 # Standalone Bioconductor R scripts
│   ├── run_limma.R            # Limma linear modeling & empirical Bayes
│   └── run_deseq2.R           # DESeq2 negative binomial GLM for RNA-seq counts
│
├── uploads/                   # Secure temporary storage for user uploads
├── results/                   # Per-session analysis outputs & generated reports
└── tests/                     # Comprehensive automated pytest suite
    ├── conftest.py
    ├── test_preprocessing.py
    ├── test_differential_expression.py
    ├── test_network_analysis.py
    ├── test_validation_ml.py
    ├── test_ai_assistant.py
    ├── test_database_and_geo.py
    ├── test_report_and_routes.py
    └── test_end_to_end.py
```

---

## 🚀 6. Installation & Execution on Windows

### Prerequisites
- **Python 3.10+** (Tested on Python 3.13)
- **R 4.4+** with `limma` and `DESeq2` (Optional; native Python fallback included)

### Step 1: Open PowerShell or Command Prompt
```powershell
cd c:\Users\rajas\OneDrive\Dokumen\MedBio_BioPrism
```

### Step 2: Activate Python Virtual Environment
```powershell
.\.venv\Scripts\Activate.ps1
```
*(Or create a new environment with `python -m venv .venv` if setting up for the first time)*

### Step 3: Install Required Packages
```powershell
.\.venv\Scripts\pip install -r requirements.txt
```

### Step 4: Launch the BioPrism Web Application
```powershell
.\.venv\Scripts\python app.py
```

Open your web browser and navigate to:
**`http://127.0.0.1:5000`**

---

## 🔑 7. API Configuration & External Services

### Google Gemini AI API (Optional)
To enable live LLM interpretations via Google Gemini:
1. Obtain an API key from [Google AI Studio](https://aistudio.google.com/).
2. Set the environment variable in your terminal before launching:
   ```powershell
   $env:GEMINI_API_KEY="your_api_key_here"
   .\.venv\Scripts\python app.py
   ```
   *(Or add `GEMINI_API_KEY=your_key` to a `.env` file)*
3. **Graceful Fallback:** If no API key is configured or if the external API is unreachable, BioPrism automatically engages its **Deterministic Scientific Interpreter**, delivering rigorous answers calculated directly from your session data with zero hallucinations. API keys are never exposed in client-side code.

### STRING Database API
BioPrism queries the official STRING REST API (`https://string-db.org/api/json/network`). If network access is offline or times out, the application seamlessly activates its local curated reference interactome database.

### R & Bioconductor Limma/DESeq2
BioPrism automatically checks for `Rscript.exe` in `C:\Program Files\R\R-4.4.2\bin\Rscript.exe` or standard system `PATH`. When available, it executes native Bioconductor `limma` and `DESeq2`. If R is not installed, it automatically engages the internal Python Welch's t-test with Benjamini-Hochberg FDR correction.

---

## 🧪 8. Scientific & Methodological Highlights

1. **Limma Empirical Bayes Shrinkage:** Moderates gene-wise sample variances by borrowing information across all transcripts, preventing false positives caused by spuriously small sample variances in small cohorts.
2. **Heteroscedasticity Safeguard:** Raw RNA-seq count matrices are not passed into homoscedastic linear models; the preprocessor auto-detects integer count matrices and guides users toward DESeq2 negative binomial modeling.
3. **Data Leakage Safeguard:** Feature scaling (`StandardScaler`) is fit strictly on training folds inside `sklearn.pipeline.Pipeline`. Test fold data is never exposed during feature scaling or model fitting.
4. **False Discovery Rate (FDR) Control:** Applies the Benjamini-Hochberg step-up procedure with enforced right-to-left cumulative minimum monotonicity to control the proportion of false positives among rejected null hypotheses.
5. **Functional Association vs Physical Binding:** Clearly highlights that STRING interactions encompass predicted functional links (co-expression, text mining) alongside experimentally established physical complexes.

---

## 🧪 9. Automated Test Suite

To run the complete automated test suite verifying all modules, APIs, database operations, and end-to-end flows:

```powershell
.\.venv\Scripts\pytest -v
```

All test modules run automatically:
- `test_preprocessing.py`: Validation, missing value imputation, duplicate collapsing, RNA-seq count inference.
- `test_differential_expression.py`: Limma/t-test execution, Benjamini-Hochberg monotonicity, GSEApy pathway enrichment.
- `test_network_analysis.py`: STRING API integration, local PPI fallback, NetworkX centralities, hub ranking.
- `test_validation_ml.py`: Logistic Regression, Random Forest, SVM, model comparison mode, LOOCV, data leakage safeguards.
- `test_ai_assistant.py`: Grounded context compilation, deterministic scientific responses, prompt safety.
- `test_database_and_geo.py`: SQLite schema, record logging, dashboard dynamic statistics, GEO cohort loader.
- `test_report_and_routes.py`: Standalone HTML report generation, all Flask web routes, download endpoints.
- `test_end_to_end.py`: End-to-end integration test spanning the complete 5-phase pipeline.

---

## ⚠️ 10. Important Scientific Disclaimers

> **Exploratory & Retrospective Nature:**
> Machine learning classification performance (e.g., cross-validated ROC-AUC) demonstrates mathematical discriminative capacity within retrospective transcriptomics cohorts. Computational biomarker prioritization serves as an exploratory hypothesis-generating step. It does **not** constitute standalone clinical diagnostic validity, disease prognosis, or therapeutic indication without prospective blinded multi-center clinical validation and orthogonal experimental assays (such as RT-qPCR, Western blot, or immunohistochemistry).

---

## 📄 License
This project is licensed under the MIT License.
