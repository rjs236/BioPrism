// BioPrism Interactive Pipeline Engine

let currentSessionId = (typeof BioPrism !== 'undefined' && BioPrism.getSessionId)
    ? BioPrism.getSessionId()
    : (document.getElementById('display-session-id') ? document.getElementById('display-session-id').innerText.trim() : null);

if (typeof BioPrism !== 'undefined' && currentSessionId) {
    BioPrism.setSessionId(currentSessionId);
    const sidEl = document.getElementById('display-session-id');
    if (sidEl) sidEl.innerText = currentSessionId;
}

// Helper: Show/Hide Spinners
function showSpinner(msg = "Processing...") {
    const sp = document.getElementById('global-spinner');
    const msgEl = document.getElementById('spinner-message');
    if (msgEl) msgEl.innerText = msg;
    if (sp) sp.style.display = 'block';
}

function hideSpinner() {
    const sp = document.getElementById('global-spinner');
    if (sp) sp.style.display = 'none';
}

// Stepper Navigation
function switchStep(step) {
    for (let i = 1; i <= 4; i++) {
        const btn = document.getElementById(`btn-step-${i}`);
        const content = document.getElementById(`step-${i}-content`);
        if (btn) {
            btn.classList.remove('active');
            if (i < step) btn.classList.add('completed');
        }
        if (content) content.style.display = 'none';
    }

    const activeBtn = document.getElementById(`btn-step-${step}`);
    const activeContent = document.getElementById(`step-${step}-content`);
    if (activeBtn) activeBtn.classList.add('active');
    if (activeContent) activeContent.style.display = 'block';

    // Resize Plotly plots if visible
    if (step === 2 && document.getElementById('volcano-plot-container')) {
        Plotly.Plots.resize('volcano-plot-container');
    } else if (step === 3 && document.getElementById('network-plot-container')) {
        Plotly.Plots.resize('network-plot-container');
    } else if (step === 4) {
        if (document.getElementById('roc-plot-container')) Plotly.Plots.resize('roc-plot-container');
        if (document.getElementById('cm-plot-container')) Plotly.Plots.resize('cm-plot-container');
    }
}

// ==================== MODULE 1: UPLOAD & QC ====================

function handleUploadResponse(data) {
    hideSpinner();
    const alertBox = document.getElementById('qc-alerts');
    alertBox.innerHTML = '';

    if (!data.success) {
        let errHtml = '<div class="alert alert-danger py-2 small mb-3"><strong>Validation Errors:</strong><ul class="mb-0 ps-3">';
        (data.errors || []).forEach(err => { errHtml += `<li>${err}</li>`; });
        errHtml += '</ul></div>';
        alertBox.innerHTML = errHtml;
        return;
    }

    // Warnings if any
    if (data.warnings && data.warnings.length > 0) {
        let warnHtml = '<div class="alert alert-warning py-2 small mb-3"><strong>QC Adjustments / Warnings:</strong><ul class="mb-0 ps-3">';
        data.warnings.forEach(w => { warnHtml += `<li>${w}</li>`; });
        warnHtml += '</ul></div>';
        alertBox.innerHTML = warnHtml;
    } else {
        alertBox.innerHTML = '<div class="alert alert-success py-2 small mb-3"><i class="bi bi-check-circle-fill me-1"></i> Dataset successfully validated and aligned!</div>';
    }

    // Update KPIs
    const info = data.info || {};
    document.getElementById('kpi-num-genes').innerText = info.num_genes || '-';
    document.getElementById('kpi-num-samples').innerText = info.num_samples || '-';
    document.getElementById('kpi-inferred-type').innerText = info.inferred_type === 'rnaseq_counts' ? 'RNA-seq Counts' : 'Microarray Log2';
    document.getElementById('kpi-range').innerText = `[${info.min_expression}, ${info.max_expression}]`;

    // Cohort Badges
    const cohortDiv = document.getElementById('kpi-cohort-badges');
    cohortDiv.innerHTML = '';
    const groups = info.groups || {};
    for (const [grp, count] of Object.entries(groups)) {
        cohortDiv.innerHTML += `<span class="badge bg-primary me-2">${grp}: ${count} samples</span>`;
    }

    // Expression Matrix Preview Table
    const exprHeader = document.getElementById('expr-preview-header');
    exprHeader.innerHTML = '';
    (info.preview_columns || []).forEach(col => {
        exprHeader.innerHTML += `<th>${col}</th>`;
    });

    const exprBody = document.getElementById('expr-preview-body');
    exprBody.innerHTML = '';
    (data.expr_preview || []).forEach(row => {
        let tr = '<tr>';
        (info.preview_columns || []).forEach(col => {
            tr += `<td>${row[col] !== undefined ? row[col] : '-'}</td>`;
        });
        tr += '</tr>';
        exprBody.innerHTML += tr;
    });

    // Metadata Preview Table
    const metaBody = document.getElementById('meta-preview-body');
    metaBody.innerHTML = '';
    (data.meta_preview || []).forEach(row => {
        metaBody.innerHTML += `<tr><td><strong>${row.Sample}</strong></td><td><span class="badge bg-secondary">${row.Group}</span></td></tr>`;
    });

    // Populate Group Selectors for Module 2
    const caseSelect = document.getElementById('de-case-group');
    const ctrlSelect = document.getElementById('de-control-group');
    if (caseSelect && ctrlSelect && info.available_groups) {
        caseSelect.innerHTML = '';
        ctrlSelect.innerHTML = '';
        info.available_groups.forEach((g, idx) => {
            const opt1 = new Option(g, g, false, idx === 0);
            const opt2 = new Option(g, g, false, idx === 1 || idx === info.available_groups.length - 1);
            caseSelect.add(opt1);
            ctrlSelect.add(opt2);
        });
    }

    // Method selection suggestion based on format
    const methodSelect = document.getElementById('de-method-select');
    if (methodSelect) {
        if (info.inferred_type === 'rnaseq_counts') {
            methodSelect.value = 'deseq2';
        } else {
            methodSelect.value = 'limma';
        }
    }

    // Toggle container
    document.getElementById('qc-placeholder').style.display = 'none';
    document.getElementById('qc-results').style.display = 'block';

    if (data.session_id) {
        currentSessionId = data.session_id;
        document.getElementById('display-session-id').innerText = data.session_id;
        if (typeof BioPrism !== 'undefined') BioPrism.setSessionId(data.session_id);
    }
    if (data.session_snapshot && typeof BioPrism !== 'undefined') {
        BioPrism.saveSessionSnapshot(data.session_id || currentSessionId, data.session_snapshot, true);
    }
    if (data.history_record && typeof BioPrism !== 'undefined') {
        BioPrism.saveHistoryRecord(data.history_record);
    }
    if (data.dataset_name) {
        document.getElementById('display-dataset-name').innerText = data.dataset_name;
    }
}

function _cacheWorkspaceState(data) {
    if (typeof BioPrism === 'undefined' || !data) return;
    if (data.session_snapshot) {
        BioPrism.saveSessionSnapshot(data.session_id || currentSessionId, data.session_snapshot, Boolean(data.reset_snapshot));
    }
    if (data.history_record) {
        BioPrism.saveHistoryRecord(data.history_record);
    }
}

async function _pipelinePost(url, payload) {
    if (typeof BioPrism !== 'undefined' && BioPrism.apiPost) {
        return await BioPrism.apiPost(url, payload);
    }
    const r = await fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
    });
    return await r.json();
}

// 1-Click Demo Loader
function loadDemoDataset(type) {
    showSpinner(`Loading ${type} benchmark demo dataset...`);
    fetch('/api/load-demo', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ type: type, session_id: currentSessionId })
    })
    .then(r => r.json())
    .then(data => handleUploadResponse(data))
    .catch(err => {
        hideSpinner();
        alert(`Failed to load demo: ${err.message}`);
    });
}

// Form Upload Handler
const uploadForm = document.getElementById('upload-form');
if (uploadForm) {
    uploadForm.addEventListener('submit', function(e) {
        e.preventDefault();
        showSpinner("Uploading, inspecting, and validating matrix parity...");
        const formData = new FormData(uploadForm);
        formData.set('session_id', currentSessionId);

        fetch('/api/upload', {
            method: 'POST',
            body: formData
        })
        .then(r => r.json())
        .then(data => handleUploadResponse(data))
        .catch(err => {
            hideSpinner();
            alert(`Upload error: ${err.message}`);
        });
    });
}

// ==================== MODULE 2: DIFFERENTIAL EXPRESSION & GSEA ====================

function runDifferentialExpression() {
    const caseGrp = document.getElementById('de-case-group').value;
    const ctrlGrp = document.getElementById('de-control-group').value;

    if (!caseGrp || !ctrlGrp) {
        alert("Please ensure valid Case and Control groups are selected.");
        return;
    }
    if (caseGrp === ctrlGrp) {
        alert("Case and Control groups must be distinct.");
        return;
    }

    const method = document.getElementById('de-method-select').value;
    const fdr = document.getElementById('de-fdr-input').value;
    const fc = document.getElementById('de-fc-input').value;

    showSpinner(`Executing ${method.toUpperCase()} differential expression model...`);

    _pipelinePost('/api/run-de', {
        session_id: currentSessionId,
        case_group: caseGrp,
        control_group: ctrlGrp,
        method: method,
        fdr_threshold: fdr,
        log2fc_threshold: fc
    })
    .then(data => {
        hideSpinner();
        if (!data.success) {
            alert(`Differential Expression Error: ${data.message}`);
            return;
        }
        _cacheWorkspaceState(data);

        const sum = data.summary || {};
        document.getElementById('kpi-de-sig').innerText = sum.significant_total || 0;
        document.getElementById('kpi-de-up').innerText = sum.upregulated || 0;
        document.getElementById('kpi-de-down').innerText = sum.downregulated || 0;
        document.getElementById('kpi-de-ns').innerText = sum.not_significant || 0;
        document.getElementById('de-engine-badge').innerText = sum.engine || 'Differential Expression';

        // Render Volcano Plot
        if (data.volcano_plot) {
            Plotly.newPlot('volcano-plot-container', data.volcano_plot.data, data.volcano_plot.layout, { responsive: true });
        }

        // Render DE Table
        const tableBody = document.getElementById('de-genes-table-body');
        tableBody.innerHTML = '';
        (data.top_genes || []).forEach(g => {
            const badgeClass = g.Regulation === 'Upregulated' ? 'badge-reg-up' :
                               g.Regulation === 'Downregulated' ? 'badge-reg-down' : 'badge-reg-ns';
            tableBody.innerHTML += `
                <tr>
                    <td><strong>${g.Gene}</strong></td>
                    <td>${Number(g.log2FoldChange).toFixed(3)}</td>
                    <td>${Number(g.AveExpr || 0).toFixed(2)}</td>
                    <td>${Number(g.p_value).toExponential(2)}</td>
                    <td>${Number(g.adj_p_val).toExponential(2)}</td>
                    <td><span class="${badgeClass}">${g.Regulation}</span></td>
                </tr>
            `;
        });

        // Enable download link
        const dlBtn = document.getElementById('btn-download-de-csv');
        if (dlBtn) {
            dlBtn.classList.remove('disabled');
            dlBtn.href = `/api/download/de_results?session_id=${currentSessionId}`;
        }
        const dlBtn2 = document.getElementById('btn-dl-de-csv-2');
        if (dlBtn2) {
            dlBtn2.classList.remove('disabled');
            dlBtn2.href = `/api/download/de_results?session_id=${currentSessionId}`;
        }
    })
    .catch(err => {
        hideSpinner();
        alert(`DE Execution Error: ${err.message}`);
    });
}

function runPathwayEnrichment() {
    const db = document.getElementById('enrichment-db-select').value;
    showSpinner(`Querying GSEApy / Biological pathways (${db})...`);

    _pipelinePost('/api/run-enrichment', { session_id: currentSessionId, gene_set: db })
    .then(data => {
        hideSpinner();
        if (!data.success) {
            alert(`Pathway Enrichment Note: ${data.message}`);
            return;
        }
        _cacheWorkspaceState(data);

        const body = document.getElementById('pathways-table-body');
        body.innerHTML = '';
        const results = data.results || [];
        if (results.length === 0) {
            body.innerHTML = '<tr><td colspan="3" class="text-center py-2 text-muted">No enriched pathways meeting threshold.</td></tr>';
            return;
        }

        results.slice(0, 15).forEach(p => {
            body.innerHTML += `
                <tr>
                    <td><strong>${p.term}</strong></td>
                    <td>${Number(p.adj_p_val).toExponential(2)}</td>
                    <td><small class="text-truncate d-inline-block" style="max-width: 140px;" title="${p.genes}">${p.genes}</small></td>
                </tr>
            `;
        });
    })
    .catch(err => {
        hideSpinner();
        alert(`Enrichment Error: ${err.message}`);
    });
}

// ==================== MODULE 3: PPI NETWORK & HUB GENES ====================

function runNetworkAnalysis() {
    const species = document.getElementById('network-species-select').value;
    const conf = document.getElementById('network-confidence-select').value;

    showSpinner("Submitting genes to STRING DB API & computing NetworkX graph centralities...");

    _pipelinePost('/api/run-network', {
        session_id: currentSessionId,
        species: species,
        confidence_score: conf
    })
    .then(data => {
        hideSpinner();
        if (!data.success) {
            alert(`Network Construction Note: ${data.message}`);
            return;
        }
        _cacheWorkspaceState(data);

        const sum = data.summary || {};
        document.getElementById('kpi-net-nodes').innerText = sum.num_nodes || 0;
        document.getElementById('kpi-net-edges').innerText = sum.num_edges || 0;
        document.getElementById('kpi-net-density').innerText = sum.density || 0;
        document.getElementById('kpi-net-components').innerText = sum.num_components || 0;
        document.getElementById('network-status-badge').innerText = sum.source_message || 'Network Built';

        // Render Plotly Network Graph
        if (data.network_json) {
            Plotly.newPlot('network-plot-container', data.network_json.data, data.network_json.layout, { responsive: true });
        }

        // Render Hub Genes Table
        const hubBody = document.getElementById('hub-genes-table-body');
        hubBody.innerHTML = '';
        const hubs = data.hub_genes || [];
        hubs.slice(0, 15).forEach(h => {
            hubBody.innerHTML += `
                <tr>
                    <td><span class="badge bg-secondary">#${h.Hub_Rank}</span></td>
                    <td><strong>${h.Gene}</strong></td>
                    <td><span class="badge bg-primary">${h.Composite_Hub_Score}</span></td>
                    <td>${h.Degree}</td>
                    <td>${h.Betweenness_Centrality}</td>
                    <td>${Number(h.log2FoldChange).toFixed(2)}</td>
                </tr>
            `;
        });

        // Pre-populate Module 4 candidate genes input
        const topHubGeneNames = hubs.slice(0, 8).map(h => h.Gene).join(', ');
        const mlInput = document.getElementById('ml-genes-input');
        if (mlInput) mlInput.value = topHubGeneNames;

        // Enable download buttons
        const hubDl = document.getElementById('btn-download-hubs-csv');
        if (hubDl) {
            hubDl.classList.remove('disabled');
            hubDl.href = `/api/download/hub_genes?session_id=${currentSessionId}`;
        }
    })
    .catch(err => {
        hideSpinner();
        alert(`Network Analysis Error: ${err.message}`);
    });
}

// ==================== MODULE 4: VALIDATION & REPORT ====================

function runMLValidation() {
    const genesStr = document.getElementById('ml-genes-input').value;
    const genes = genesStr.split(',').map(s => s.trim().toUpperCase()).filter(s => s.length > 0);

    if (genes.length === 0) {
        alert("Please provide at least 1 candidate gene for classification.");
        return;
    }

    const clf = document.getElementById('ml-classifier-select').value;
    const splits = document.getElementById('ml-cv-select').value;
    const useIndep = document.getElementById('check-independent-val').checked;

    showSpinner("Executing Stratified Cross-Validation with strict leakage prevention...");

    _pipelinePost('/api/run-validation', {
        session_id: currentSessionId,
        selected_genes: genes,
        classifier: clf,
        cv_splits: splits,
        use_independent_cohort: useIndep
    })
    .then(data => {
        hideSpinner();
        if (!data.success) {
            alert(`Validation Error: ${data.message}`);
            return;
        }
        _cacheWorkspaceState(data);

        const m = data.metrics || {};
        document.getElementById('kpi-ml-auc').innerText = m.roc_auc !== undefined ? m.roc_auc.toFixed(3) : '-';
        document.getElementById('kpi-ml-acc').innerText = m.accuracy !== undefined ? `${(m.accuracy * 100).toFixed(1)}%` : '-';
        document.getElementById('kpi-ml-prec').innerText = m.precision !== undefined ? `${(m.precision * 100).toFixed(1)}%` : '-';
        document.getElementById('kpi-ml-rec').innerText = m.recall !== undefined ? `${(m.recall * 100).toFixed(1)}%` : '-';
        document.getElementById('kpi-ml-f1').innerText = m.f1_score !== undefined ? m.f1_score.toFixed(3) : '-';

        // Render Model Comparison Table if comparison data exists
        const compContainer = document.getElementById('model-comparison-container');
        const compBody = document.getElementById('model-comparison-table-body');
        if (m.comparison && m.comparison.length > 0) {
            compBody.innerHTML = '';
            m.comparison.forEach(c => {
                compBody.innerHTML += `
                    <tr>
                        <td><strong>${c.Model}</strong></td>
                        <td><span class="badge bg-primary">${Number(c.ROC_AUC).toFixed(3)}</span></td>
                        <td>${(Number(c.Accuracy) * 100).toFixed(1)}%</td>
                        <td>${(Number(c.Precision) * 100).toFixed(1)}%</td>
                        <td>${(Number(c.Recall) * 100).toFixed(1)}%</td>
                        <td>${Number(c.F1_Score).toFixed(3)}</td>
                    </tr>
                `;
            });
            if (compContainer) compContainer.style.display = 'block';
        } else if (compContainer) {
            compContainer.style.display = 'none';
        }

        // Render ROC Curve
        if (data.roc_plot_json) {
            Plotly.newPlot('roc-plot-container', data.roc_plot_json.data, data.roc_plot_json.layout, { responsive: true });
        }

        // Render Confusion Matrix
        if (data.cm_plot_json) {
            Plotly.newPlot('cm-plot-container', data.cm_plot_json.data, data.cm_plot_json.layout, { responsive: true });
        }

        // Render Gene Coefficients
        const featBody = document.getElementById('ml-features-table-body');
        featBody.innerHTML = '';
        (m.gene_importance || []).forEach(f => {
            featBody.innerHTML += `
                <tr>
                    <td><strong>${f.Gene}</strong></td>
                    <td>${f.Coefficient !== undefined ? f.Coefficient : f.Importance}</td>
                    <td><small class="text-muted">${f.Direction}</small></td>
                </tr>
            `;
        });

        // Independent Cohort Results
        const indepBox = document.getElementById('indep-val-box');
        const indepDetails = document.getElementById('indep-val-details');
        if (m.independent_evaluation && m.independent_evaluation.success) {
            const ie = m.independent_evaluation;
            indepDetails.innerHTML = `Tested on ${ie.num_samples} independent validation samples. Accuracy: <strong>${ie.accuracy ? (ie.accuracy * 100).toFixed(1) + '%' : 'N/A'}</strong>, ROC-AUC: <strong>${ie.roc_auc || 'N/A'}</strong>. Predicted labels: <code>${ie.predicted_classes.join(', ')}</code>`;
            indepBox.style.display = 'flex';
        } else {
            indepBox.style.display = 'none';
        }

        // Enable download link
        const dlMetrics = document.getElementById('btn-dl-metrics-csv');
        if (dlMetrics) {
            dlMetrics.classList.remove('disabled');
            dlMetrics.href = `/api/download/validation_metrics?session_id=${currentSessionId}`;
        }
    })
    .catch(err => {
        hideSpinner();
        alert(`Validation Error: ${err.message}`);
    });
}

function generateReport() {
    showSpinner("Compiling standalone Capstone HTML analysis report...");
    const snap = typeof BioPrism !== 'undefined' ? BioPrism.getSessionSnapshot(currentSessionId) : null;
    fetch('/api/generate-report', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ session_id: currentSessionId, session_snapshot: snap })
    })
    .then(r => r.json())
    .then(data => {
        hideSpinner();
        if (data.success && (data.preview_url || data.report_url)) {
            window.open(`/reports?session_id=${encodeURIComponent(currentSessionId)}`, '_blank');
        } else {
            alert("Could not generate report.");
        }
    })
    .catch(err => {
        hideSpinner();
        alert(`Report Generation Error: ${err.message}`);
    });
}

function resetPipeline() {
    if (confirm("Reset current pipeline session and start fresh?")) {
        window.location.href = '/pipeline';
    }
}
