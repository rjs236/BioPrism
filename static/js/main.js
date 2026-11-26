/**
 * BioPrism - Core Frontend Utilities & Cross-Page State Management
 */
const BioPrism = {
    // Get active session ID from URL or localStorage and keep both synchronized
    getSessionId: function() {
        const urlParams = new URLSearchParams(window.location.search);
        const urlSid = urlParams.get('session_id');
        if (urlSid) {
            localStorage.setItem('bioprism_session_id', urlSid);
            return urlSid;
        }
        let sid = localStorage.getItem('bioprism_session_id');
        if (!sid) {
            sid = Math.random().toString(36).substring(2, 10);
            localStorage.setItem('bioprism_session_id', sid);
        }
        return sid;
    },

    setSessionId: function(sid) {
        if (sid) {
            localStorage.setItem('bioprism_session_id', sid);
        }
    },

    // Ensure session_id is appended to URL if navigating
    syncUrlSession: function() {
        const urlParams = new URLSearchParams(window.location.search);
        const sid = this.getSessionId();
        if (!urlParams.get('session_id')) {
            const currentUrl = new URL(window.location.href);
            currentUrl.searchParams.set('session_id', sid);
            window.history.replaceState({}, '', currentUrl);
        }
        return sid;
    },

    // Propagate active session_id to internal navigation and session download links
    updateInternalLinks: function(sid) {
        if (!sid) return;
        const staticFiles = ['template_expression', 'template_metadata', 'demo_expression', 'demo_metadata'];
        document.querySelectorAll('a[href^="/"]').forEach(a => {
            const rawHref = a.getAttribute('href');
            if (!rawHref || rawHref.startsWith('//') || rawHref.startsWith('/static/')) return;
            // Skip template/demo downloads
            if (staticFiles.some(tf => rawHref.includes('/api/download/' + tf))) return;
            try {
                const u = new URL(rawHref, window.location.origin);
                // Preserve explicit session_id if link already targets a specific history run (unless it was an empty server placeholder)
                const existingSid = u.searchParams.get('session_id');
                if (!existingSid || existingSid === window.BIOPRISM_SERVER_SESSION_ID) {
                    u.searchParams.set('session_id', sid);
                    a.setAttribute('href', u.pathname + u.search + u.hash);
                }
            } catch (e) {}
        });
    },

    // Persist lightweight session snapshot in localStorage (resilient to serverless cold starts)
    saveSessionSnapshot: function(sid, snapshot) {
        if (!sid || !snapshot) return;
        const key = 'bioprism_session_' + sid;
        try {
            let toStore = Object.assign({}, snapshot);
            let serialized = JSON.stringify(toStore);
            // If > 2.2 MB, compact bulky plot traces/raw CSV before caching
            if (serialized.length > 2200000) {
                toStore = Object.assign({}, snapshot, {
                    volcano_plot: null,
                    pathway_bar_plot: null,
                    network_json: null
                });
                if (toStore.de_summary) {
                    toStore.de_summary = Object.assign({}, toStore.de_summary, { heatmap_plot: null, pca_plot: null });
                }
                if (toStore.ml_metrics) {
                    toStore.ml_metrics = Object.assign({}, toStore.ml_metrics, { roc_plot: null, confusion_matrix_plot: null, expression_comparison_plot: null });
                }
                serialized = JSON.stringify(toStore);
            }
            localStorage.setItem(key, serialized);
        } catch (err) {
            // Quota exceeded fallback: prune oldest cached sessions
            try {
                Object.keys(localStorage).forEach(k => {
                    if (k.startsWith('bioprism_session_') && k !== key) {
                        localStorage.removeItem(k);
                    }
                });
                localStorage.setItem(key, JSON.stringify({
                    session_id: snapshot.session_id,
                    dataset_name: snapshot.dataset_name,
                    preproc_info: snapshot.preproc_info,
                    de_summary: snapshot.de_summary,
                    pathway_results: snapshot.pathway_results,
                    network_summary: snapshot.network_summary,
                    hub_genes: snapshot.hub_genes,
                    ml_metrics: snapshot.ml_metrics
                }));
            } catch (e2) {}
        }
    },

    getSessionSnapshot: function(sid) {
        if (!sid) return null;
        try {
            const raw = localStorage.getItem('bioprism_session_' + sid);
            return raw ? JSON.parse(raw) : null;
        } catch (e) {
            return null;
        }
    },

    // Client-side history synchronization (merges cleanly with server SQLite provenance)
    getDeletedHistoryIds: function() {
        try {
            return JSON.parse(localStorage.getItem('bioprism_deleted_history') || '[]');
        } catch (e) {
            return [];
        }
    },

    saveHistoryRecord: function(rec) {
        if (!rec || !rec.id) return;
        try {
            // Unmark if re-executed
            let deleted = this.getDeletedHistoryIds().filter(id => id !== rec.id);
            localStorage.setItem('bioprism_deleted_history', JSON.stringify(deleted));

            let list = this.getHistoryRecords(true);
            const idx = list.findIndex(item => item.id === rec.id);
            if (idx >= 0) {
                list[idx] = Object.assign({}, list[idx], rec);
            } else {
                list.unshift(rec);
            }
            // Keep up to 40 most recent runs locally
            localStorage.setItem('bioprism_history_runs', JSON.stringify(list.slice(0, 40)));
        } catch (e) {}
    },

    getHistoryRecords: function(includeDeleted = false) {
        try {
            let list = JSON.parse(localStorage.getItem('bioprism_history_runs') || '[]');
            if (!includeDeleted) {
                const deleted = new Set(this.getDeletedHistoryIds());
                list = list.filter(r => !deleted.has(r.id));
            }
            return list;
        } catch (e) {
            return [];
        }
    },

    deleteHistoryRecord: function(runId) {
        if (!runId) return;
        try {
            let list = this.getHistoryRecords(true).filter(r => r.id !== runId);
            localStorage.setItem('bioprism_history_runs', JSON.stringify(list));
            let deleted = this.getDeletedHistoryIds();
            if (!deleted.includes(runId)) deleted.push(runId);
            localStorage.setItem('bioprism_deleted_history', JSON.stringify(deleted));
        } catch (e) {}
    },

    clearHistoryRecords: function() {
        try {
            localStorage.removeItem('bioprism_history_runs');
            localStorage.removeItem('bioprism_deleted_history');
        } catch (e) {}
    },

    // Auto-restore server session when navigating to workspace pages across serverless containers
    restoreServerSessionIfNeeded: async function(sid) {
        const snap = this.getSessionSnapshot(sid);
        if (!snap || !snap.preproc_info) return;

        // Update topbar dataset label even on non-reload pages
        const label = document.getElementById('activeDatasetLabel');
        if (label && (!label.innerText || label.innerText.trim() === 'No Dataset Loaded' || label.innerText.trim() === 'Custom Upload')) {
            if (snap.dataset_name) label.innerText = snap.dataset_name;
        }

        const workspacePaths = [
            '/preprocessing', '/differential-expression', '/pathway-enrichment',
            '/network-analysis', '/biomarker-validation', '/ai-assistant', '/pipeline'
        ];
        const path = window.location.pathname;
        if (!workspacePaths.includes(path)) return;

        // If server already rendered this session with data, no reload needed
        if (window.BIOPRISM_SERVER_HAS_SESSION && window.BIOPRISM_SERVER_SESSION_ID === sid) return;

        const guardKey = 'bioprism_restored_' + path + '_' + sid;
        if (sessionStorage.getItem(guardKey)) return;
        sessionStorage.setItem(guardKey, '1');

        try {
            const resp = await fetch('/api/session/restore', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ session_id: sid, session_snapshot: snap })
            });
            const r = await resp.json();
            if (r.success) {
                window.location.href = path + '?session_id=' + encodeURIComponent(sid);
            }
        } catch (e) {}
    },

    // Download helper that restores session snapshot via POST so exports always work on serverless
    bindSnapshotDownloads: function(sid) {
        document.querySelectorAll('a[href*="/api/download/"]').forEach(link => {
            const href = link.getAttribute('href') || '';
            if (
                href.includes('template_') ||
                href.includes('demo_') ||
                link.id === 'btnDownloadHtmlReport'
            ) {
                return;
            }
            link.addEventListener('click', async function(e) {
                if (this.hasAttribute('disabled') || this.classList.contains('disabled')) {
                    e.preventDefault();
                    return;
                }
                const activeSid = BioPrism.getSessionId();
                const snap = BioPrism.getSessionSnapshot(activeSid);
                if (!snap) return; // Fall back to normal GET
                e.preventDefault();
                try {
                    const targetUrl = new URL(this.getAttribute('href'), window.location.origin);
                    targetUrl.searchParams.set('session_id', activeSid);
                    const resp = await fetch(targetUrl.pathname + targetUrl.search, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ session_id: activeSid, session_snapshot: snap })
                    });
                    if (!resp.ok) {
                        BioPrism.notify("Please execute the corresponding analysis step before exporting.", "error");
                        return;
                    }
                    const blob = await resp.blob();
                    let filename = 'bioprism_export.csv';
                    const disp = resp.headers.get('Content-Disposition');
                    if (disp && disp.includes('filename=')) {
                        filename = disp.split('filename=')[1].replace(/["']/g, '').trim();
                    }
                    const dlUrl = URL.createObjectURL(blob);
                    const tempA = document.createElement('a');
                    tempA.href = dlUrl;
                    tempA.download = filename;
                    document.body.appendChild(tempA);
                    tempA.click();
                    tempA.remove();
                    setTimeout(() => URL.revokeObjectURL(dlUrl), 1000);
                } catch (err) {
                    BioPrism.notify("Download failed: " + err.message, "error");
                }
            });
        });
    },

    // Show modern toast / alert notification
    notify: function(message, type = 'info') {
        const alertBox = document.createElement('div');
        alertBox.className = `alert alert-${type === 'error' ? 'danger' : type} alert-dismissible fade show position-fixed top-0 end-0 m-4 shadow-lg`;
        alertBox.style.zIndex = '99999';
        alertBox.style.minWidth = '300px';
        alertBox.innerHTML = `
            <div class="d-flex align-items-center gap-2">
                <i class="bi bi-${type === 'success' ? 'check-circle-fill' : type === 'error' ? 'exclamation-octagon-fill' : 'info-circle-fill'}"></i>
                <div>${message}</div>
            </div>
            <button type="button" class="btn-close" data-bs-dismiss="alert" aria-label="Close"></button>
        `;
        document.body.appendChild(alertBox);
        setTimeout(() => {
            alertBox.classList.remove('show');
            setTimeout(() => alertBox.remove(), 250);
        }, 4500);
    },

    // Render Plotly plot from JSON specification
    renderPlotly: function(elementId, plotSpec) {
        const el = document.getElementById(elementId);
        if (!el || !plotSpec) return;
        const data = plotSpec.data || [];
        const layout = plotSpec.layout || {};
        layout.autosize = true;
        layout.responsive = true;
        Plotly.newPlot(elementId, data, layout, { responsive: true, displayModeBar: true });
    },

    // Fetch JSON helper with automatic snapshot persistence & cold-start session restoration
    apiPost: async function(url, payload) {
        const activeSid = (payload && payload.session_id) || this.getSessionId();
        const sendReq = async (bodyObj) => {
            const res = await fetch(url, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(bodyObj)
            });
            const text = await res.text();
            try {
                return { status: res.status, data: JSON.parse(text) };
            } catch (parseErr) {
                console.error("Non-JSON API response from " + url + ":", text);
                return { status: res.status, data: { success: false, message: `Server error (${res.status}): Please check parameters.` } };
            }
        };

        try {
            let { data } = await sendReq(payload || {});
            // If server container is stateless/cold, attach client snapshot and retry once
            if (data && data.needs_restore && activeSid) {
                const snap = this.getSessionSnapshot(activeSid);
                if (snap) {
                    const retryPayload = Object.assign({}, payload, { session_snapshot: snap });
                    ({ data } = await sendReq(retryPayload));
                }
            }

            if (data && data.session_snapshot && activeSid) {
                this.saveSessionSnapshot(data.session_id || activeSid, data.session_snapshot);
            }
            if (data && data.history_record) {
                this.saveHistoryRecord(data.history_record);
            }
            return data;
        } catch (err) {
            console.error("API POST error:", err);
            return { success: false, message: err.toString() };
        }
    }
};

document.addEventListener('DOMContentLoaded', () => {
    const sid = BioPrism.syncUrlSession();
    BioPrism.updateInternalLinks(sid);
    BioPrism.bindSnapshotDownloads(sid);
    BioPrism.restoreServerSessionIfNeeded(sid);

    // Animate progress bar as DOM becomes interactive
    const pBar = document.getElementById('preloaderProgressBar');
    if (pBar) {
        pBar.style.width = '65%';
    }
});

// Handle Opening Splash Loader Dismissal on Window Load
window.addEventListener('load', () => {
    const preloader = document.getElementById('bioprismPreloader');
    const pBar = document.getElementById('preloaderProgressBar');
    const pStatus = document.getElementById('preloaderStatusText');

    if (pBar) pBar.style.width = '100%';
    if (pStatus) {
        pStatus.innerHTML = '<span class="text-success"><i class="bi bi-check2-circle me-1"></i>Analysis Environment Ready</span>';
    }

    setTimeout(() => {
        if (preloader) {
            preloader.classList.add('fade-out');
            setTimeout(() => {
                preloader.style.display = 'none';
            }, 460);
        }
    }, 320);
});

// Safety fallback: dismiss preloader after 1.6s in case of delayed network assets
setTimeout(() => {
    const preloader = document.getElementById('bioprismPreloader');
    if (preloader && !preloader.classList.contains('fade-out')) {
        preloader.classList.add('fade-out');
        setTimeout(() => {
            preloader.style.display = 'none';
        }, 460);
    }
}, 1600);
