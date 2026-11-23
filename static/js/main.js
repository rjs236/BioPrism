/**
 * BioPrism - Core Frontend Utilities & State Management
 */
const BioPrism = {
    // Get active session ID from URL or localStorage
    getSessionId: function() {
        const urlParams = new URLSearchParams(window.location.search);
        let sid = urlParams.get('session_id') || localStorage.getItem('bioprism_session_id');
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

    // Fetch JSON helper with defensive error parsing
    apiPost: async function(url, payload) {
        try {
            const res = await fetch(url, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });
            const text = await res.text();
            try {
                return JSON.parse(text);
            } catch (parseErr) {
                console.error("Non-JSON API response from " + url + ":", text);
                return { success: false, message: `Server error (${res.status}): Please check parameters.` };
            }
        } catch (err) {
            console.error("API POST error:", err);
            return { success: false, message: err.toString() };
        }
    }
};

document.addEventListener('DOMContentLoaded', () => {
    BioPrism.syncUrlSession();
    
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

