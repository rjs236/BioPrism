import os
import sys

# Ensure the root application directory is in Python's search path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from app import app

class VercelPathFix:
    """WSGI middleware to strip /api/index destination path added by Vercel rewrites."""
    def __init__(self, wsgi_app):
        self.wsgi_app = wsgi_app

    def __call__(self, environ, start_response):
        path = environ.get('PATH_INFO', '')
        for prefix in ['/api/index.py', '/api/index']:
            if path.startswith(prefix):
                new_path = path[len(prefix):]
                if not new_path.startswith('/'):
                    new_path = '/' + new_path
                environ['PATH_INFO'] = new_path
                break
        return self.wsgi_app(environ, start_response)

app.wsgi_app = VercelPathFix(app.wsgi_app)

