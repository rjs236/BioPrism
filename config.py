import os
import shutil
import tempfile

BASE_DIR = os.path.abspath(os.path.dirname(__file__))

class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY', 'bioprism-capstone-secret-key-2025')
    SAMPLE_DATA_FOLDER = os.path.join(BASE_DIR, 'static', 'sample_data')
    R_SCRIPTS_FOLDER = os.path.join(BASE_DIR, 'r_scripts')
    
    # Auto-detect Rscript path
    _custom_r_path = r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
    if os.path.exists(_custom_r_path):
        RSCRIPT_PATH = _custom_r_path
    else:
        RSCRIPT_PATH = shutil.which("Rscript") or "Rscript"
        
    MAX_CONTENT_LENGTH = 64 * 1024 * 1024  # 64 MB max upload
    ALLOWED_EXTENSIONS = {'csv', 'tsv', 'txt'}
    
    # SQLite Database & dynamic storage paths (supports standard local & Vercel serverless tempdir)
    if os.environ.get('VERCEL'):
        _tmp = tempfile.gettempdir()
        DATABASE_PATH = os.path.join(_tmp, 'bioprism.db')
        UPLOAD_FOLDER = os.path.join(_tmp, 'uploads')
        RESULTS_FOLDER = os.path.join(_tmp, 'results')
    else:
        DATABASE_PATH = os.path.join(BASE_DIR, 'bioprism.db')
        UPLOAD_FOLDER = os.path.join(BASE_DIR, 'uploads')
        RESULTS_FOLDER = os.path.join(BASE_DIR, 'results')
