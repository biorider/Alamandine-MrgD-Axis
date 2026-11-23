"""Central paths for the Phase 0 pipeline. Override any path with an environment variable.
Run every script from the repository root:  python scripts/<name>.py"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
def P(*a): return os.path.join(ROOT, *a)

# User-supplied / locally held inputs (not redistributed; see data/manifests)
USER_DATA  = os.environ.get('PHASE0_USER_DATA',  r'C:\Users\u9urt\Desktop\Makale proje\Alamandine-study\data')          # DepMap files
NEW_AML    = os.environ.get('PHASE0_NEW_AML',    r'C:\Users\u9urt\Desktop\New-AML-Data')                                  # GSE116256/207104/241989 RAW tars
TCGA_STAR  = os.environ.get('PHASE0_TCGA_STAR',  r'C:\Users\u9urt\Desktop\Makale proje\ML-Performance\data\raw\TCGA_STAR')  # GDC STAR augmented counts (151 files)
TCGA_FILES = os.environ.get('PHASE0_TCGA_FILES', r'C:\Users\u9urt\Desktop\Makale proje\ML-Performance\data\raw\TCGA_LAML_files.json')

# Open data downloaded by 00_download_open_data.py
OPEN    = os.environ.get('PHASE0_OPEN', P('data', 'open'))
INTERIM = P('data', 'interim')
TABLES  = P('results', 'tables')
FIGS    = P('results', 'figures')
for d in (OPEN, INTERIM, TABLES, FIGS): os.makedirs(d, exist_ok=True)

GTF = os.path.join(OPEN, 'gencode.v36.annotation.gtf.gz')
RAS8 = ['AGT', 'REN', 'ACE', 'ACE2', 'AGTR1', 'AGTR2', 'MAS1', 'MRGPRD']
