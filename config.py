import os

# ---- Mycocosm configuration ----
MYCOCOSM_FUNGI_URL = "https://mycocosm.jgi.doe.gov/fungi/fungi.info.html"
PORTALS_TABLE_FILENAME = "mycocosm_fungi_data.csv"

# ---- Paths ----
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'local_data')