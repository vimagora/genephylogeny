import os

# ---- Mycocosm configuration ----
## Genome table
MYCOCOSM_FUNGI_URL = "https://mycocosm.jgi.doe.gov/fungi/fungi.info.html"


## API configuration
BASE_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
FILES_PER_PAGE = 50  # Set to 50 files per page
REQUEST_DELAY = 1  # Delay in seconds between requests

# ---- Paths ----
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'local_data')
JSON_DIR = os.path.join(DATA_DIR, 'json_files')

PORTALS_TABLE_PATH = os.path.join(DATA_DIR, "mycocosm_fungi_data.csv")
ORGANISM_IDS_PATH = os.path.join(DATA_DIR, 'selectedorganism_ids.csv')
MYCOCOSM_FILELIST_PATH = os.path.join(DATA_DIR, 'mycocosm_data.csv')
