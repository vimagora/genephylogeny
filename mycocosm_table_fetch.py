import os
from config import MYCOCOSM_FUNGI_URL, PORTALS_TABLE_FILENAME, DATA_DIR
from utils.web_utils import download_mycocosm_fungi_table

# Run the function
portals_csv_path = os.path.join(DATA_DIR, PORTALS_TABLE_FILENAME)

download_mycocosm_fungi_table(MYCOCOSM_FUNGI_URL, portals_csv_path)