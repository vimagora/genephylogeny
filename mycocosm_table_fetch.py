import os
from config import MYCOCOSM_FUNGI_URL, PORTALS_TABLE_PATH
from utils.web_utils import download_mycocosm_fungi_table

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    download_mycocosm_fungi_table(MYCOCOSM_FUNGI_URL, PORTALS_TABLE_PATH)