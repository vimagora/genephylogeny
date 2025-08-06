import os
import csv
from config import JSON_DIR, ORGANISM_IDS_PATH
from utils.web_utils import fetch_all_files, parse_and_export
from local_data.credentials import JGI_API_TOKEN

# Authentication headers
headers = {
    "accept": "application/json",
    "Authorization": JGI_API_TOKEN
}    

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    organism_ids = []
    with open(ORGANISM_IDS_PATH, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if "portal" in row:
                organism_ids.append(row["portal"])
    print(organism_ids)

    os.makedirs(JSON_DIR, exist_ok=True)

    for organism_id in organism_ids:
        # Check if any JSON files exist for this organism to avoid re-fetching
        if not any(f.startswith(f"all_files_{organism_id}_page_") for f in os.listdir(JSON_DIR)):
            fetch_all_files(organism_id, headers)
        else:
            print(f"🗂️ Using cached JSON for {organism_id}...")

    parse_and_export(organism_ids)