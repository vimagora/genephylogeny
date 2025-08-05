import os
import json
import requests
import time
import csv
import pandas as pd
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from config import JGI_API_BASE_URL, FILES_PER_PAGE, JSON_FOLDER, REQUEST_DELAY, ALL_FILES_METADATA_PATH

def download_mycocosm_fungi_table(url, output_csv_file):
    """
    Downloads the HTML table from the given MycoCosm URL and saves it as a CSV file.

    Args:
        url (str): The URL of the web page containing the table.
        output_csv_file (str): The name of the CSV file to save the data to.
    """
    try:
        # Fetch the HTML content
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, 'html.parser')

        # Find the first table (can adjust if needed)
        table = soup.find('table')
        headers = [th.get_text(strip=True) for th in table.find_all('th')]

        rows = []
        for tr in table.find_all('tr')[1:]:  # skip header
            cells = tr.find_all('td')
            row_data = []
            row_links = []

            for cell in cells:
                # Get text
                cell_text = cell.get_text(strip=True)
                row_data.append(cell_text)

                # Get link if it exists
                link_tag = cell.find('a')
                if link_tag and link_tag.get('href'):
                    full_url = urljoin(url, link_tag['href'])
                    row_links.append(full_url)
                else:
                    row_links.append(None)

            # Combine text and links into a flat list
            combined_row = row_data + row_links
            rows.append(combined_row)

        # Add headers for links
        link_headers = [f"{col}_link" for col in headers]
        final_headers = headers + link_headers

        # Create DataFrame, clean up the unnecessary columns and characters, and save
        df = pd.DataFrame(rows, columns=final_headers)
        
        df['portal'] = df['Name_link'].str.replace('https://mycocosm.jgi.doe.gov/', '')
        df['reference'] = df['Published_link']
        df = df.loc[:, ~df.columns.str.endswith('_link')]

        df.to_csv(output_csv_file, index=False, encoding='utf-8')

        print(f"\n✅ Successfully saved table with links to {output_csv_file}")

    except Exception as e:
        print(f"❌ Error: {e}")


def fetch_all_files(organism_id, header):
    """
    Fetch the json listing of all files for a given organism ID from JGI API.

    Args:
        organism_id (str): The ID of the organism to fetch files for.
    """
    print(f"Fetching all files for {organism_id} from JGI...")
    params = {
        "organism": organism_id,
        "api_version": 2,
        "a": "false",  # Exclude archived
        "h": "false",  # Exclude hidden
        "d": "asc",    # Sort ascending
        "p": 1,        # Page number (start at 1)
        "x": FILES_PER_PAGE,  # Files per page
        "t": "simple"
    }

    # Create the JSON folder if it doesn't exist
    os.makedirs(JSON_FOLDER, exist_ok=True)

    page = 1
    total_files = 0

    try:
        while True:
            params["p"] = page  # Update page number
            print(f"Fetching page {page} for {organism_id}...")
            response = requests.get(JGI_API_BASE_URL, params=params, headers=header)
            response.raise_for_status()
            data = response.json()

            # Get files from the current page
            current_files = data.get("organisms", [{}])[0].get("files", [])
            if not current_files:  # Break if no files are returned
                print(f"No more files found for {organism_id}. Stopping pagination.")
                break

            total_files += len(current_files)
            print(f"Found {len(current_files)} files on page {page} for {organism_id}.")

            # Save current page to a JSON file in the JSON folder
            page_filename = os.path.join(JSON_FOLDER, f"all_files_{organism_id}_page_{page}.json")
            with open(page_filename, "w") as f:
                json.dump(data, f, indent=2)
            print(f"✅ Saved page {page} to {page_filename}")

            page += 1
            time.sleep(REQUEST_DELAY)  # Delay between requests

        if total_files == 0:
            print(f"❌ No files found for {organism_id} across any pages.")
        else:
            print(f"✅ Total {total_files} files saved for {organism_id}.")
        return total_files > 0

    except requests.exceptions.HTTPError as e:
        print(f"HTTP error for {organism_id}: {e} - {response.status_code}: {response.text}")
        return False
    except requests.exceptions.RequestException as e:
        print(f"Request error for {organism_id}: {e}")
        return False
    except Exception as e:
        print(f"General error for {organism_id}: {e}")
        return False
    
def parse_and_export(organism_ids):
    found = []
    # Ensure the JSON folder exists
    if not os.path.exists(JSON_FOLDER):
        print(f"No {JSON_FOLDER} folder found. Please run fetch_all_files() first.")
        return

    for organism_id in organism_ids:
        json_files = [f for f in os.listdir(JSON_FOLDER) if f.startswith(f"all_files_{organism_id}_page_") and f.endswith(".json")]
        if not json_files:
            print(f"No JSON files found for {organism_id} in {JSON_FOLDER}. Skipping.")
            found.append({
                "organism": organism_id,
                "file_name": "NO FILES FOUND",
                "file_id": "",
                "_id": "",
                "file_status": "",
                "md5sum": "",
                "file_date": "",
                "ncbi_taxon_id": "",
                "jat_label": "",
                "ncbi_taxon_class": "",
                "ncbi_taxon_family": "",
                "ncbi_taxon_order": "",
                "ncbi_taxon_genus": "",
                "ncbi_taxon_species": "",
                "file_type": "",
                "portal_display_location": ""
            })
            continue

        print(f"Processing files for {organism_id}...")
        organism_file_count = 0
        
        for json_file in sorted(json_files):  # Sort to process pages in order
            json_file_path = os.path.join(JSON_FOLDER, json_file)
            print(f"Parsing {json_file}...")
            with open(json_file_path, "r") as f:
                data = json.load(f)

            files = data.get("organisms", [{}])[0].get("files", [])
            organism_file_count += len(files)
            print(f"Found {len(files)} files in {json_file}.")

            for file in files:
                metadata = file.get("metadata", {})
                ncbi_taxon = metadata.get("ncbi_taxon", {})
                portal = metadata.get("portal", {})
                
                found.append({
                    "organism": organism_id,
                    "file_name": file.get("file_name"),
                    "file_id": file.get("file_id"),
                    "_id": file.get("_id"),
                    "file_status": file.get("file_status"),
                    "md5sum": file.get("md5sum"),
                    "file_date": file.get("file_date"),
                    "ncbi_taxon_id": metadata.get("ncbi_taxon_id", ""),
                    "jat_label": metadata.get("jat_label", ""),
                    "ncbi_taxon_class": ncbi_taxon.get("ncbi_taxon_class", ""),
                    "ncbi_taxon_family": ncbi_taxon.get("ncbi_taxon_family", ""),
                    "ncbi_taxon_order": ncbi_taxon.get("ncbi_taxon_order", ""),
                    "ncbi_taxon_genus": ncbi_taxon.get("ncbi_taxon_genus", ""),
                    "ncbi_taxon_species": ncbi_taxon.get("ncbi_taxon_species", ""),
                    "file_type": file.get("file_type", ""),
                    "portal_display_location": portal.get("display_location", "")
                })
                    
        # If no match found across all pages
        if organism_file_count == 0:
            print(f"⚠️ No files in parsed pages for {organism_id}.")
            found.append({
                "organism": organism_id,
                "file_name": "NO FILES FOUND",
                "file_id": "",
                "_id": "",
                "file_status": "",
                "md5sum": "",
                "file_date": "",
                "ncbi_taxon_id": "",
                "jat_label": "",
                "ncbi_taxon_class": "",
                "ncbi_taxon_family": "",
                "ncbi_taxon_order": "",
                "ncbi_taxon_genus": "",
                "ncbi_taxon_species": "",
                "file_type": "",
                "portal_display_location": ""
            })

        # Export to CSV
        output_file = ALL_FILES_METADATA_PATH
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "organism", "file_name", "file_id", "_id", "file_status", "md5sum",
                "file_date", "ncbi_taxon_id", "jat_label", "ncbi_taxon_class",
                "ncbi_taxon_family", "ncbi_taxon_order", "ncbi_taxon_genus",
                "ncbi_taxon_species", "file_type", "portal_display_location"
            ])
            writer.writeheader()
            writer.writerows(found)
    