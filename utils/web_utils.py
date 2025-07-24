import requests
from bs4 import BeautifulSoup
import pandas as pd
from urllib.parse import urljoin

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