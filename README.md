# genephylogeny

**genephylogeny** is a Python toolkit for building putative phylogenies of fungal orthologs using publicly available proteomes and orthology analysis. It automates the process of fetching, wrangling, and cleaning fungal proteome data from the JGI MycoCosm portal, preparing it for downstream phylogenetic analysis.

---

## Features

- **Automated Data Fetching:** Downloads fungal genome and proteome metadata from MycoCosm.
- **File Listing & Metadata Extraction:** Retrieves and parses file listings for selected organisms via the JGI API.
- **Proteome Extraction & Renaming:** Extracts compressed proteome files and standardizes FASTA headers.
- **Sequence Filtering:** Removes sequences outside user-defined length thresholds.
- **Phylogeny Table Construction:** Builds tables summarizing taxonomic and file completeness information.

---

## Workflow Overview

1. **Fetch Fungal Portal Table**
   - `mycocosm_table_fetch.py`: Downloads the main fungal portal table from MycoCosm.

2. **Fetch File Listings for Selected Organisms**
   - `mycocosm_filelist_fetch.py`: Uses the JGI API to fetch file listings for each organism in your selection.

3. **Wrangle File Metadata**
   - `mycocosm_filelist_wrangle.py`: Processes file metadata and builds phylogeny summary tables.
   - Manual selection of interesting species should be stored in `proteome_list_orthofinder.csv`.

4. **Extract and Rename Proteomes**
   - `proteome_file_process.py` and `proteome_file_process_custom.py`: Extracts compressed proteome files and renames FASTA headers for consistency.

5. **Filter Sequences by Length**
   - `proteome_file_cleanup.py`: Removes sequences that are too short or too long.

---

## Setup

1. **Clone the Repository**
   ```sh
   git clone https://github.com/yourusername/genephylogeny.git
   cd genephylogeny
   ```

2. **Install Dependencies**
   ```sh
   pip install -r requirements.txt
   ```
   *Required packages include:* `pandas`, `biopython`, `requests`, `beautifulsoup4`

3. **Configure Credentials**
   - Edit `local_data/credentials.py` and add your JGI API token.

4. **Edit Configuration**
   - Adjust paths and settings in `config.py` as needed.

---

## Usage

Run each step in order, or adapt the workflow to your needs:

```sh
python mycocosm_table_fetch.py
python mycocosm_filelist_fetch.py
python mycocosm_filelist_wrangle.py
python proteome_file_cleanup.py
```

Intermediate and output files are saved in the `local_data` directory and its subfolders.

---

## File Structure

- `config.py` — Central configuration for paths and URLs.
- `local_data/` — All downloaded and processed data.
- `utils/` — Utility scripts for web requests and data wrangling.
- `mycocosm_table_fetch.py` — Downloads the fungal portal table.
- `mycocosm_filelist_fetch.py` — Fetches file listings for selected organisms.
- `mycocosm_filelist_wrangle.py` — Processes and summarizes file metadata.
- `proteome_file_cleanup.py` — Filters and cleans proteome FASTA files.

---

## Notes

- **Credentials:** Never share your JGI API token or upload your modified `credentials.py`.
- **Data Volume:** Downloading and processing all fungal proteomes may require significant disk space and time.
- **Customization:** You can adjust sequence length thresholds and other parameters in the relevant scripts.

---

## License

MIT License. See `LICENSE` for details.

---

## Citation

If you use genephylogeny in your research, please