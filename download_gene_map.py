"""
Download gene mapping data from STRING database.

File: download_gene_map.py
Author: Akaash Venkat, Audi Liu
Updated: 2026 (Python 3 modernization)

This module handles downloading and processing gene mapping data from the STRING
protein interaction database, organizing genes into groups based on connectivity.
"""

import logging
import os
import re
import time
from pathlib import Path
from typing import Dict, List, Tuple

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuration
INFO_DIR = Path("info_files")
INPUT_DIR = Path("input_files")
SVG_DIR = Path("svg_files")

GENE_DATABASE_FILE = INFO_DIR / "gene_database.txt"
GENE_GROUP_FILE = INFO_DIR / "gene_group.txt"
INTERMEDIATE_GENES_FILE = INFO_DIR / "intermediate_genes.txt"
UNIDENTIFIABLE_GENE_FILE = INFO_DIR / "unidentifiable_genes.txt"
CHANGED_NAME_GENE_FILE = INFO_DIR / "changed_name_genes.txt"

# Global data structures
GENE_LIST: List[Tuple[str, Dict[str, float]]] = []
UNIDENTIFIABLE_LIST: List[str] = []
CHANGED_NAME: Dict[str, str] = {}
GROUP: Dict[str, str] = {}
B_D_PAIR: Dict[str, str] = {}

# Constants
STRING_BASE_URL = "https://string-db.org/"
SPECIES = "Homo sapiens"
CONNECTION_LIMIT = "500"
WAIT_TIMEOUT = 20


def setup_driver() -> webdriver.Chrome:
    """Initialize and return a Chrome WebDriver with optimized options."""
    options = Options()
    options.add_argument('--no-sandbox')
    options.add_argument('--disable-dev-shm-usage')
    return webdriver.Chrome(options=options)


def read_database() -> None:
    """Parse gene database file into GENE_LIST."""
    if not GENE_DATABASE_FILE.exists():
        logger.warning(f"Database file not found: {GENE_DATABASE_FILE}")
        return

    try:
        with open(GENE_DATABASE_FILE, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Parse: "GENE - NEIGHBOR1(confidence), NEIGHBOR2(confidence), ..."
                parts = line.split(' - ', 1)
                if len(parts) != 2:
                    continue

                gene_name = parts[0].strip()
                neighbors_str = parts[1]

                neighbors = {}
                for neighbor_entry in neighbors_str.split(', '):
                    match = re.match(r'(\w+)\(([\d.]+)\)', neighbor_entry.strip())
                    if match:
                        neighbor_name, confidence = match.groups()
                        neighbors[neighbor_name] = float(confidence)

                GENE_LIST.append((gene_name, neighbors))
        logger.info(f"Loaded {len(GENE_LIST)} genes from database")
    except Exception as e:
        logger.error(f"Error reading database: {e}")


def read_unidentifiable() -> None:
    """Read list of genes that couldn't be identified."""
    if not UNIDENTIFIABLE_GENE_FILE.exists():
        return

    try:
        with open(UNIDENTIFIABLE_GENE_FILE, 'r') as f:
            for line in f:
                line = line.strip()
                if line and "cannot be found" not in line:
                    UNIDENTIFIABLE_LIST.append(line)
        logger.info(f"Loaded {len(UNIDENTIFIABLE_LIST)} unidentifiable genes")
    except Exception as e:
        logger.error(f"Error reading unidentifiable genes: {e}")


def read_changed_names() -> None:
    """Read list of genes that were renamed."""
    if not CHANGED_NAME_GENE_FILE.exists():
        return

    try:
        with open(CHANGED_NAME_GENE_FILE, 'r') as f:
            for line in f:
                line = line.strip()
                if line and "renamed" not in line:
                    if '=>' in line:
                        orig, new = line.split('=>')
                        CHANGED_NAME[orig.strip()] = new.strip()
        logger.info(f"Loaded {len(CHANGED_NAME)} renamed genes")
    except Exception as e:
        logger.error(f"Error reading changed names: {e}")


def write_database() -> None:
    """Write GENE_LIST back to database file."""
    try:
        with open(GENE_DATABASE_FILE, 'w') as f:
            for gene_name, neighbors in sorted(GENE_LIST):
                neighbors_str = ', '.join(
                    f"{name}({value})"
                    for name, value in sorted(neighbors.items())
                )
                f.write(f"{gene_name} - {neighbors_str}\n\n")
        logger.info("Database file updated")
    except Exception as e:
        logger.error(f"Error writing database: {e}")


def write_gene_groups() -> None:
    """Write gene group classifications to file."""
    try:
        with open(GENE_GROUP_FILE, 'w') as f:
            groups = {
                'A': "Input gene that has direct connection with another input gene",
                'B': "Input gene that is indirectly connected to another input gene, via an intermediate gene",
                'C': "Input gene that is not directly or indirectly connected to another input gene",
                'D': "Intermediate gene that connects Group B genes with Group A or other Group B genes"
            }

            for group_id, description in groups.items():
                f.write(f"Group {group_id}: {description}\n---\n")
                cluster = get_list_for_group(group_id)
                for gene in cluster:
                    f.write(f"{gene}\n")
                f.write("\n\n\n")
        logger.info("Gene groups file updated")
    except Exception as e:
        logger.error(f"Error writing gene groups: {e}")


def write_intermediate_genes() -> None:
    """Write B-D gene pairings."""
    try:
        if B_D_PAIR:
            with open(INTERMEDIATE_GENES_FILE, 'w') as f:
                f.write("The following pairings (B Genes : D Genes) indicate that the B Gene requires "
                       "its respective D Gene to serve as an intermediate gene.\n\n")
                for b_gene, d_gene in sorted(B_D_PAIR.items()):
                    f.write(f"{b_gene} : {d_gene}\n")
            logger.info(f"Wrote {len(B_D_PAIR)} intermediate gene pairs")
        else:
            INTERMEDIATE_GENES_FILE.unlink(missing_ok=True)
    except Exception as e:
        logger.error(f"Error writing intermediate genes: {e}")


def write_unidentifiable() -> None:
    """Write unidentifiable genes."""
    try:
        unique_unidentifiable = list(dict.fromkeys(UNIDENTIFIABLE_LIST))
        if unique_unidentifiable:
            with open(UNIDENTIFIABLE_GENE_FILE, 'w') as f:
                f.write("The following genes cannot be found on the online STRING database:\n\n")
                for gene in unique_unidentifiable:
                    f.write(f"{gene}\n")
            logger.info(f"Wrote {len(unique_unidentifiable)} unidentifiable genes")
        else:
            UNIDENTIFIABLE_GENE_FILE.unlink(missing_ok=True)
    except Exception as e:
        logger.error(f"Error writing unidentifiable genes: {e}")


def write_changed_names() -> None:
    """Write renamed genes."""
    try:
        if CHANGED_NAME:
            with open(CHANGED_NAME_GENE_FILE, 'w') as f:
                f.write("The following genes have been renamed per the STRING database:\n\n")
                for orig, new in CHANGED_NAME.items():
                    f.write(f"{orig} => {new}\n")
            logger.info(f"Wrote {len(CHANGED_NAME)} renamed genes")
        else:
            CHANGED_NAME_GENE_FILE.unlink(missing_ok=True)
    except Exception as e:
        logger.error(f"Error writing changed names: {e}")


def initialize_connections() -> None:
    """Initialize all genes with group 'C'."""
    for gene_name, _ in GENE_LIST:
        GROUP[gene_name] = "C"


def identify_group_a() -> None:
    """Identify genes with direct connections (Group A)."""
    gene_names = {gene[0] for gene in GENE_LIST}

    for gene_name, neighbors in GENE_LIST:
        if any(neighbor in gene_names for neighbor in neighbors.keys()):
            GROUP[gene_name] = "A"

    logger.info(f"Identified {sum(1 for g in GROUP.values() if g == 'A')} Group A genes")


def identify_group_b() -> None:
    """Identify genes with indirect connections (Group B)."""
    for i, (gene_name, neighbors) in enumerate(GENE_LIST):
        if GROUP[gene_name] == "A":
            continue

        best_match = None
        best_score = -1

        for j, (other_gene, other_neighbors) in enumerate(GENE_LIST):
            if i == j:
                continue

            for inter_gene, confidence in neighbors.items():
                if inter_gene in other_neighbors:
                    score = min(confidence, other_neighbors[inter_gene])
                    if score > best_score:
                        best_score = score
                        best_match = (score, inter_gene, other_gene)

        if best_match and GROUP[gene_name] == "C":
            GROUP[gene_name] = "B"
            GROUP[best_match[1]] = "D"
            B_D_PAIR[gene_name] = best_match[1]

    logger.info(f"Identified {sum(1 for g in GROUP.values() if g == 'B')} Group B genes")
    logger.info(f"Identified {sum(1 for g in GROUP.values() if g == 'D')} Group D genes")


def get_list_for_group(group_id: str) -> List[str]:
    """Get all genes in a specific group."""
    return sorted([gene for gene, group in GROUP.items() if group == group_id])


def find_neighbor(input_gene: str) -> Dict[str, float] | int | str:
    """
    Query STRING database for gene neighbors.
    
    Returns:
        - dict: neighbor genes with confidence scores
        - -1: gene not found
        - str: corrected gene name
    """
    driver = setup_driver()
    try:
        wait = WebDriverWait(driver, WAIT_TIMEOUT)

        # Navigate to STRING
        driver.get(STRING_BASE_URL)
        wait.until(EC.element_to_be_clickable((By.ID, "search"))).click()

        # Fill in search form
        wait.until(EC.presence_of_element_located((By.ID, "primary_input:single_identifier"))).send_keys(input_gene)
        wait.until(EC.presence_of_element_located((By.ID, "species_text_single_identifier"))).send_keys(SPECIES)

        # Submit
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='input_form_single_identifier']/div[4]/a")
        )).click()

        time.sleep(5)
        page_source = driver.page_source

        # Check if protein found
        if "Sorry, STRING did not find a protein" in page_source:
            logger.warning(f"Gene not found: {input_gene}")
            return -1

        # Handle disambiguation
        if "Please select one" in page_source:
            wait.until(EC.element_to_be_clickable(
                (By.XPATH, "//*[@id='proceed_form']/div[1]/div/div[2]/a[2]")
            )).click()

        time.sleep(15)

        # Get correct gene name
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='bottom_page_selector_settings']")
        )).click()
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='bottom_page_selector_legend']")
        )).click()

        time.sleep(5)
        page_source = driver.page_source
        match = re.search(r'<td class="td_name middle_row first_row last_row"[^>]*>">([^<]+)</td>', page_source)
        if match:
            correct_name = match.group(1)
            if input_gene != correct_name:
                logger.info(f"Gene renamed: {input_gene} -> {correct_name}")
                return correct_name

        # Configure connection limit
        time.sleep(15)
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='bottom_page_selector_table']")
        )).click()
        wait.until(EC.element_to_be_clickable(
            (By.ID, "bottom_page_selector_settings")
        )).click()
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='standard_parameters']/div/div[1]/div[3]/div[2]/div[2]/div[1]/label")
        )).click()
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//select[@name='limit']/option[text()='custom value']")
        )).click()

        custom_input = wait.until(EC.presence_of_element_located((By.ID, "custom_limit_input")))
        custom_input.clear()
        custom_input.send_keys(CONNECTION_LIMIT)

        time.sleep(5)
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='standard_parameters']/div/div[1]/div[5]/a")
        )).click()

        time.sleep(20)

        # Fetch connection data
        wait.until(EC.element_to_be_clickable(
            (By.ID, "bottom_page_selector_table")
        )).click()
        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='bottom_page_selector_legend']")
        )).click()

        gene_connectors = {}
        connectors = wait.until(EC.presence_of_all_elements_located(
            (By.CLASS_NAME, "linked_item_row")
        ))

        for connector in connectors:
            text = connector.text.split('\n')
            if len(text) >= 2:
                neighbor = text[0].strip()
                confidence = float(text[-1].strip())
                gene_connectors[neighbor] = confidence

        logger.info(f"Found {len(gene_connectors)} neighbors for {input_gene}")
        return gene_connectors

    except Exception as e:
        logger.error(f"Error querying gene {input_gene}: {e}")
        return {}
    finally:
        driver.quit()


def parse_input() -> None:
    """Parse input gene list and fetch neighbor data."""
    input_file = INPUT_DIR / "original_gene_list.txt"
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return

    try:
        with open(input_file, 'r') as f:
            input_genes = list(dict.fromkeys(
                line.strip() for line in f if line.strip()
            ))

        existing_genes = {gene[0] for gene in GENE_LIST}

        for gene in input_genes:
            if gene in existing_genes:
                continue

            if gene in CHANGED_NAME or gene in UNIDENTIFIABLE_LIST:
                continue

            neighbors = find_neighbor(gene)

            if neighbors == -1:
                UNIDENTIFIABLE_LIST.append(gene)
            elif isinstance(neighbors, str):
                # Gene name was corrected
                CHANGED_NAME[gene] = neighbors
                gene = neighbors
                neighbors = find_neighbor(gene)
                if neighbors not in (-1, gene):
                    GENE_LIST.append((gene, neighbors))
            else:
                GENE_LIST.append((gene, neighbors))

            GENE_LIST.sort()
            write_database()

        logger.info(f"Processed {len(input_genes)} input genes")

    except Exception as e:
        logger.error(f"Error parsing input: {e}")


def download_svg(gene_list: List[str]) -> None:
    """Download SVG visualization from STRING."""
    if len(gene_list) < 2:
        logger.warning("Need at least 2 genes for SVG visualization")
        return

    driver = setup_driver()
    try:
        wait = WebDriverWait(driver, WAIT_TIMEOUT)

        driver.get(STRING_BASE_URL)
        wait.until(EC.element_to_be_clickable((By.ID, "search"))).click()
        wait.until(EC.element_to_be_clickable((By.ID, "multiple_identifiers"))).click()

        gene_input = wait.until(EC.presence_of_element_located(
            (By.ID, "primary_input:multiple_identifiers")
        ))
        gene_input.send_keys("\n".join(gene_list))

        wait.until(EC.presence_of_element_located(
            (By.ID, "species_text_multiple_identifiers")
        )).send_keys(SPECIES)

        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='input_form_multiple_identifiers']/div[5]/a")
        )).click()

        time.sleep(10)

        page_source = driver.page_source
        if "appear to match your input" in page_source:
            wait.until(EC.element_to_be_clickable(
                (By.XPATH, "//*[@id='proceed_form']/div[1]/div/div[2]/a[3]")
            )).click()

        time.sleep(20)
        
        # Configure display options
        wait.until(EC.element_to_be_clickable(
            (By.ID, "bottom_page_selector_table")
        )).click()
        time.sleep(5)

        wait.until(EC.element_to_be_clickable(
            (By.ID, "bottom_page_selector_settings")
        )).click()
        time.sleep(5)

        wait.until(EC.element_to_be_clickable(
            (By.ID, "confidence")
        )).send_keys(" ")
        time.sleep(10)

        wait.until(EC.element_to_be_clickable(
            (By.ID, "block_structures")
        )).send_keys(" ")
        time.sleep(10)

        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='standard_parameters']/div/div[1]/div[5]/a")
        )).click()
        time.sleep(15)

        wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='bottom_page_selector_legend']")
        )).click()
        time.sleep(10)

        wait.until(EC.element_to_be_clickable(
            (By.ID, "bottom_page_selector_table")
        )).click()
        time.sleep(25)

        # Download SVG
        download_button = wait.until(EC.element_to_be_clickable(
            (By.XPATH, "//*[@id='bottom_page_selector_table_container']/div/div[2]/div/div[3]/div[2]/a")
        ))
        download_button.click()
        time.sleep(30)

        logger.info("SVG downloaded successfully")

    except Exception as e:
        logger.error(f"Error downloading SVG: {e}")
    finally:
        driver.quit()


def main() -> None:
    """Main execution function."""
    try:
        # Setup directories
        for dir_path in [INFO_DIR, INPUT_DIR, SVG_DIR]:
            dir_path.mkdir(exist_ok=True)

        # Initialize files
        for file_path in [GENE_DATABASE_FILE, UNIDENTIFIABLE_GENE_FILE, CHANGED_NAME_GENE_FILE]:
            file_path.touch(exist_ok=True)

        # Read existing data
        read_database()
        read_unidentifiable()
        read_changed_names()

        # Process input genes
        parse_input()

        # Classify genes
        initialize_connections()
        identify_group_a()
        identify_group_b()

        # Write results
        write_gene_groups()
        write_intermediate_genes()
        write_unidentifiable()
        write_changed_names()

        # Collect all genes
        all_genes = (
            get_list_for_group("A") +
            get_list_for_group("B") +
            get_list_for_group("C") +
            get_list_for_group("D")
        )

        # Download visualization
        if all_genes:
            download_svg(all_genes)

        logger.info("Gene mapping download completed successfully")

    except Exception as e:
        logger.error(f"Fatal error in main: {e}")
        raise


if __name__ == "__main__":
    main()
