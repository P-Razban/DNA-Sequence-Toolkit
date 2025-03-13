from functools import total_ordering
from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import argparse
import logging
import sqlite3
import re

class MiRDBHandler:
    def __init__(self, fasta_path, organism, threshold):
        self.url = 'http://www.mirdb.org/custom.html'
        self.fasta_sequences = self._load_fasta(fasta_path)
        self.organism = organism
        self.query_type = 'mRNA Target Sequence'
        self.threshold = threshold

    def _load_fasta(self, fasta_path) -> list:
        try:
            with open(fasta_path, 'r') as fasta_file:
                return list(SeqIO.parse(fasta_file, 'fasta'))
        except IOError:
            logging.exception("Error reading the provided FASTA file.")
            raise

class WebDriverHandler:
    def __init__(self, show_browser):
        self.driver = self._initialize_driver(show_browser)

    def _initialize_driver(self, show_browser):
        driver_path = 'chromedriver.exe'
        options = Options()
        if not show_browser:
            options.add_argument('--headless')
        return webdriver.Chrome(options=options, executable_path=driver_path)

    def select_option(self, element_name, option_text):
        dropdown = Select(self.driver.find_element(By.NAME, element_name))
        dropdown.select_by_visible_text(option_text)

    def input_sequence(self, sequence):
        input_box = self.driver.find_element(By.NAME, 'customSub')
        input_box.send_keys(sequence)

    def proceed_to_results(self):
        try:
            submit_button = self.driver.find_element(By.XPATH, '/html/body/table[2]/tbody/tr/td[3]/form/table/tbody/tr[5]/td/input[1]')
            submit_button.click()
        except Exception:
            logging.exception("Error while proceeding to results.")

    def wait_for_results(self, timeout=30):
        try:
            result_ready = EC.presence_of_element_located((By.XPATH, '/html/body/form/input[2]'))
            WebDriverWait(self.driver, timeout).until(result_ready)
            self.driver.find_element(By.XPATH, '/html/body/form/input[2]').click()
        except TimeoutException:
            logging.exception("Timeout waiting for results.")

class HTMLScraper:
    def __init__(self):
        self.soup = None

    def parse_html(self, page_source):
        self.soup = BeautifulSoup(page_source, 'html.parser')

    def filter_by_score(self, cutoff):
        matching_rows = []
        table_rows = self.soup.find('table', id='table1').find('tbody').find_all('tr')
        for idx, row in enumerate(table_rows[1:], start=1):
            try:
                columns = row.find_all('td')
                if int(columns[2].text) >= cutoff:
                    matching_rows.append(idx)
            except (AttributeError, ValueError):
                logging.exception("Error parsing score column.")
        return matching_rows

    def extract_details(self):
        try:
            table_data = self.soup.find_all('td')
            score = table_data[7].text if table_data[7].text.isdigit() else table_data[9].text
            locations = [int(x) for x in table_data[9].text.replace(' ', '').split(',')]
            return int(score), locations
        except (AttributeError, ValueError):
            logging.exception("Error extracting details.")
            return None, None

    def extract_seed_info(self):
        try:
            seeds = self.soup.find_all('font', {'color': '#0000FF'})
            clean_seeds = [re.sub(r"^\d+\s|\s\d+\s|\s\d+$", "", s.text.strip()) for s in seeds]
            return len(seeds), clean_seeds
        except AttributeError:
            logging.exception("Error extracting seed information.")
            return 0, []

    def extract_mirna_name(self):
        try:
            return self.soup.find_all('a', href=True)[1].font.text
        except AttributeError:
            logging.exception("Error extracting miRNA name.")
            return None

class ResultDatabase:
    def __init__(self, db_path):
        self.connection = sqlite3.connect(db_path)
        self._create_table()

    def _create_table(self):
        query = (
            """
            CREATE TABLE IF NOT EXISTS targets (
                sequence_id TEXT,
                score INTEGER,
                seed_count INTEGER,
                mirna_name TEXT,
                clean_seeds BLOB,
                locations BLOB
            );
            """
        )
        self._execute_query(query)

    def _execute_query(self, query, params=None):
        with self.connection:
            cursor = self.connection.cursor()
            if params:
                cursor.execute(query, params)
            else:
                cursor.execute(query)

    def insert_result(self, sequence_id, score, seed_count, mirna_name, clean_seeds, locations):
        query = (
            """
            INSERT INTO targets (sequence_id, score, seed_count, mirna_name, clean_seeds, locations)
            VALUES (?, ?, ?, ?, ?, ?);
            """
        )
        self._execute_query(query, (sequence_id, score, seed_count, mirna_name, str(clean_seeds), str(locations)))

    def export_to_csv(self, output_file):
        data_frame = pd.read_sql_query("SELECT * FROM targets", self.connection)
        data_frame.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="miRDB Automation Tool")
    parser.add_argument('fasta', type=str, help="Input FASTA file")
    parser.add_argument('output', type=str, help="Output CSV file")
    parser.add_argument('species', type=str, choices=["Human", "Mouse", "Rat", "Chicken", "Dog"], help="Target species")
    parser.add_argument('-c', '--cutoff', type=int, default=80, help="Score cutoff")
    parser.add_argument('-v', '--visible', action='store_true', help="Show browser during execution")
    args = parser.parse_args()

    logging.basicConfig(filename="mirdb_log.log", format="%(asctime)s - %(message)s", level=logging.INFO)

    mirdb_handler = MiRDBHandler(args.fasta, args.species, args.cutoff)
    driver_handler = WebDriverHandler(args.visible)
    scraper = HTMLScraper()
    results_db = ResultDatabase(f"{args.output}.db")

    try:
        for idx, record in enumerate(mirdb_handler.fasta_sequences):
            if 100 <= len(record.seq) <= 30000:
                driver_handler.driver.get(mirdb_handler.url)
                driver_handler.select_option('searchSpecies', mirdb_handler.organism)
                driver_handler.select_option('subChoice', mirdb_handler.query_type)
                driver_handler.input_sequence(str(record.seq))
                driver_handler.proceed_to_results()
                driver_handler.wait_for_results()

                page_source = driver_handler.driver.page_source
                scraper.parse_html(page_source)

                for row_index in scraper.filter_by_score(mirdb_handler.threshold):
                    details_buttons = driver_handler.driver.find_elements(By.NAME, '.submit')
                    details_buttons[row_index].click()

                    detail_page = driver_handler.driver.page_source
                    scraper.parse_html(detail_page)

                    score, locations = scraper.extract_details()
                    seed_count, clean_seeds = scraper.extract_seed_info()
                    mirna_name = scraper.extract_mirna_name()

                    results_db.insert_result(record.id, score, seed_count, mirna_name, clean_seeds, locations)
                    driver_handler.driver.back()
            else:
                logging.warning(f"Sequence {record.id} is out of length range.")

        results_db.export_to_csv(args.output)
    finally:
        driver_handler.driver.quit()
