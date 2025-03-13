# DNA-Sequence-Toolkit
A Python toolkit for extracting, processing, aligning, and comparing DNA sequences, including miRNA extraction from mirDB.
# 🧬 DNA Sequence Toolkit

This repository contains a Python toolkit for **extracting, processing, aligning, and comparing DNA sequences**, including **miRNA extraction from mirDB**.

## 🚀 Features
✔ Extract DNA subsequences from FASTA files
✔ Combine multiple FASTA files into one
✔ Perform multiple sequence alignment (MSA) using ClustalW
✔ Compare aligned sequences and analyze variations
✔ Extract miRNA data from miRDB

🛠️ Installation & Requirements
🔹 Dependencies
Ensure you have the following installed:

Python 

    version 3

Packages and libraries

    pandas
    beautifulsoup4
    Biopython
    ClustalW2 (for alignment)
    geckodriver
    openpyxl
    selenium

chromedriver (miRDB extraction)

🔹 Installation
1. Clone the repository

Copy
git clone https://github.com/P-Razban/DNA-Sequence-Toolkit.git

cd DNA-Sequence-Toolkit

2. Install dependencies

Copy
pip install -r requirements.txt

(Also you can manually install Biopython, Pandas and other requrments:)
Copy
pip install biopython pandas


🧬 Usage Guide
1️⃣ Extract & Combine DNA Sequences
🔹 Extract a subsequence from a FASTA file

Copy
python extract_combine_fasta.py --input input.fasta --output extract_subsequence.fasta --start 1 --end 100

🔹 Combine multiple FASTA files into one

Copy
python extract_combine_fasta.py --combine folder_path --output combine_fasta_files.fasta

2️⃣ Perform Multiple Sequence Alignment (MSA)

Copy
python align_sequences.py --input combine_fasta_files.fasta --output aligned.fasta --clustal_path "C:\Program Files (x86)\ClustalW2\clustalw2.exe"


3️⃣ Compare Aligned DNA Sequences

Copy
python compare_sequences.py --csv input_data.csv --align aligned.fasta --output results.csv

4️⃣ Extract miRNA Data from miRDB

Copy
python mirdb_search.py --mirna mirna_list.fasta --output mirdb_results.csv {organism} -c cutoff -v visible
organism = Human,Rat,Mouse,Chicken,Dog
default cutoff=80
-v = Show browser during execution

📂 Project Structure

DNA-Sequence-Toolkit/
│── scripts/
│   ├── mirdb_search.py        # Extract miRNA from mirDB website
│   ├── extract_combine_fasta.py  # Extract & combine FASTA sequences
│   ├── align_sequences.py     # Perform multiple sequence alignment
│   ├── compare_sequences.py   # Compare aligned sequences
│── README.md
│── requirements.txt  # List dependencies (Biopython, Pandas, etc.)
│── .gitignore        # Ignore unnecessary files

📜 License
This project is open-source and available under the MIT License.