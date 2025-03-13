import os
import argparse
from Bio.Align.Applications import ClustalwCommandline

def run_msa(input_file, output_file, clustal_exe):
    """Perform Multiple Sequence Alignment (MSA) using ClustalW."""
    
    # Check if ClustalW executable exists
    if not os.path.isfile(clustal_exe):
        print(f"Error: ClustalW executable not found at {clustal_exe}")
        return

    # Create a ClustalW command line
    clustal = ClustalwCommandline(clustal_exe, infile=input_file, outfile=output_file)

    try:
        stdout, stderr = clustal()

        if stderr:
            print(f"Error running ClustalW: {stderr}")
        else:
            print(f"Alignment successful. Output saved in {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform multiple sequence alignment using ClustalW.")
    parser.add_argument("input_file", type=str, help="Path to input FASTA file")
    parser.add_argument("output_file", type=str, help="Path to save aligned sequences")
    parser.add_argument("clustal_exe", type=str, help="Path to ClustalW executable")

    args = parser.parse_args()

    run_msa(args.input_file, args.output_file, args.clustal_exe)