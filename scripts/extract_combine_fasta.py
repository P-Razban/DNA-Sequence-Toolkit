import os
import argparse
from Bio import SeqIO

def extract_subsequence(input_file, output_fasta_file, start_pos, end_pos):
    """Extracts subsequences from a FASTA file and writes them to an output file."""
    with open(output_fasta_file, "w") as f_out:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            f_out.write(f">{seq_record.description}\n")
            subsequence = seq_record.seq[start_pos:end_pos]
            f_out.write(f"{subsequence}\n\n")
    print(f"Subsequences extracted and saved to {output_fasta_file}.")

def combine_fasta_files(input_directory, combined_output_file):
    """Combines multiple FASTA files from a directory into a single output FASTA file."""
    with open(combined_output_file, 'w') as combined_file:
        for filename in sorted(os.listdir(input_directory)):
            file_path = os.path.join(input_directory, filename)
            if os.path.isfile(file_path) and filename.endswith(('.fasta', '.fa')):
                try:
                    with open(file_path, 'r') as input_file:
                        combined_file.write(input_file.read() + "\n")
                except Exception as e:
                    print(f"Error reading file {filename}: {e}")
    print(f"Combining completed. The combined file is saved as {combined_output_file}.")

def main():
    parser = argparse.ArgumentParser(description="FASTA Sequence Processing Tool")
    subparsers = parser.add_subparsers(dest="command", help="Choose an operation")
    
    # Extract Subsequence Parser
    extract_parser = subparsers.add_parser("extract", help="Extract subsequences from a FASTA file")
    extract_parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    extract_parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    extract_parser.add_argument("-s", "--start", type=int, required=True, help="Start position (1-based index)")
    extract_parser.add_argument("-e", "--end", type=int, required=True, help="End position (1-based index)")
    
    # Combine FASTA Files Parser
    combine_parser = subparsers.add_parser("combine", help="Combine multiple FASTA files into one")
    combine_parser.add_argument("-d", "--directory", required=True, help="Directory containing FASTA files")
    combine_parser.add_argument("-o", "--output", required=True, help="Output combined FASTA file")
    
    args = parser.parse_args()
    
    if args.command == "extract":
        extract_subsequence(args.input, args.output, args.start, args.end)
    elif args.command == "combine":
        combine_fasta_files(args.directory, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()