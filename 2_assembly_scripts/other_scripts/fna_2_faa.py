import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def convert_fna_to_faa(input_file, output_file):
    # Open FNA file
    with open(input_file, "r") as fna_file:
        # Open FAA file for writing
        with open(output_file, "w") as faa_file:
            # Parse FNA file and write FAA file
            for record in SeqIO.parse(fna_file, "fasta"):
                # Translate nucleotide sequence to amino acid sequence
                amino_seq = record.seq.translate()
                
                # Write to FAA file
                faa_file.write(">" + record.id + "\n")
                faa_file.write(str(amino_seq) + "\n")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert FNA file to FAA file")
    parser.add_argument("input_file", help="Input FNA file path")
    parser.add_argument("output_file", help="Output FAA file path")
    args = parser.parse_args()

    # Convert FNA to FAA
    convert_fna_to_faa(args.input_file, args.output_file)