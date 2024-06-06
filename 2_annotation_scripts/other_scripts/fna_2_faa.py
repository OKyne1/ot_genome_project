from Bio.Seq import Seq
from Bio import SeqIO
import os

def translate_fna_to_faa(fna_file, faa_file):
    with open(faa_file, 'w') as faa_output:
        for record in SeqIO.parse(fna_file, 'fasta'):
            protein_seq = record.seq.translate(to_stop=True)  # Translate nucleotide sequence to protein sequence
            faa_output.write(f">{record.id}\n{protein_seq}\n")

def convert_all_fna_to_faa(directory):
    for fna_file in os.listdir(directory):
        if fna_file.endswith('.fna'):
            fna_path = os.path.join(directory, fna_file)
            base_name = os.path.splitext(fna_file)[0]
            faa_file = f"{base_name}.faa"
            faa_path = os.path.join(directory, faa_file)
            translate_fna_to_faa(fna_path, faa_path)
            print(f"Converted {fna_path} to {faa_path}")

if __name__ == "__main__":
    directory = '.'  # Set to the current directory or specify the directory containing the .fna files
    convert_all_fna_to_faa(directory)
