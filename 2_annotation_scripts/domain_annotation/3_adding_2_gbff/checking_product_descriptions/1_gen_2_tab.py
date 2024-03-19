# This python code takes a gbff file and outputs a table with locus tag, gene and product as column titles
# Example usage: python3 gen_2_tab.py ../../1_gbff_2_faa/boryong.gbff out.txt

import argparse
import os
from Bio import SeqIO

def print_gene_info(genbank_file, output_file):
    with open(output_file, 'w') as f:
        f.write("locus_tag\tgene\tproduct\n")
        record = SeqIO.read(genbank_file, 'genbank')
        for feature in record.features:
            if feature.type == 'CDS':
                locus = feature.qualifiers.get('locus_tag', ['NA'])[0]
                gene = feature.qualifiers.get('gene', ['NA'])[0]
                product = feature.qualifiers.get('product', ['NA'])[0]
                f.write(f"{locus}\t{gene}\t{product}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Print locus tag, gene, and product of genes from GenBank files into a table format.")
    parser.add_argument("genbank_files", nargs='+', help="Input GenBank files")
    args = parser.parse_args()

    for genbank_file in args.genbank_files:
        output_file = os.path.splitext(genbank_file)[0] + ".txt"
        print_gene_info(genbank_file, output_file)
        print(f"Gene information extracted from {genbank_file} and saved to {output_file}")


