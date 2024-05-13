import sys
import os
from Bio import SeqIO

def parse_genbank(genbank_files):
    for genbank_file in genbank_files:
        output_file = os.path.splitext(os.path.basename(genbank_file))[0] + ".txt"
        output_path = os.path.join(os.getcwd(), output_file)
        with open(output_path, 'w') as f_out:
            # Write header
            f_out.write("Gene\tProduct\tLocus\tStart\tEnd\tComplement\n")
            
            # Parse the GenBank file
            for record in SeqIO.parse(genbank_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        gene = feature.qualifiers.get("gene", ["NA"])[0]
                        product = feature.qualifiers.get("product", ["NA"])[0]
                        locus = feature.qualifiers.get("locus_tag", ["NA"])[0]
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        # Check if the feature location is on the complement strand
                        complement = "Yes" if feature.location.strand == -1 else "No"
                        # Write to file
                        f_out.write(f"{gene}\t{product}\t{locus}\t{start}\t{end}\t{complement}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python gbff_2_txt.py <genbank_file1> <genbank_file2> ...")
        sys.exit(1)
    
    genbank_files = sys.argv[1:]
    parse_genbank(genbank_files)