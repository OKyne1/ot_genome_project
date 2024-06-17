import argparse
import os
from Bio import SeqIO

def extract_gene_annotations(genbank_file, output_file):
    with open(genbank_file, "r") as input_handle, open(output_file, "w") as output_handle:
        # Write the header line
        output_handle.write("Gene\tStart\tEnd\tStrand\tProduct\tLocus Tag\tNote\n")
        
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "gene" or feature.type == "CDS":
                    gene_info = feature.qualifiers.get("gene", ["Unknown Gene"])[0]
                    start = feature.location.start.position + 1  # Adding 1 to convert from 0-based to 1-based
                    end = feature.location.end.position  # End position is inclusive, so no need to add 1
                    strand = feature.location.strand
                    product = feature.qualifiers.get("product", ["none"])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", ["none"])[0]
                    note = feature.qualifiers.get("note", ["none"])[0]
                    
                    annotation = f"{gene_info}\t{start}\t{end}\t{strand}\t{product}\t{locus_tag}\t{note}\n"
                    output_handle.write(annotation)

def main():
    parser = argparse.ArgumentParser(description="Convert a GenBank file into a text file containing gene annotations.")
    parser.add_argument("genbank_file", help="The path to the input GenBank file")
    args = parser.parse_args()
    
    genbank_file = args.genbank_file
    base_name = os.path.splitext(genbank_file)[0]
    output_file = f"{base_name}.txt"
    
    extract_gene_annotations(genbank_file, output_file)
    print(f"Gene annotations have been written to {output_file}")

if __name__ == "__main__":
    main()
