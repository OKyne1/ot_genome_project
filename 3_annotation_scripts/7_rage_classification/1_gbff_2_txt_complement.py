# Import necessary modules
import sys
import os
from Bio import SeqIO

################################### creating temporary files to work with (instead of gbff) ######################################
# Usage: python gbff_2_txt_complement.py <genbank_file1> <genbank_file2> ...
# Environment: rage
# This script is to create a txt file with gene name, product, locus, start, end, complement and contig information
##################################################################################################################################

def parse_genbank(genbank_files):
    # Iterate over each GenBank file provided
    for genbank_file in genbank_files:
        # Create output file name based on the input file name
        output_file = os.path.splitext(os.path.basename(genbank_file))[0] + ".txt"
        output_path = os.path.join(os.getcwd(), output_file)
        
        # Open the output file for writing
        with open(output_path, 'w') as f_out:
            # Write header line to the output file
            f_out.write("Gene\tProduct\tLocus\tStart\tEnd\tComplement\tContig\n")
            
            # Parse the GenBank file using Biopython's SeqIO
            for record in SeqIO.parse(genbank_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        # Extract relevant information from the feature
                        gene = feature.qualifiers.get("gene", ["NA"])[0]
                        product = feature.qualifiers.get("product", ["NA"])[0]
                        locus = feature.qualifiers.get("locus_tag", ["NA"])[0]
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        # Check if the feature is on the complement strand
                        complement = "Yes" if feature.location.strand == -1 else "No"
                        # Get the contig ID
                        contig = record.id
                        # Write the extracted information to the output file
                        f_out.write(f"{gene}\t{product}\t{locus}\t{start}\t{end}\t{complement}\t{contig}\n")

if __name__ == "__main__":
    # Check if at least one GenBank file is provided
    if len(sys.argv) < 2:
        print("Usage: python gbff_2_txt.py <genbank_file1> <genbank_file2> ...")
        sys.exit(1)
    
    # Get the list of GenBank files from command line arguments
    genbank_files = sys.argv[1:]
    # Parse the provided GenBank files
    parse_genbank(genbank_files)