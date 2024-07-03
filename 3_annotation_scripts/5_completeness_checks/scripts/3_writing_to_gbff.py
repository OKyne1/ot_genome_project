from Bio import SeqIO
import os

######################################################## writing to gbff ############################################################
# Function: writing the complete anks to the gbff file
# Usage: python 3_writing_to_gbff.py input_file.gbff ank_file.txt
# Environment: rage
#####################################################################################################################################

def process_gbff(input_file, locus_tag_file):
    # Read locus tags and values from the text file
    locus_tags_to_modify = {}
    with open(locus_tag_file, "r") as f:
        for line in f:
            locus_tag, percent_aligned, completeness, name = line.strip().split("|")
            if completeness == "complete":
                name_suffix = name.split("_")[-1]
                locus_tags_to_modify[locus_tag] = name_suffix

    # Parse the input GenBank file using Biopython's SeqIO
    with open(input_file, "r") as f:
        records = list(SeqIO.parse(f, "genbank"))

    # Iterate over each record and each feature to modify 'product' descriptions
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"][0]
                # Check if the locus tag needs modification based on completeness
                if locus_tag in locus_tags_to_modify:
                    if "product" in feature.qualifiers:
                        name_suffix = locus_tags_to_modify[locus_tag]
                        # Append ", complete gene_name" to the 'product' description
                        feature.qualifiers["product"][0] += f", complete {name_suffix}"

    # Create an output filename by appending "_completeness_checked.gbff"
    output_filename = os.path.splitext(input_file)[0] + "_completeness_checked.gbff"

    # Write the modified records back to a new GenBank file
    with open(output_filename, "w") as f:
        SeqIO.write(records, f, "genbank")

if __name__ == "__main__":
    import argparse

    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Append ', complete gene_name' to product description of genes in a GenBank file.")
    parser.add_argument("input_file", help="Input GenBank file")
    parser.add_argument("locus_tag_file", help="Text file containing locus tags")
    args = parser.parse_args()

    # Call the function to process the GenBank file with locus tags
    process_gbff(args.input_file, args.locus_tag_file)
