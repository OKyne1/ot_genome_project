from Bio import SeqIO  # Import SeqIO module from Biopython for handling sequence files
import os  # Import os module for file path operations

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
            # Split each line of the text file by '|' delimiter
            locus_tag, percent_aligned, completeness, name = line.strip().split("|")
            if completeness == "complete":
                name_suffix = name.split("_")[-1]  # Extract the last part of the name after splitting by '_'
                locus_tags_to_modify[locus_tag] = name_suffix  # Store locus tag and name suffix in dictionary

    # Read GenBank records from input file
    with open(input_file, "r") as f:
        records = list(SeqIO.parse(f, "genbank"))

    # Process each record and update CDS features if locus tag matches
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"][0]  # Extract locus tag of CDS feature
                if locus_tag in locus_tags_to_modify:
                    if "product" in feature.qualifiers:
                        name_suffix = locus_tags_to_modify[locus_tag]  # Get name suffix from dictionary
                        # Append ', Ank{name_suffix}' to the existing product description
                        feature.qualifiers["product"][0] += f", Ank{name_suffix}"

    # Define output filename based on input file name
    output_filename = os.path.splitext(input_file)[0] + "_anked.gbff"
    
    # Write updated records to output GenBank file
    with open(output_filename, "w") as f:
        SeqIO.write(records, f, "genbank")

if __name__ == "__main__":
    import argparse

    # Command line argument parsing
    parser = argparse.ArgumentParser(description="Append ', complete gene_name' to product description of genes in a GenBank file.")
    parser.add_argument("input_file", help="Input GenBank file")
    parser.add_argument("locus_tag_file", help="Text file containing locus tags")
    args = parser.parse_args()

    # Call process_gbff function with input file and locus tag file
    process_gbff(args.input_file, args.locus_tag_file)