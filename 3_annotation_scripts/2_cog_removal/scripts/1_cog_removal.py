import argparse  # For handling command-line arguments
import os  # For handling file paths
from Bio import SeqIO  # For reading and writing sequence files in various formats

############################################## Removal of dodgy COG names #############################################################
# python script_name.py file1.gbff file2.gbff ... filen.gbff
# this outputs a files with this structure modified_(original_name).gbff
# This script removes the value in /gene= for genes with a /db_xref value matching the ones in the file cog_list.txt
# These were identified by 
### 1. identifying all the novel gene names given by bakta
### 2. looking for their support from the different db_xref values
### 3. Identification of the genes with low support/no support in the other db's but annotated from COG --> quite a few had this issue
# It is possible that this approach results in excessive gene name removal, my tests haven't shown this but it may be worth considering
# In the boryong genome this removes 134/1742 genes
#######################################################################################################################################


def read_dbxref_file(dbxref_file):
    """Reads the dbxref values from the given text file."""
    # Open the file and read all lines, stripping any whitespace, and store them in a set
    with open(dbxref_file, 'r') as file:
        dbxref_values = {line.strip() for line in file}
    return dbxref_values

def get_script_directory():
    """Returns the directory containing the script."""
    # Return the directory where the script is located
    return os.path.dirname(os.path.realpath(__file__))

def process_genbank_file(gb_file, dbxref_values):
    """Process the given GenBank file to remove 'gene' qualifier if any db_xref matches."""
    modified_records = []  # List to hold modified records

    # Parse the GenBank file
    for record in SeqIO.parse(gb_file, "genbank"):
        # Iterate through each feature in the record
        for feature in record.features:
            # Check if the feature has a 'db_xref' qualifier
            if 'db_xref' in feature.qualifiers:
                # Iterate through each 'db_xref' in the feature
                for dbxref in feature.qualifiers['db_xref']:
                    # Check if any value in dbxref_values matches the dbxref
                    if any(value in dbxref for value in dbxref_values):
                        # If there's a match, and the feature has a 'locus_tag', print a message and delete the 'gene' qualifier
                        if 'locus_tag' in feature.qualifiers:
                            locus_tag = feature.qualifiers['locus_tag'][0]
                            print(f"Removing 'gene' from feature with locus tag '{locus_tag}' in {gb_file}")
                            del feature.qualifiers['gene']
        # Add the modified record to the list
        modified_records.append(record)

    # Save the modified records back to a new GenBank file
    output_file = f"modified_{gb_file}"
    with open(output_file, "w") as output_handle:
        SeqIO.write(modified_records, output_handle, "genbank")
    print(f"Processed {gb_file}. Output saved to {output_file}")

def main():
    # Set up argument parser for command-line arguments
    parser = argparse.ArgumentParser(description="Process GenBank files to remove gene qualifiers based on db_xref values.")
    parser.add_argument("genbank_files", nargs='+', help="GenBank files to process.")
    args = parser.parse_args()

    # Determine the path to the dbxref file, which is located in the same directory as the script
    dbxref_file = os.path.join(get_script_directory(), "cog_list.txt")
    # Read the dbxref values from the file
    dbxref_values = read_dbxref_file(dbxref_file)

    # Process each GenBank file provided as a command-line argument
    for gb_file in args.genbank_files:
        process_genbank_file(gb_file, dbxref_values)

if __name__ == "__main__":
    main()  # Call the main function when the script is executed