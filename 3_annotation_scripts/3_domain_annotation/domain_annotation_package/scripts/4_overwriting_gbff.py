import argparse  # For parsing command-line arguments
import os  # For file path operations
from Bio import SeqIO  # For reading and writing GenBank files

#################################### Updating gbff with hmmer details ####################################################################
# Function to update product descriptions in GenBank files based on domain information from text files
# Usage: python 4_overwritting_gbff.py -gen genbank_file1.gbff genbank_file2.gbff ... -txt text_file1.txt text_file2.txt ...
# Environmet: rage
##########################################################################################################################################

def update_product_from_txt(genbank_file, txt_file, domain):
    # Read locus tags and product descriptions from the txt file into a dictionary
    locus_to_product = {}
    with open(txt_file, 'r') as txt_handle:
        for line in txt_handle:
            locus_tag, product = line.strip().split('\t')
            locus_to_product[locus_tag] = product

    # Update product descriptions and remove gene names in the GenBank file
    with open(genbank_file, 'r') as gbff_handle:
        records = list(SeqIO.parse(gbff_handle, 'genbank'))

    for record in records:
        for feature in record.features:
            if feature.type == 'gene' and 'gene' in feature.qualifiers and 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag'][0]
                if locus_tag in locus_to_product:
                    gene_names = feature.qualifiers['gene']
                    # Filter gene names based on domain-specific criteria
                    filtered_gene_names = [name for name in gene_names if name.lower() in ('trag', 'bamd')]
                    if filtered_gene_names:
                        feature.qualifiers['gene'] = filtered_gene_names
                    else:
                        del feature.qualifiers['gene']  # Remove gene names if filtered list is empty
            elif feature.type == 'CDS' and 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] in locus_to_product:
                gene_names = feature.qualifiers.get('gene', [])
                gene_names_lower = [name.lower() for name in gene_names] if gene_names else []
                if domain == 'tpr' and ('trag' in gene_names_lower or 'bamd' in gene_names_lower):
                    continue  # Skip modification for traG or bamD in the 'CDS' section
                    # traG and bamD are 2 proteins with TPRs and so would be overwritten without this. 
                    # Some other spurious annotations overwritten are pilW as to our knowledge it is absent in Orientia
                filtered_gene_names = [name for name in gene_names if name.lower() in ('trag', 'bamd')]
                if filtered_gene_names:
                    feature.qualifiers['gene'] = filtered_gene_names
                else:
                    if 'gene' in feature.qualifiers:
                        del feature.qualifiers['gene']  # Remove gene name if not 'trag' or 'bamd'

                locus_tag = feature.qualifiers['locus_tag'][0]
                if locus_tag in locus_to_product:
                    if 'product' in feature.qualifiers:  # Check if product is already set
                        feature.qualifiers['product'] = [locus_to_product[locus_tag]]  # Replace product
                    else:
                        feature.qualifiers['product'] = [locus_to_product[locus_tag]]  # Set product if not present

    # Write the updated records back to the GenBank file
    with open(genbank_file, 'w') as gbff_handle:
        SeqIO.write(records, gbff_handle, 'genbank')

# Function to determine domain type based on filename
def determine_domain_from_filename(filename):
    if 'ank' in filename:
        return 'ank'
    elif 'tpr' in filename:
        return 'tpr'
    else:
        return 'ank'  # Default to 'ank' if domain not found

# Main section to parse command-line arguments and execute update process
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update product descriptions in GenBank files from text files.")
    parser.add_argument("-gen", "--genbank_files", nargs='+', help="Input GenBank file paths", required=True)
    parser.add_argument("-txt", "--text_files", nargs='+', help="Input text files with locus tags and product descriptions", required=True)
    args = parser.parse_args()

    # Iterate through each GenBank file and each text file
    for genbank_file in args.genbank_files:
        for txt_file in args.text_files:
            domain = determine_domain_from_filename(txt_file)  # Determine domain type from text file name
            update_product_from_txt(genbank_file, txt_file, domain)  # Update product descriptions in GenBank file
            print(f"Product descriptions and gene names updated successfully in {genbank_file} with {txt_file}. Domain used: {domain}")