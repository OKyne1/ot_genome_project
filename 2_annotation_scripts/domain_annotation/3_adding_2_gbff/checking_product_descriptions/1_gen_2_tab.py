# This python code takes a gbff file and outputs a table with locus tag, gene and product as column titles
# Example usage: python3 gen_2_tab.py ../../1_gbff_2_faa/boryong.gbff out.txt

# Import necessary modules
import argparse  # For parsing command-line arguments
from Bio import SeqIO  # For handling sequence data in various formats

# Define a function to print gene information from a GenBank file into a table format
def print_gene_info(genbank_file, output_file):
    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Write the header row with column names
        f.write("locus_tag\tgene\tproduct\n")
        # Read the GenBank file and parse it
        record = SeqIO.read(genbank_file, 'genbank')
        # Loop through each feature in the GenBank record
        for feature in record.features:
            # Check if the feature type is CDS (coding sequence)
            if feature.type == 'CDS':
                # Extract locus tag, gene, and product information from the feature
                # If the information is not available, use 'NA' as default
                locus = feature.qualifiers.get('locus_tag', ['NA'])[0]
                gene = feature.qualifiers.get('gene', ['NA'])[0]
                product = feature.qualifiers.get('product', ['NA'])[0]
                # Write the gene information into the output file as a tab-separated row
                f.write(f"{locus}\t{gene}\t{product}\n")

# Entry point of the script
if __name__ == "__main__":
    # Create an argument parser object with a description
    parser = argparse.ArgumentParser(description="Print locus tag, gene, and product of genes from a GenBank file into a table format.")
    # Add command-line arguments for input and output file paths
    parser.add_argument("genbank_file", help="Input GenBank file path")
    parser.add_argument("output_file", help="Output text file path")
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the print_gene_info function with the provided input and output file paths
    print_gene_info(args.genbank_file, args.output_file)

