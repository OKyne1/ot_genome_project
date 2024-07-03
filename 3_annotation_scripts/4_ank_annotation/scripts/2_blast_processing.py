import sys
import os

####################################################### blast processing #####################################################################
# Usage: python3 blast_processing.py blast_outputs.txt (any number)
# The method uses the alignment of sequences and then using the amount of these sequences aligned as a percentage of the query 
# This exclude bits of the subject which are longer than the query and aligned but in tests this didn't affect the results
# Functionally this makes little difference and is deemed acceptable. But ideally, the number of positions alignmed would be calcated instead.
##############################################################################################################################################

import re  # Regular expression module for pattern matching

def process_blast_output(blast_file):
    # Dictionary containing gene lengths for specific genes
    gene_lengths = {
        "03": 214,
        "08": 389,
        "11": 228,
        "10": 551,
        "12": 494,
        "20": 508,
        "24": 936
    }

    gene_matches = {}  # Dictionary to store gene matches and their attributes

    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')  # Split each line into columns based on tab delimiter
            match = re.search(r'Ank(\d{2})', columns[1])  # Find 2 digits following "Ank" in the second column
            if match:
                gene_name_match = match.group(1)  # Extract the matched gene name
            else:
                continue  # If no match is found, skip processing this line

            gene_name_output = columns[0]  # First column contains the gene name
            start = int(columns[6])  # Start position of the alignment
            end = int(columns[7])    # End position of the alignment

            gene_length = gene_lengths.get(gene_name_match)  # Get the length of the gene from the dictionary
            if gene_length is None:
                continue  # If gene length is not defined, skip processing this gene

            old_end = end

            if end > gene_length:
                end = gene_length  # Adjust end position if it exceeds gene length

            if start < 1:
                start = 1  # Adjust start position if it is less than 1

            if old_end != end:
                print(f"For gene '{gene_name_output}', adjusted end value from {old_end} to {end}")

            gene_list = [0] * gene_length  # Initialize a list to mark positions aligned

            for i in range(start - 1, end):
                gene_list[i] = 1  # Mark positions aligned within the gene sequence

            matches = gene_list.count(1)  # Count the number of positions aligned
            per_match = matches / gene_length  # Calculate percentage of positions aligned
            completeness = 'complete' if per_match >= 0.95 else 'truncated'  # Determine completeness based on alignment percentage

            if gene_name_output in gene_matches:
                existing_entry = gene_matches[gene_name_output]
                existing_list = existing_entry['matches']
                for i in range(start - 1, end):
                    if i < len(existing_list):
                        existing_list[i] = 1  # Update existing alignment marks
                existing_entry['matches'] = existing_list
                existing_entry['per_match'] = existing_list.count(1) / gene_length
                existing_entry['completeness'] = 'complete' if existing_entry['per_match'] >= 0.95 else 'truncated'
            else:
                gene_matches[gene_name_output] = {
                    'matches': gene_list,
                    'per_match': per_match,
                    'completeness': completeness,
                    'matched': gene_name_match
                }

    return gene_matches  # Return dictionary containing processed gene matches

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py blast_output1.txt blast_output2.txt ...")
        sys.exit(1)  # Exit program with error if no blast output files provided

    for blast_file in sys.argv[1:]:
        gene_matches = process_blast_output(blast_file)  # Process each blast output file
        
        # Extract the base filename without extension
        base_filename = os.path.splitext(os.path.basename(blast_file))[0]
        output_filename = f"{base_filename}_complete.txt"  # Create output filename
        
        # Write the output to the file
        with open(output_filename, 'w') as output_file:
            for gene, data in gene_matches.items():
                output_file.write(f"{gene}|Percent Match: {data['per_match']}|{data['completeness']}|{data['matched']}\n")  # Write gene information to output file

if __name__ == "__main__":
    main()