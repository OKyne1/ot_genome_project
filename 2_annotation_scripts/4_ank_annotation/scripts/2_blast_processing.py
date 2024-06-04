import sys
import os

# Blast commands: blastn -query genome.fna -db path/to/rage_complete_db/complete_rage -out genome_name.txt -outfmt 6 -num_alignments 1 -evalue 1e-30

# Usage: python3 blast_processing.py blast_outputs.txt (any number)

# Currently I've used a slightly dirty method to ensure their is no double counting. We want to know the amount of positions aligned. 
# My method uses the alignment of sequences and then using the amount of these sequences aligned as a percentage of the query (excluding bits of the subject which are longer than the query and aligned).
# Functionally this makes little difference and is deemed acceptable. But ideally, the number of positions alignmed would be calcated instead.

import re

def process_blast_output(blast_file):
    gene_lengths = {
        "03": 214,
        "08": 389,
        "11": 228,
        "10": 551,
        "12": 494,
        "20": 508,
        "24": 936
    }

    gene_matches = {}

    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            match = re.search(r'Ank(\d{2})', columns[1])  # Find 2 digits following "Ank"
            if match:
                gene_name_match = match.group(1)
            else:
                continue  # If no match is found, skip this line

            gene_name_output = columns[0]
            start = int(columns[6])
            end = int(columns[7])

            gene_length = gene_lengths.get(gene_name_match)
            if gene_length is None:
                continue

            old_end = end

            if end > gene_length:
                end = gene_length

            if start < 1:
                start = 1

            if old_end != end:
                print(f"For gene '{gene_name_output}', adjusted end value from {old_end} to {end}")

            gene_list = [0] * gene_length

            for i in range(start - 1, end):
                gene_list[i] = 1

            matches = gene_list.count(1)
            per_match = matches / gene_length
            completeness = 'complete' if per_match >= 0.95 else 'truncated'

            if gene_name_output in gene_matches:
                existing_entry = gene_matches[gene_name_output]
                existing_list = existing_entry['matches']
                for i in range(start - 1, end):
                    if i < len(existing_list):
                        existing_list[i] = 1
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

    return gene_matches



def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py blast_output1.txt blast_output2.txt ...")
        sys.exit(1)

    for blast_file in sys.argv[1:]:
        gene_matches = process_blast_output(blast_file)
        
        # Extract the base filename without extension
        base_filename = os.path.splitext(os.path.basename(blast_file))[0]
        output_filename = f"{base_filename}_complete.txt"
        
        # Write the output to the file
        with open(output_filename, 'w') as output_file:
            for gene, data in gene_matches.items():
                output_file.write(f"{gene}|Percent Match: {data['per_match']}|{data['completeness']}|{data['matched']}\n")

if __name__ == "__main__":
    main()
