import sys
import os
# Usage: python3 blast_processing.py blast_outputs.txt (any number)

# Currently I've used a slightly dirty method to ensure their is no double counting. We want to know the amount of positions aligned. 
# My method uses the alignment of sequences and then using the amount of these sequences aligned as a percentage of the query (excluding bits of the subject which are longer than the query and aligned).
# Functionally this makes little difference and is deemed acceptable. But ideally, the number of positions alignmed would be calcated instead.

def process_blast_output(blast_file):
    # Extracting this info from a fasta: awk '/^>/ {header = "\"" substr($0, 2) "\""; getline; print header ": " length($0)","}' all_full_tra_proteins.fna # code used to make the dict
    # Define the gene lengths
    gene_lengths = {
        "Boryong_01064_SpoT_Full": 712,
        "GILLIAM_01018_SpoT_Full": 712,
        "Ikeda_00298_SpoT_Full": 712,
        "Karp_01184_SpoT_Full": 712,
        "KATO_01297_SpoT_Full": 712,
        "TA686_02046_SpoT_Full": 712,
        "UT176_01039_SpoT_Full": 713,
        "UT76HP_01102_SpoT_Full": 711
    }

    gene_matches = {}

    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            gene_name_match = columns[1]
            gene_name_output = columns[0]
            start = int(columns[6])
            end = int(columns[7])

            # Ensure start is less than or equal to end
            if start > end:
                start, end = end, start

            # Get the gene length from the dictionary
            gene_length = gene_lengths.get(gene_name_match)
            if gene_length is None:
                continue

            # Store the old end value
            old_end = end

            # Ensure the end value does not exceed the gene length
            if end > gene_length:
                end = gene_length

            # Ensure start is within valid range
            if start < 1:
                start = 1

            # Print old and new end values if they are different
            if old_end != end:
                print(f"For gene '{gene_name_output}', adjusted end value from {old_end} to {end}")

            # Create a list of zeros of the gene length
            gene_list = [0] * gene_length

            # Mark the matching positions as 1s
            for i in range(start - 1, end):
                gene_list[i] = 1

            # Calculate per_match and completeness
            matches = gene_list.count(1)
            per_match = matches / gene_length
            completeness = 'complete' if per_match >= 0.95 else 'truncated'

            # Update gene_matches dictionary
            if gene_name_output in gene_matches:
                # Modify existing entry
                existing_entry = gene_matches[gene_name_output]
                existing_list = existing_entry['matches']
                for i in range(start - 1, end):
                    if i < len(existing_list):
                        existing_list[i] = 1
                existing_entry['matches'] = existing_list
                existing_entry['per_match'] = existing_list.count(1) / gene_length
                existing_entry['completeness'] = 'complete' if existing_entry['per_match'] >= 0.95 else 'truncated'
            else:
                # Create new entry
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
