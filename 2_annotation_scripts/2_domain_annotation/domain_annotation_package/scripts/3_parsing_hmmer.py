# I modified this code to make it take as many inputs as i like, give output files names input_parsed.txt and domain should be specified based on the name of the input.

import argparse
import os

def process_input(input_file, output_file, domain):
    unique_entries = set()  # To store unique entries in column 1
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            columns = line.strip().split()
            if len(columns) >= 12:
                col1 = columns[0]
                col11_value = columns[11]
                if domain == 'ank':
                    col2 = f"Ankyrin repeat protein with {col11_value} ankyrin repeats"
                elif domain == 'tpr':
                    col2 = "Tetratricopeptide repeat protein"
                    # Check if col1 is already written to output
                    if col1 in unique_entries:
                        continue  # Skip writing duplicate rows
                    else:
                        unique_entries.add(col1)  # Add col1 to set of unique entries
                else:
                    raise ValueError("Invalid domain type specified")
                outfile.write(f"{col1}\t{col2}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process input file(s) and output new file(s) with two columns.')
    parser.add_argument('input_files', nargs='+', type=str, help='Input file(s) name(s)')

    args = parser.parse_args()

    for input_file in args.input_files:
        # Determine the domain type based on the input file name
        if input_file.endswith('ank.txt'):
            domain = 'ank'
        elif input_file.endswith('tpr.txt'):
            domain = 'tpr'
        else:
            raise ValueError(f"Could not determine domain type from input file name: {input_file}")

        # Generate output file name by appending _parsed.txt to the input file name
        output_file = os.path.splitext(input_file)[0] + '_parsed.txt'

        process_input(input_file, output_file, domain)

