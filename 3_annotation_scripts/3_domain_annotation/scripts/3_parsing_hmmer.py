import argparse  # For handling command-line arguments
import os  # For file path manipulations

###################################################### hmmer parsing #####################################################
# Function: This code processes hmmer output and identifies the locus_tag for genes which are identified as ank or tpr
# It then generates a simplified txt file with the locus tag and domain description for the hmmer search
# Usage: python3 3_parsing_hmmer.py hmmer_outputs/*.txt
##########################################################################################################################

def process_input(input_file, output_file, domain):
    """
    Processes input file to generate output file with two columns based on domain.

    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
        domain (str): Domain type ('ank' or 'tpr') based on input file name.
    """
    unique_entries = set()  # To store unique entries in column 1
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):  # Skip lines starting with '#'
                continue
            columns = line.strip().split()  # Split line into columns
            if len(columns) >= 12:
                col1 = columns[0]  # First column value
                col11_value = columns[11]  # Value from column 12 (assuming zero-indexed)
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
                outfile.write(f"{col1}\t{col2}\n")  # Write col1 and col2 to output file

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

        # Process input file to generate output file
        process_input(input_file, output_file, domain)
