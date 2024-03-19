import argparse

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
    parser = argparse.ArgumentParser(description='Process input file and output a new file with two columns.')
    parser.add_argument('input_file', type=str, help='Input file name')
    parser.add_argument('output_file', type=str, help='Output file name')
    parser.add_argument('-dom', '--domain', choices=['ank', 'tpr'], default='ank', help='Domain type for column 2')

    args = parser.parse_args()

    process_input(args.input_file, args.output_file, args.domain)
