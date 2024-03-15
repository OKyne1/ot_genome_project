# This scrip does processing on the hmmersearch output when the --domtblout command is used. However, when using --tblout this is not required.

import argparse  # Import the argparse module to handle command-line arguments
import os  # Import the os module for file operations

def parse_input_file(input_file):
    # Initialize data structures
    data = {}  # Dictionary to store values from column 1 mapped to values from query name
    unique_col1_values = set()  # Set to store unique values from column 1

    # Open the input file for reading
    with open(input_file, 'r') as f:
        # Iterate over each line in the file
        for line in f:
            # Skip lines starting with '#' (comments)
            if not line.startswith('#'):
                # Split the line by whitespace
                parts = line.strip().split()
                # Ensure the line has enough columns
                if len(parts) >= 20:
                    # Extract values from columns 1 and query name
                    col1_value = parts[0]
                    col2_value = parts[3]  # Extract data from query name
                    # Add the value from column 1 to the set of unique values
                    unique_col1_values.add(col1_value)
                    # Store the value from query name in the dictionary,
                    # mapped to the value from column 1
                    if col1_value not in data:
                        data[col1_value] = col2_value
                    # If the value from column 1 is already in the dictionary,
                    # but with a different value from query name, raise an error
                    elif data[col1_value] != col2_value:
                        raise ValueError(f"Error: Different values found for {col1_value} in column 2")

    # Return the set of unique values from column 1 and the dictionary
    return unique_col1_values, data


def generate_output_file(unique_col1_values, data, input_file):
    # Generate the output file name by replacing '.txt' with '_count.txt'
    output_file = os.path.splitext(input_file)[0] + '_count.txt'
    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Iterate over the unique values from column 1
        for col1_value in unique_col1_values:
            # Retrieve the value from column 2 corresponding to the value from column 1
            col2_value = data.get(col1_value)
            # Calculate the number of times the value from column 1 appears in the input file
            col3_value = sum(1 for line in open(input_file) if col1_value in line)
            # Write the values to the output file, separated by tabs
            f.write(f"{col1_value}\t{col2_value}\t{col3_value}\n")


def main():
    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(description="Process .txt file and generate output in table format")
    # Add the input file argument
    parser.add_argument("input_files", nargs='+', help="Input .txt files")
    # Parse the command-line arguments
    args = parser.parse_args()

    try:
        # Iterate over each input file
        for input_file in args.input_files:
            # Parse the input file
            unique_col1_values, data = parse_input_file(input_file)
            # Generate the output file
            generate_output_file(unique_col1_values, data, input_file)
        # Print a success message
        print("Output files generated successfully.")
    # Catch any exceptions and print an error message
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
