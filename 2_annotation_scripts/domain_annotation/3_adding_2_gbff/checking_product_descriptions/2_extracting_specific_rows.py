import argparse

def find_matching_rows(file1, file2, output_file):
    # Read locus tags from the first file into a set
    with open(file1, 'r') as f1:
        locus_tags = set(line.strip() for line in f1)

    # Find matching rows in file2 and write them to the output file
    with open(file2, 'r') as f2, open(output_file, 'w') as fout:
        # Write the header from file2 to the output file
        fout.write(next(f2))
        
        # Iterate through the rows in file2
        for line in f2:
            # Extract the locus tag from the current row
            locus_tag = line.split()[0]
            # Check if the locus tag is in the set from file1
            if locus_tag in locus_tags:
                # Write the matching row to the output file
                fout.write(line)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Find matching rows in file2 based on locus tags from file1.")
    parser.add_argument("file1", help="Path to the first text file with locus tags.")
    parser.add_argument("file2", help="Path to the second text file with rows to be filtered.")
    parser.add_argument("output_file", help="Path to the output text file.")
    args = parser.parse_args()

    # Call the function to find and write matching rows
    find_matching_rows(args.file1, args.file2, args.output_file)
