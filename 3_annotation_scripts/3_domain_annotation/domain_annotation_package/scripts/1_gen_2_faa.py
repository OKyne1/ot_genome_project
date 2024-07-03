import re  # Regular expressions for pattern matching
import argparse  # For handling command-line arguments
import os  # For file path manipulations

#################################### Generate an faa file for processing ###################################################
# Function: This code processes bakta outputs into an .faa file with the locus tag and coding sequence.
# usage: python 1_gen_2_faa.py input_file1.gbff input_file2.gbff ... input_fileN.gbff
############################################################################################################################

def extract_locus_translation(input_file):
    """
    Extracts locus_tag and translation information for each CDS from a .gbff file.

    Args:
        input_file (str): Path to the input .gbff file.

    Returns:
        list: List of strings containing locus_tag and translation information for each CDS.
    """
    locus_translations = []  # List to store extracted locus_tag and translation info

    with open(input_file, 'r') as f:
        cds_info = []  # Temporary list to store current CDS information
        translation_started = False  # Flag to indicate if translation sequence is being read
        for line in f:
            if line.strip().startswith("/locus_tag="):  # Check for locus_tag line
                # Extract locus_tag value using regex
                locus_tag = re.search(r'="(.+?)"', line.strip()).group(1)
            elif line.strip().startswith("/translation="):  # Check for translation start line
                translation_started = True  # Set flag to true
                # Extract the beginning of the translation sequence
                translation = re.search(r'="(.+)', line.strip()).group(1)
                while not translation.endswith('"'):  # Continue reading until the end of the translation sequence
                    line = next(f).strip()  # Read next line and strip whitespace
                    translation += line  # Append to translation sequence
                translation = translation.strip('"')  # Remove the ending quote
                cds_info.append(f'>{locus_tag}')  # Append locus_tag in FASTA format
                cds_info.append(f'{translation}')  # Append translation sequence
                # Join locus_tag and translation and add to the list
                locus_translations.append('\n'.join(cds_info))
                cds_info = []  # Reset CDS information list
                translation_started = False  # Reset flag
            elif translation_started:
                # Continue appending to translation sequence if the flag is set
                translation += " " + line.strip()

    return locus_translations  # Return the list of extracted locus_tag and translation info

def save_to_faa(input_file):
    """
    Saves locus_tag and translation information to a .faa file.

    Args:
        input_file (str): Path to the input .gbff file.
    """
    # Generate output filename by replacing the input file extension with .faa
    output_file = os.path.splitext(input_file)[0] + ".faa"
    # Extract locus_tag and translation information
    locus_translations = extract_locus_translation(input_file)
    # Write the extracted information to the output .faa file
    with open(output_file, 'w') as f:
        for cds_info in locus_translations:
            f.write(cds_info + '\n')  # Write each CDS info to the file
    # Print completion message with output file path
    print(f"Extraction completed. Output saved to {output_file}")

def main():
    # Set up argument parser for command-line arguments
    parser = argparse.ArgumentParser(description="Extract locus_tag and translation for each CDS from GBFF files")
    parser.add_argument("input_files", nargs='+', help="Input .gbff files")  # Accept multiple input files
    args = parser.parse_args()  # Parse arguments

    # Process each input file
    for input_file in args.input_files:
        save_to_faa(input_file)  # Extract and save locus_tag and translation info

if __name__ == "__main__":
    main()  # Call the main function when the script is executed