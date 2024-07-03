import re  # Regular expression library for pattern matching
import argparse  # Argument parsing library
import os  # Operating system functionality

#################################### Generate an faa file for processing ###################################################
# Function: This code processes bakta outputs into an .faa file with the locus tag and coding sequence.
# usage: python 1_gen_2_faa.py input_file1.gbff input_file2.gbff ... input_fileN.gbff
# Originally from the domain annotation part
############################################################################################################################

def extract_locus_translation(input_file):
    """
    Extracts locus_tag and translation information for each CDS from a .gbff file.

    Args:
        input_file (str): Path to the input .gbff file.

    Returns:
        list: List of strings containing locus_tag and translation information for each CDS.
    """
    locus_translations = []  # List to store locus_tag and translation information

    with open(input_file, 'r') as f:
        cds_info = []  # Temporary list to store each CDS information
        translation_started = False  # Flag to track if translation has started
        for line in f:
            if line.strip().startswith("/locus_tag="):
                locus_tag = re.search(r'="(.+?)"', line.strip()).group(1)  # Extract locus_tag using regex
            elif line.strip().startswith("/translation="):
                translation_started = True  # Set flag to True indicating translation has started
                translation = re.search(r'="(.+)', line.strip()).group(1)  # Extract translation using regex
                while not translation.endswith('"'):
                    line = next(f).strip()  # Read subsequent lines until end of translation
                    translation += line  # Append lines to translation
                translation = translation.strip('"')  # Remove surrounding quotes
                cds_info.append(f'>{locus_tag}')  # Add locus_tag as header
                cds_info.append(f'{translation}')  # Add translation sequence
                locus_translations.append('\n'.join(cds_info))  # Join header and sequence, add to list
                cds_info = []  # Clear cds_info for next CDS
                translation_started = False  # Reset translation flag
            elif translation_started:
                translation += " " + line.strip()  # Append lines to translation if translation has started

    return locus_translations  # Return list of locus_tag and translation information

def save_to_faa(input_file):
    """
    Saves locus_tag and translation information to a .faa file.

    Args:
        input_file (str): Path to the input .gbff file.
    """
    output_file = os.path.splitext(input_file)[0] + ".faa"  # Generate output filename with .faa extension
    locus_translations = extract_locus_translation(input_file)  # Extract locus_tag and translation information
    with open(output_file, 'w') as f:
        for cds_info in locus_translations:
            f.write(cds_info + '\n')  # Write locus_tag and translation information to .faa file
    print(f"Extraction completed. Output saved to {output_file}")  # Print completion message with output filename

def main():
    parser = argparse.ArgumentParser(description="Extract locus_tag and translation for each CDS from GBFF files")
    parser.add_argument("input_files", nargs='+', help="Input .gbff files")
    args = parser.parse_args()

    for input_file in args.input_files:
        save_to_faa(input_file)  # Process each input .gbff file and save locus_tag and translation to .faa file

if __name__ == "__main__":
    main()