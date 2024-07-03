# Function: This code processes bakta outputs into an .faa file with the locus tag and coding sequence. This was generated for subsequent domain name annotation.

import re
import argparse
import os

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
    locus_translations = []  

    with open(input_file, 'r') as f:
        cds_info = []  
        translation_started = False  
        for line in f:
            if line.strip().startswith("/locus_tag="):
                locus_tag = re.search(r'="(.+?)"', line.strip()).group(1)
            elif line.strip().startswith("/translation="):
                translation_started = True  
                translation = re.search(r'="(.+)', line.strip()).group(1)
                while not translation.endswith('"'):
                    line = next(f).strip()  
                    translation += line  
                translation = translation.strip('"')  
                cds_info.append(f'>{locus_tag}')
                cds_info.append(f'{translation}')
                locus_translations.append('\n'.join(cds_info))
                cds_info = []  
                translation_started = False  
            elif translation_started:
                translation += " " + line.strip()  

    return locus_translations  

def save_to_faa(input_file):
    """
    Saves locus_tag and translation information to a .faa file.

    Args:
        input_file (str): Path to the input .gbff file.
    """
    output_file = os.path.splitext(input_file)[0] + ".faa"  
    locus_translations = extract_locus_translation(input_file)
    with open(output_file, 'w') as f:
        for cds_info in locus_translations:
            f.write(cds_info + '\n')  
    print(f"Extraction completed. Output saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract locus_tag and translation for each CDS from GBFF files")
    parser.add_argument("input_files", nargs='+', help="Input .gbff files")
    args = parser.parse_args()

    for input_file in args.input_files:
        save_to_faa(input_file)

if __name__ == "__main__":
    main()
