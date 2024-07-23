import argparse
import os

# Inputs: .bed and xxx.txt
# Outputs: xxx_complete_RAGEs.bed

################################################## Identification of Complete RAGEs #################################################################
# Usage: python 3_boundaries_to_complete.py input_file.bed input_file.txt
# Environment: rage
#
# Uses regions defined by an input bed file
# Then checks if the corresponding regions in the txt file (generated from the gbff file) contain the required genes to be defined as a complete RAGE
#
# Possible modifications: 
# Gene requirements can be changed --> tra_genes is a list of genes that are required once in each complete RAGE (one day may want to add non-tra genes to this list)
# Number of traD required --> modify this "specific_count_required = 2" (though I can't see why you'd want to)
# Adding transposases or cargo genes --> modify transposase or cargo lists respectively
#####################################################################################################################################################


def parse_bed_file(bed_file):
    bed_dict = {}
    # Open the BED file and read line by line
    with open(bed_file, 'r') as bed:
        for line in bed:
            # Split each line into fields
            fields = line.strip().split('\t')
            bed_line = line.strip()  # Entire line as key
            region_id = fields[0]  # Assuming the first field is a unique identifier
            # Store the region ID and coordinates in the dictionary
            bed_dict[bed_line] = (region_id, int(fields[1]), int(fields[2]))
    return bed_dict

def parse_txt_file(txt_file, bed_dict):
    result_dict = {}
    # Open the TXT file and read line by line
    with open(txt_file, 'r') as txt:
        for line in txt:
            # Split each line into fields
            fields = line.strip().split('\t')
            start_txt = int(fields[3])
            end_txt = int(fields[4])
            # Compare TXT coordinates with BED coordinates
            for bed_line, (region_id, start_bed, end_bed) in bed_dict.items():
                if start_bed <= start_txt <= end_bed and start_bed <= end_txt <= end_bed:
                    region = (fields[0], fields[1], fields[2], fields[3], fields[4], fields[5])
                    if bed_line in result_dict:
                        result_dict[bed_line].append(region)
                    else:
                        result_dict[bed_line] = [region]
    return result_dict

def check_files_for_strings(file_data, tra_genes, transposase, cargo):
    # Ensure strings are in lowercase for case-insensitive comparison
    lowercase_tra_genes = [s.lower() for s in tra_genes]
    lowercase_transposase = [s.lower() for s in transposase]
    lowercase_cargo = [s.lower() for s in cargo]
    
    # Define the specific string and the required count
    specific_string = 'complete trad'.lower()
    # Seting the threshold of traD requirement as >=2 (due to the 2 traD genes required for a complete RAGE)
    specific_count_required = 2

    passed_regions = []
    
    for file_key, file_content in file_data.items():
        # Initialize counters for each string for this region
        tra_genes_counts = {s: 0 for s in lowercase_tra_genes}
        transposase_counts = {s: 0 for s in lowercase_transposase}
        cargo_counts = {s: 0 for s in lowercase_cargo}
        
        specific_string_count = 0

        # Check for matches in the file content
        for entry in file_content:
            line_lower = ' '.join(map(str.lower, entry))
            for s in lowercase_tra_genes:
                if s in line_lower:
                    tra_genes_counts[s] += 1
                    if s == specific_string:
                        specific_string_count += 1
            for s in lowercase_transposase:
                if s in line_lower:
                    transposase_counts[s] += 1
            for s in lowercase_cargo:
                if s in line_lower:
                    cargo_counts[s] += 1
        
        # Check if all tra_genes strings are present at least once
        all_tra_genes_present = all(count > 0 for count in tra_genes_counts.values())
        transposase_present = any(count > 0 for count in transposase_counts.values())
        cargo_present = any(count > 0 for count in cargo_counts.values())

        # Additional check for specific string count
        if all_tra_genes_present and transposase_present and cargo_present and specific_string_count >= specific_count_required:
            passed_regions.append(file_key)
    
    return passed_regions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BED and TXT files.")
    parser.add_argument("bed_file", help="Input BED file")
    parser.add_argument("txt_file", help="Input TXT file")
    args = parser.parse_args()

    # Parse BED file
    bed_dict = parse_bed_file(args.bed_file)
    
    # Parse TXT file
    result_dict = parse_txt_file(args.txt_file, bed_dict)
    
    # Gene lists
    # tra_genes required once in each complete rage (with the exception of traD which requires 2 copies)
    tra_genes = ['complete traa', 'complete trab', 'complete trac', 'complete trad', 'complete trae', 'complete traf', 'complete trag', 'complete trah', 'complete trai', 'complete trak', 'complete tral', 'complete tran', 'complete trau', 'complete trav', 'complete traw', 'complete integrase', 'complete trbc']
    # complete RAGEs require at least one of these transposases
    transposase = ['transposase', 'tpn', 'is630', 'is110', 'is5', 'isot6']
    # RAGE cargo genes, 1 is required per complete RAGE
    cargo = ['spot', 'ppGpp', 'synthetase', 'hydrolase', 'methyltransferase', 'helicase', 'histidine kinase', 'mrp', 'atp-binding', 'endonuclease', 'hnh', 'membrane protein', 'ankyrin', 'ank', 'tetratricopeptide', 'tpr', 'hypothetical']  # Initially just going to take the proteins from Jeanne's paper - may add more from the list_combined file in the future
    
    # Check if strings are present in the regions
    passed_regions = check_files_for_strings(result_dict, tra_genes, transposase, cargo)
    
    # Extract base name of the TXT file and append "_complete_RAGEs.bed" to it
    output_file_name = os.path.splitext(os.path.basename(args.txt_file))[0] + "_complete_RAGEs.bed"
    
    # Write passed regions to the output file
    with open(output_file_name, "w") as bed_output:
        for region in passed_regions:
            bed_output.write(region + "\n")
    print(f"Output written to {output_file_name}")