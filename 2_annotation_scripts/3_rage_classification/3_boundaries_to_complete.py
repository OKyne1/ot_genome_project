import argparse
import os

# Inputs: .bed and xxx.txt
# Outputs: xxx_complete_RAGEs.bed

def parse_bed_file(bed_file):
    bed_dict = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            fields = line.strip().split('\t')
            bed_line = line.strip()  # Entire line as key
            region_id = fields[0]  # Assuming the first field is a unique identifier
            bed_dict[bed_line] = (region_id, int(fields[1]), int(fields[2]))
    return bed_dict

def parse_txt_file(txt_file, bed_dict):
    result_dict = {}
    with open(txt_file, 'r') as txt:
        for line in txt:
            fields = line.strip().split('\t')
            start_txt = int(fields[3])
            end_txt = int(fields[4])
            for bed_line, (region_id, start_bed, end_bed) in bed_dict.items():
                if start_bed <= start_txt <= end_bed and start_bed <= end_txt <= end_bed:
                    region = (fields[0], fields[1], fields[2], fields[3], fields[4], fields[5])
                    if bed_line in result_dict:
                        result_dict[bed_line].append(region)
                    else:
                        result_dict[bed_line] = [region]
    return result_dict

def check_files_for_strings(file_data, tra_genes, tran, trah, transposase, cargo):
    # Ensure strings are in lowercase for case-insensitive comparison
    lowercase_tra_genes = [s.lower() for s in tra_genes]
    lowercase_tran = [s.lower() for s in tran]
    lowercase_trah = [s.lower() for s in trah]
    lowercase_transposase = [s.lower() for s in transposase]
    lowercase_cargo = [s.lower() for s in cargo]
    
    passed_regions = []
    
    for file_key, file_content in file_data.items():
        # Initialize counters for each string
        tra_genes_counts = {s: 0 for s in lowercase_tra_genes}
        tran_counts = {s: 0 for s in lowercase_tran}
        trah_counts = {s: 0 for s in lowercase_trah}
        transposase_counts = {s: 0 for s in lowercase_transposase}
        cargo_counts = {s: 0 for s in lowercase_cargo}
        
        # Check for matches in the file content
        for entry in file_content:
            line_lower = ' '.join(map(str.lower, entry))
            for s in lowercase_tra_genes:
                if s in line_lower:
                    tra_genes_counts[s] += 1
            for s in lowercase_tran:
                if s in line_lower:
                    tran_counts[s] += 1
            for s in lowercase_trah:
                if s in line_lower:
                    trah_counts[s] += 1
            for s in lowercase_transposase:
                if s in line_lower:
                    transposase_counts[s] += 1
            for s in lowercase_cargo:
                if s in line_lower:
                    cargo_counts[s] += 1
        
        # Check if all strings are present at least once
        all_tra_genes_present = all(count > 0 for count in tra_genes_counts.values())
        
        # Check if at least one string is present in additional lists
        tran_present = any(count > 0 for count in tran_counts.values())
        trah_present = any(count > 0 for count in trah_counts.values())
        transposase_present = any(count > 0 for count in transposase_counts.values())
        cargo_present = any(count > 0 for count in cargo_counts.values())
        
        if all_tra_genes_present and tran_present and trah_present and transposase_present and cargo_present:
            passed_regions.append(file_key)
    
    return passed_regions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BED and TXT files.")
    parser.add_argument("bed_file", help="Input BED file")
    parser.add_argument("txt_file", help="Input TXT file")
    args = parser.parse_args()

    # Parse bed file
    bed_dict = parse_bed_file(args.bed_file)
    
    # Parse txt file
    result_dict = parse_txt_file(args.txt_file, bed_dict)
    
    # Example usage with dictionary data
    tra_genes = ['traa', 'trab', 'trac', 'trad', 'trae', 'traf', 'trag', 'trai', 'trak', 'tral', 'tran', 'trau', 'trav', 'traw']
    tran = ['tran\t', 'tran '] # Not contained in tra_genes as there are mishits otherwise
    trah = ['trah\t', 'trah '] # Not contained in tra_genes as there are mishits otherwise
    transposase = ['transposase', 'tpn', 'is630', 'is110', 'is5', 'isot6']
    cargo = ['spot', 'ppGpp', 'synthetase', 'hydrolase', 'methyltransferase', 'helicase', 'histidine kinase', 'mrp', 'atp-binding', 'endonuclease', 'hnh', 'membrane protein', 'ankyrin', 'ank', 'tetratricopeptide', 'tpr', 'hypothetical'] # Initially just going to take the proteins from Jeanne's paper - may add more from the list_combined file in the future
    
    passed_regions = check_files_for_strings(result_dict, tra_genes, tran, trah, transposase, cargo)
    
    # Extract base name of the TXT file and append "_complete_RAGEs.bed" to it
    output_file_name = os.path.splitext(os.path.basename(args.txt_file))[0] + "_complete_RAGEs.bed"
    
    # Write passed regions to the output file
    with open(output_file_name, "w") as bed_output:
        for region in passed_regions:
            bed_output.write(region + "\n")
