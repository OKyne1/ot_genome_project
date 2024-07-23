import argparse
import os
import glob

################################################# rage boundaries ##########################################################
# Usage: python 2_rage_boundaries.py <input_file1> <input_file2> ...
# Environment: rage
# This script finds regions bounded by an integase (running into the region) and ends with either a dnase or an integrase
# The script is a bit more complicated as depending on the boundary type the last gene may or may not be included
# More information can be found written up in a READ.md file
# The boundaries identified are written to bed files.
############################################################################################################################

def split_file_by_col7(input_file):
    # Dictionary to store lines based on the value of the 7th column
    split_files = {}
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            col7_value = cols[6]  # Get the value of the 7th column
            if col7_value not in split_files:
                split_files[col7_value] = []
            split_files[col7_value].append(line)
    
    # List to store the names of split files
    split_file_names = []
    for col7_value, lines in split_files.items():
        # Create a filename based on the original filename and the 7th column value
        split_filename = f"{os.path.splitext(input_file)[0]}-{col7_value}-split.txt"
        split_file_names.append(split_filename)
        with open(split_filename, 'w', encoding='utf-8') as split_file:
            split_file.writelines(lines)  # Write the lines to the split file
    return split_file_names

def process_split_files(split_files):
    # List to store the names of output BED files
    output_files = []

    for split_file in split_files:
        base_name = os.path.splitext(os.path.basename(split_file))[0]
        # Extract the desired value between the first and second '-' in the filename
        extracted_value = base_name.split('-')[1]
        output_file = f"{base_name}-boundary.bed"
        output_files.append(output_file)

        with open(split_file, 'r') as f:
            lines = f.readlines()

        with open(output_file, 'w') as f_out:
            line_index = 0
            while line_index < len(lines):
                # Initialize variables for processing boundaries
                boundaries = []
                start_integrase = None
                start_dnaa = None
                end_integrase = None
                end_dnaa = None
                start_direction_complement = None
                end_direction_complement = None

                restart_loop = False
                while line_index < len(lines) and not restart_loop:
                    line = lines[line_index]
                    parts = line.strip().split('\t')
                    gene = parts[0].lower()
                    product = parts[1].lower()

                    # Check for integrase gene/product and set flags
                    if ('integrase' in gene or 'integrase' in product) and (start_dnaa is None and start_integrase is None):
                        start_integrase = True
                        boundaries.append(parts[3])
                        boundaries.append(parts[4])
                        if parts[5].lower() == 'yes':
                            start_direction_complement = True
                        else:
                            start_direction_complement = False

                    # Check for dnaa gene/product and set flags
                    elif ('dnaa' in gene or 'dnaa' in product) and (start_dnaa is None and start_integrase is None):
                        start_dnaa = True
                        boundaries.append(parts[3])
                        boundaries.append(parts[4])

                    # Process end dnaa gene/product
                    elif ('dnaa' in gene or 'dnaa' in product) and (start_dnaa or start_integrase):
                        end_dnaa = True
                        if start_integrase and not start_direction_complement:
                            boundaries.append(parts[3])
                            boundaries.append(parts[4])
                            f_out.write(f"{extracted_value}\t{boundaries[0]}\t{boundaries[-1]}\n")
                            # Reset variables
                            start_integrase = None
                            start_dnaa = None
                            end_integrase = None
                            end_dnaa = None
                            start_direction_complement = None
                            end_direction_complement = None
                            boundaries.clear()
                            restart_loop = True
                        elif start_integrase and start_direction_complement:
                            restart_loop = True
                        elif start_dnaa:
                            restart_loop = True

                    # Process end integrase gene/product
                    elif ('integrase' in gene or 'integrase' in product) and (start_dnaa or start_integrase):
                        end_integrase = True
                        if parts[5].lower() == 'yes':
                            end_direction_complement = True
                        else:
                            end_direction_complement = False 
                        if start_dnaa:
                            if end_direction_complement:
                                boundaries.append(parts[3])
                                boundaries.append(parts[4])
                                f_out.write(f"{extracted_value}\t{boundaries[0]}\t{boundaries[-1]}\n")
                                # Reset variables
                                start_integrase = None
                                start_dnaa = None
                                end_integrase = None
                                end_dnaa = None
                                start_direction_complement = None
                                end_direction_complement = None
                                boundaries.clear()
                                restart_loop = True
                            elif not end_direction_complement:
                                restart_loop = True

                        elif start_integrase:
                            if start_direction_complement and end_direction_complement:
                                boundaries.append(parts[3])
                                boundaries.append(parts[4])
                                f_out.write(f"{extracted_value}\t{boundaries[2]}\t{boundaries[-1]}\n")
                                # Reset variables
                                start_integrase = None
                                start_dnaa = None
                                end_integrase = None
                                end_dnaa = None
                                start_direction_complement = None
                                end_direction_complement = None
                                boundaries.clear()
                                restart_loop = True
                            elif not start_direction_complement and not end_direction_complement:
                                f_out.write(f"{extracted_value}\t{boundaries[0]}\t{boundaries[-1]}\n")
                                restart_loop = True
                            elif not start_direction_complement and end_direction_complement:
                                boundaries.append(parts[3])
                                boundaries.append(parts[4])
                                f_out.write(f"{extracted_value}\t{boundaries[0]}\t{boundaries[-1]}\n")
                                # Reset variables
                                start_integrase = None
                                start_dnaa = None
                                end_integrase = None
                                end_dnaa = None
                                start_direction_complement = None
                                end_direction_complement = None
                                boundaries.clear()
                                restart_loop = True
                            elif start_direction_complement and not end_direction_complement:
                                restart_loop = True

                    # Continue processing boundaries
                    elif start_dnaa or start_integrase:
                        boundaries.append(parts[3])
                        boundaries.append(parts[4])

                    line_index += 1

                if restart_loop:
                    line_index -= 1

    return output_files

def combine_bed_files(output_files):
    # Dictionary to store files to be combined based on the base name
    combined_files = {}
    
    for file in output_files:
        base_name = file.split('-')[0]
        if base_name not in combined_files:
            combined_files[base_name] = []
        combined_files[base_name].append(file)
    
    # Combine the contents of the files with the same base name
    for base_name, files in combined_files.items():
        with open(f"{base_name}_rage_boundaries.bed", 'w') as combined_file:
            for file in files:
                with open(file, 'r') as f:
                    combined_file.write(f.read())

def remove_temp_files(files):
    # Remove temporary files
    for file in files:
        if os.path.exists(file):
            os.remove(file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process files and generate output")
    parser.add_argument("input_files", nargs="+", help="Input files")
    args = parser.parse_args()

    split_files = []
    output_files = []

    # Split input files by the 7th column
    for input_file in args.input_files:
        split_files.extend(split_file_by_col7(input_file))

    # Process the split files
    output_files.extend(process_split_files(split_files))
    
    # Combine BED files with the same base name
    combine_bed_files(output_files)

    # Remove temporary split files and individual BED files
    remove_temp_files(split_files + output_files)