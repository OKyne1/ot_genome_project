# Import necessary modules
import argparse
import os

################################################ rage derived regions ###############################################################################
#Matching Items in List Data:
#The region (line) from the table file must contain an item from the list_data in either the first or second column (cols[0] or cols[1]). The item comparison is case-insensitive, and any single quotes in the list items are ignored.
#
#No Matching Items in Exclusion Data:
### The region must not contain any item from the exclusion_data in either the first or second column. This exclusion check is also case-insensitive.
#
#Region Boundaries:
### If the above conditions are met, the code collects boundary information from the region and subsequent regions until it reaches a stopping condition:
### It skips lines that do not match any items in the list_data or that match items in the exclusion_data.
### It continues to collect boundaries from matching lines until either:
### A line is found that matches an exclusion item.
### A line is found that does not match any list item.
### The code reaches a limit of skipping one unmatched line before stopping the boundary collection.
#
#Sufficient Boundary Data:
#The boundaries list must contain more than four elements before it is considered valid for writing to the BED file. This ensures that both start and end positions are captured adequately.
####################################################################################################################################################

# Function to read the list file and return its contents as a list
def read_list_file(list_file):
    with open(list_file, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f]

# Function to process the table file and filter based on the list and exclusion data
def process_table_file(table_file, list_data, exclusion_data, output_file, name_value):
    boundaries = []
    boundaries.clear()  # Clear the boundaries list to ensure it's empty
    with open(table_file, 'r', encoding='utf-8') as f:
        prev_line_skipped = False
        for line in f:
            cols = line.strip().split('\t')  # Split the line into columns
            for item in list_data:
                item_lower = item.lower().replace("'", "")
                # Check if the item is in the first or second column
                if item_lower in cols[0].lower() or item_lower in cols[1].lower():
                    # Check for exclusions
                    if not any(exclusion.lower() in cols[0].lower() or exclusion.lower() in cols[1].lower() for exclusion in exclusion_data):
                        boundaries.append(cols[6])
                        boundaries.append(cols[3])
                        boundaries.append(cols[4])
                        skip_count = 0
                        next_line = f.readline().strip().split('\t')
                        while next_line:
                            if len(next_line) >= 2:
                                line_skipped = False
                                for next_item in list_data:
                                    next_item_lower = next_item.lower().replace("'", "")
                                    # Check if the next item is in the first or second column of the next line
                                    if next_item_lower in next_line[0].lower() or next_item_lower in next_line[1].lower():
                                        if not any(exclusion.lower() in next_line[0].lower() or exclusion.lower() in next_line[1].lower() for exclusion in exclusion_data):
                                            boundaries.append(next_line[4])
                                            skip_count = 0
                                        else:
                                            skip_count += 1
                                            line_skipped = True
                                        break
                                else:
                                    if skip_count < 1 and not line_skipped:
                                        skip_count += 1
                                    else:
                                        break
                                next_line = f.readline().strip().split('\t')
                            else:
                                break
                        if len(boundaries) > 4:
                            print(f"{boundaries[0]}\t{boundaries[1]}\t{boundaries[-1]}", file=output_file)
                        boundaries.clear()
                        prev_line_skipped = False
                    break
            else:
                prev_line_skipped = True

# Function to split the input file by the value in the 7th column
def split_file_by_col7(input_file):
    split_files = {}
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            col7_value = cols[6]  # Get the value of the 7th column
            if col7_value not in split_files:
                split_files[col7_value] = []
            split_files[col7_value].append(line)
    split_file_names = []
    for col7_value, lines in split_files.items():
        split_filename = f"{os.path.splitext(input_file)[0]}_{col7_value}.txt"
        split_file_names.append(split_filename)
        with open(split_filename, 'w', encoding='utf-8') as split_file:
            split_file.writelines(lines)  # Write the lines to the split file
    return split_file_names

# Function to combine the output files into a single file
def combine_output_files(output_files, combined_filename):
    with open(combined_filename, 'w', encoding='utf-8') as combined_file:
        for output_file in output_files:
            with open(output_file, 'r', encoding='utf-8') as f:
                combined_file.writelines(f.readlines())

# Main execution block
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process table and list files.")
    parser.add_argument("-tab", help="Table text file", nargs='+', required=True)
    parser.add_argument("-list", help="List of RAGE proteins", required=True)
    parser.add_argument("-exclusion", help="List of excluded proteins", required=True)
    args = parser.parse_args()

    list_data = read_list_file(args.list)
    exclusion_data = read_list_file(args.exclusion)

    for table_file in args.tab:
        split_files = split_file_by_col7(table_file)  # Split the table file
        output_files = []
        for split_file in split_files:
            output_filename = os.path.splitext(split_file)[0] + "_processed.bed"
            output_files.append(output_filename)
            with open(output_filename, 'w', encoding='utf-8') as rage:
                process_table_file(split_file, list_data, exclusion_data, rage, os.path.splitext(os.path.basename(split_file))[0])
        combined_filename = os.path.splitext(table_file)[0] + "_rage_derived.bed"
        combine_output_files(output_files, combined_filename)  # Combine the output files
        for split_file in split_files:
            os.remove(split_file)  # Remove the temporary split files
        for output_file in output_files:
            os.remove(output_file)  # Remove the temporary output files