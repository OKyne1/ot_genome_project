import argparse
import os

# Function to read the list file and return its contents as a list
def read_list_file(list_file):
    with open(list_file, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f]

def process_table_file(table_file, list_data, exclusion_data, output_file, name_value):
    # Initialize an empty list to store boundaries
    boundaries = []
    # Open the table file for reading
    with open(table_file, 'r', encoding='utf-8') as f:
        # Flag to track if the previous line was skipped
        prev_line_skipped = False
        
        # Loop through each line in the file
        for line in f:
            # Split the line into columns using tab as delimiter
            cols = line.strip().split('\t')

            # Check if any item from the list is present in the current line
            for item in list_data:
                # Convert item to lowercase and remove apostrophes
                item_lower = item.lower().replace("'", "")
                # Check if the item is present in either of the first two columns of the line
                if item_lower in cols[0].lower() or item_lower in cols[1].lower():
                    # Check if any exclusion item is present in the line
                    if not any(exclusion.lower() in cols[0].lower() or exclusion.lower() in cols[1].lower() for exclusion in exclusion_data):
                        # Add start and end boundaries to the boundaries list
                        boundaries.append(cols[3])
                        boundaries.append(cols[4])
                        # Counter to track the number of lines skipped within the region
                        skip_count = 0
                        # Read the next line and split it into columns
                        next_line = f.readline().strip().split('\t')
                        
                        # Loop to process subsequent lines within the same region
                        while next_line:
                            if len(next_line) >= 2:
                                # Flag to track if the current line has been skipped
                                line_skipped = False
                                # Check if any item from the list is present in the current line
                                for next_item in list_data:
                                    next_item_lower = next_item.lower().replace("'", "")
                                    # If the item is present, check if it should be excluded
                                    if next_item_lower in next_line[0].lower() or next_item_lower in next_line[1].lower():
                                        if not any(exclusion.lower() in next_line[0].lower() or exclusion.lower() in next_line[1].lower() for exclusion in exclusion_data):
                                            # Add end boundary to the boundaries list
                                            boundaries.append(next_line[4])
                                            # Reset skip_count if a match is found
                                            skip_count = 0
                                        else:
                                            # Increment skip_count and set line_skipped if line should be skipped
                                            skip_count += 1
                                            line_skipped = True
                                        # Break after finding the first match for the next item
                                        break
                                else:
                                    # If less than one line has been skipped and the line has not been skipped already, skip this line
                                    if skip_count < 1 and not line_skipped:
                                        skip_count += 1
                                    else:
                                        # Break if already skipped twice or line has been skipped
                                        break
                                # Read the next line and split it into columns
                                next_line = f.readline().strip().split('\t')
                            else:
                                # Break if next_line does not contain enough elements
                                break
                        
                        # If boundaries list contains more than three elements, print the region details to the output file
                        if len(boundaries) > 3:
                            print(f"{name_value}\t{boundaries[0]}\t{boundaries[-1]}", file=output_file)
                        # Clear the boundaries list for the next region
                        boundaries.clear()
                        # Reset prev_line_skipped since a match is found
                        prev_line_skipped = False
                    # Break after processing the current line
                    break
            else:
                # Set prev_line_skipped to True if the current line is skipped
                prev_line_skipped = True

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Process table and list files.")
    parser.add_argument("-tab", help="Table text file", nargs='+', required=True)
    parser.add_argument("-list", help="List text file", required=True)
    parser.add_argument("-exclusion", help="Exclusion list text file", required=True)
    args = parser.parse_args()

    # Read the list files
    list_data = read_list_file(args.list)
    exclusion_data = read_list_file(args.exclusion)

    for table_file in args.tab:
        # Create output filename based on table filename
        output_filename = os.path.splitext(os.path.basename(table_file))[0] + "_rage.bed"

        # Open the rage file for writing
        with open(output_filename, 'w', encoding='utf-8') as rage:
            # Process the table file
            process_table_file(table_file, list_data, exclusion_data, rage, os.path.splitext(os.path.basename(table_file))[0])
