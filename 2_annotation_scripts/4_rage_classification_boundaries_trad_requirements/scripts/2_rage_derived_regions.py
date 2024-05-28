import argparse
import os

# Function to read the list file and return its contents as a list
def read_list_file(list_file):
    with open(list_file, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f]

def process_table_file(table_file, list_data, exclusion_data, output_file, name_value):
    boundaries = []
    boundaries.clear()
    with open(table_file, 'r', encoding='utf-8') as f:
        prev_line_skipped = False
        for line in f:
            cols = line.strip().split('\t')
            for item in list_data:
                item_lower = item.lower().replace("'", "")
                if item_lower in cols[0].lower() or item_lower in cols[1].lower():
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

def split_file_by_col7(input_file):
    split_files = {}
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            col7_value = cols[6]
            if col7_value not in split_files:
                split_files[col7_value] = []
            split_files[col7_value].append(line)
    split_file_names = []
    for col7_value, lines in split_files.items():
        split_filename = f"{os.path.splitext(input_file)[0]}_{col7_value}.txt"
        split_file_names.append(split_filename)
        with open(split_filename, 'w', encoding='utf-8') as split_file:
            split_file.writelines(lines)
    return split_file_names

def combine_output_files(output_files, combined_filename):
    with open(combined_filename, 'w', encoding='utf-8') as combined_file:
        for output_file in output_files:
            with open(output_file, 'r', encoding='utf-8') as f:
                combined_file.writelines(f.readlines())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process table and list files.")
    parser.add_argument("-tab", help="Table text file", nargs='+', required=True)
    parser.add_argument("-list", help="List of RAGE proteins", required=True)
    parser.add_argument("-exclusion", help="List of excluded proteins", required=True)
    args = parser.parse_args()

    list_data = read_list_file(args.list)
    exclusion_data = read_list_file(args.exclusion)

    for table_file in args.tab:
        split_files = split_file_by_col7(table_file)
        output_files = []
        for split_file in split_files:
            output_filename = os.path.splitext(split_file)[0] + "_processed.bed"
            output_files.append(output_filename)
            with open(output_filename, 'w', encoding='utf-8') as rage:
                process_table_file(split_file, list_data, exclusion_data, rage, os.path.splitext(os.path.basename(split_file))[0])
        combined_filename = os.path.splitext(table_file)[0] + "_rage_derived.bed"
        combine_output_files(output_files, combined_filename)
        for split_file in split_files:
            os.remove(split_file)
        for output_file in output_files:
            os.remove(output_file)
