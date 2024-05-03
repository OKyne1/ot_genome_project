import argparse
import os

def process_file(input_files):
    for input_file in input_files:
        with open(input_file, 'r') as f:
            lines = f.readlines()

        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = f"{base_name}_rage_boundaries.bed"

        with open(output_file, 'w') as f_out:
            line_index = 0
            while line_index < len(lines):
                boundaries = []
                start_integrase = None
                start_dnaa = None
                end_integrase = None
                end_dnaa = None
                start_direction_complement = None
                end_direction_complement = None

                restart_loop = False
                while not restart_loop and line_index < len(lines):
                    line = lines[line_index]
                    parts = line.strip().split('\t')
                    gene = parts[0].lower()
                    product = parts[1].lower()

                    # Process the line with the current important values
                    # Modify the important values as needed
                    if ('integrase' in gene or 'integrase' in product) and (start_dnaa is None and start_integrase is None):
                        start_integrase = True
                        boundaries.append(parts[3])
                        boundaries.append(parts[4])
                        if parts[5].lower() == 'yes':
                            start_direction_complement = True
                        else:
                            start_direction_complement = False

                    elif ('dnaa' in gene or 'dnaa' in product) and (start_dnaa is None and start_integrase is None):
                        start_dnaa = True
                        boundaries.append(parts[3])
                        boundaries.append(parts[4])

                    elif ('dnaa' in gene or 'dnaa' in product) and (start_dnaa or start_integrase):
                        end_dnaa = True
                        if start_integrase and not start_direction_complement:
                            boundaries.append(parts[3])
                            boundaries.append(parts[4])
                            f_out.write(f"{base_name}\t{boundaries[0]}\t{boundaries[-1]}\n")
                            start_integrase = None
                            start_dnaa = None
                            end_integrase = None
                            end_dnaa = None
                            start_direction_complement = None
                            end_direction_complement = None
                            boundaries.clear()
                        elif start_integrase and start_direction_complement:
                            # Modify important values here
                            restart_loop = True
                        elif start_dnaa:
                            # Modify important values here
                            restart_loop = True

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
                                f_out.write(f"{base_name}\t{boundaries[0]}\t{boundaries[-1]}\n")
                                start_integrase = None
                                start_dnaa = None
                                end_integrase = None
                                end_dnaa = None
                                start_direction_complement = None
                                end_direction_complement = None
                                boundaries.clear()
                            elif not end_direction_complement:
                                # Modify important values here
                                restart_loop = True

                        elif start_integrase:
                            if start_direction_complement and end_direction_complement:
                                boundaries.append(parts[3])
                                boundaries.append(parts[4])
                                f_out.write(f"{base_name}\t{boundaries[2]}\t{boundaries[-1]}\n")
                                start_integrase = None
                                start_dnaa = None
                                end_integrase = None
                                end_dnaa = None
                                start_direction_complement = None
                                end_direction_complement = None
                                boundaries.clear()
                            elif not start_direction_complement and not end_direction_complement:
                                f_out.write(f"{base_name}\t{boundaries[0]}\t{boundaries[-1]}\n")
                                # Modify important values here
                                restart_loop = True
                            elif not start_direction_complement and end_direction_complement:
                                boundaries.append(parts[3])
                                boundaries.append(parts[4])
                                f_out.write(f"{base_name}\t{boundaries[0]}\t{boundaries[-1]}\n")
                                start_integrase = None
                                start_dnaa = None
                                end_integrase = None
                                end_dnaa = None
                                start_direction_complement = None
                                end_direction_complement = None
                                boundaries.clear()
                            elif start_direction_complement and not end_direction_complement:
                                # Modify important values here
                                restart_loop = True

                    elif start_dnaa or start_integrase:
                        boundaries.append(parts[3])
                        boundaries.append(parts[4])

                    # Move to the next line
                    line_index += 1

                # If restart_loop is True, reset the line index to reprocess the current line
                if restart_loop:
                    line_index -= 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process files and generate output")
    parser.add_argument("input_files", nargs="+", help="Input files")
    args = parser.parse_args()

    process_file(args.input_files)
