import argparse
import os
from Bio import SeqIO

def read_dbxref_file(dbxref_file):
    """Reads the dbxref values from the given text file."""
    with open(dbxref_file, 'r') as file:
        dbxref_values = {line.strip() for line in file}
    return dbxref_values

def get_script_directory():
    """Returns the directory containing the script."""
    return os.path.dirname(os.path.realpath(__file__))

def process_genbank_file(gb_file, dbxref_values):
    """Process the given GenBank file to remove 'gene' qualifier if any db_xref matches."""
    modified_records = []

    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if 'db_xref' in feature.qualifiers:
                for dbxref in feature.qualifiers['db_xref']:
                    if any(value in dbxref for value in dbxref_values):
                        if 'locus_tag' in feature.qualifiers:
                            locus_tag = feature.qualifiers['locus_tag'][0]
                            print(f"Removing 'gene' from feature with locus tag '{locus_tag}' in {gb_file}")
                            del feature.qualifiers['gene']
        modified_records.append(record)

    # Save the modified records back to a new GenBank file
    output_file = f"modified_{gb_file}"
    with open(output_file, "w") as output_handle:
        SeqIO.write(modified_records, output_handle, "genbank")
    print(f"Processed {gb_file}. Output saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Process GenBank files to remove gene qualifiers based on db_xref values.")
    parser.add_argument("genbank_files", nargs='+', help="GenBank files to process.")
    args = parser.parse_args()

    dbxref_file = os.path.join(get_script_directory(), "cog_list.txt")
    dbxref_values = read_dbxref_file(dbxref_file)

    for gb_file in args.genbank_files:
        process_genbank_file(gb_file, dbxref_values)

if __name__ == "__main__":
    main()
