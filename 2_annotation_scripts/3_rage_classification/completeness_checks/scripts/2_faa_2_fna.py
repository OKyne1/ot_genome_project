import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def translate_faa_to_fna(faa_file):
    # Determine the output file name with .fna extension
    fna_file = faa_file.rsplit(".", 1)[0] + ".fna"
    # Open FAA file for reading
    with open(faa_file, "r") as faa_handle:
        # Open FNA file for writing
        with open(fna_file, "w") as fna_handle:
            # Iterate over each sequence record in the FAA file
            for record in SeqIO.parse(faa_handle, "fasta"):
                # Translate amino acid sequence to nucleotides
                nucleotide_seq = record.seq.back_transcribe()
                # Create a new SeqRecord with translated sequence
                nucleotide_record = SeqRecord(nucleotide_seq, id=record.id, description=record.description)
                # Write the translated sequence to the FNA file
                SeqIO.write(nucleotide_record, fna_handle, "fasta")
    print(f"Successfully converted {faa_file} to {fna_file}")

# Get command line inputs (FAA file paths)
faa_files = sys.argv[1:]

# Check if at least one FAA file is provided
if not faa_files:
    print("Please provide at least one FAA file.")
    sys.exit(1)

# Convert each FAA file to FNA
for faa_file in faa_files:
    translate_faa_to_fna(faa_file)
