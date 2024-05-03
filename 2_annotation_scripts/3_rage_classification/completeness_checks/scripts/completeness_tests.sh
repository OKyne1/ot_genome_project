#!/bin/bash
# What needs to be in the environment? blast+ and biopython...
# Activate the correct environment
script_dir=$(dirname "$(realpath "$0")")
# Make fna file from the gbff, give it the same name
python3 1_gen_2_faa.py "$@"
python3	2_faa_2_fna.py *.faa
for faa_file in *.faa; do
# Generate txt file name by replacing .fna with .txt
    txt_file="${faa_file%.faa}.txt"
    # Run blast for each .fna file
    blastn -query "$faa_file" -db "$script_dir"/rage_complete_db/complete_rage -out "$txt_file" -outfmt 6 -num_alignments 1 -evalue 1e-30
done
python3 "$script_dir"/1_blast_processing.py *.txt
# need to make a way to match up these things
for gbff_file in *.gbff; do
    txt_file="${gbff_file%.fna}_blast_complete.txt"
    # Run the Python script with the gbff and txt file names
    python3 "$script_dir"/2_formating_names.py "$gbff_file" "$txt_file"
done
rm *.fna *_blast.txt *_blast_complete.txt
echo "The script seems to be finished. Hopefully it worked."
# Input files: gbff_files db_path_and_name 
# genome.gbff --> genome.fna
# genome.fna + output_name.txt --> genome_blast.txt
# genome_blast.txt --> genome_blast_complete.txt
# genome.gbff genome_completeness_checked.gbff genome_blast_complete.txt --> genome_completeness_checked.gbff
############################# I think this should work but the blast search needs to be modified for faa files 
############################# My alignment stuff isn't going to work with this either
############################# Should I just translate the files to fna files? -> probably very easy
## Need to test the conversion file and the script as a whole.