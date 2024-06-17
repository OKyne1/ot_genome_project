#!/bin/bash
# Activate the correct environment
# Requires biopython, blast and their dependents. 
# The exact versions used can be found in the yml file, however, it should just work in an environment with these installed.
# Usage: bash path/gene_completeness_master_script.sh *.gbff (currently this MUST be run from the directory containing the gbff files)
#script_dir=$(dirname "$(realpath "$0")") # Stopped using this, as dropbox messed up my file paths
script_dir2=$(dirname "$BASH_SOURCE") # Had to use the relative path not the real path as dropbox is f**king up the paths.
echo "$script_dir2"
current_dir=$(pwd)
echo "$current_dir"
# Make faa file from the gbff, give it the same name
python3 "$script_dir2"/1_gen_2_faa.py "$@"
#python3	"$script_dir"/2_faa_2_fna.py *.faa
for faa_file in *.faa; do
# Generate txt file name by replacing .fna with .txt
    txt_file="${faa_file%.faa}.txt"
    # Run blast for each .fna file
    blastp -query "$faa_file" -db "$script_dir2"/rage_complete_db/aa_db/complete_rage_aa -out "$current_dir"/"$txt_file" -outfmt 6 -num_alignments 1 -evalue 1e-30
    # working script: blastp -query boryong.faa -db ../../scripts/rage_complete_db/aa_db/complete_rage_aa -out ./boryong.txt -outfmt 6 -num_alignments 1 -evalue 1e-30
done
python3 "$script_dir2"/2_blast_processing.py *.txt
# need to make a way to match up these things
for gbff_file in *.gbff; do
    txt_file="${gbff_file%.gbff}_complete.txt"
    # Run the Python script with the gbff and txt file names
    python3 "$script_dir2"/3_writing_to_gbff.py "$gbff_file" "$txt_file"
done
echo "The script seems to be finished. Hopefully it worked."
# do some file cleaning up