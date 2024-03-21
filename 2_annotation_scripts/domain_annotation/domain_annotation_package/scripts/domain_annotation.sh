#!/bin/bash

# Cleaned up the comments - I'm not sure if this will add characters breaking the code - check next run.
# May want to make it cluster compatable

# Usage: bash domain_annotation.sh <gbff_files>
script_dir=$(dirname "$(realpath "$0")")

## genbank to fasta format
python3 "$script_dir"/1_gen_2_faa.py "$@"
## outputs .faa files

echo "script 1 complete"

mkdir faa_files
mv *.faa faa_files/

# Step 2 
## Domain finding
bash "$script_dir"/2_hmmersearch.sh "$script_dir"/../hmm_files/*.hmm faa_files/*.faa
## output: .txt files, structure of name will always be gbffname_hmmname.txt. hmm name can either be tpr or ank.

echo "script 2 complete"

mkdir hmmer_outputs
mv *.txt hmmer_outputs/

# Step 3 
## Code should takes any number of input txt files and specify dom based on these
python3 "$script_dir"/3_parsing_hmmer.py hmmer_outputs/*.txt
## output: txt file

echo "script 3 complete"

# Step 4
## gbff modification
python3 "$script_dir"/4_overwriting_gbff.py -gen "$@" -txt hmmer_outputs/*_parsed.txt

echo "script 4 complete"
