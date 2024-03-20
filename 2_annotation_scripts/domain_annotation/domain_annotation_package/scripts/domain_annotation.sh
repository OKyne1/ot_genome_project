#!/bin/bash
# Usage: bash domain_annotation.sh <gbff_files>
# Specify the hmm files from the script?
# Need to work out how to do the relative paths to these different scripts
script_dir=$(dirname "$(realpath "$0")")
## genbank to fasta format
python3 "$script_dir"/1_gen_2_faa.py "$@"
## The script will be called using: bash domain_annotations.sh <genbank files>
## outputs .faa files
echo "script 1 complete"
mkdir faa_files
mv *.faa faa_files/
# Step 2 ## Domain finding
bash "$script_dir"/2_hmmersearch.sh "$script_dir"/../hmm_files/*.hmm faa_files/*.faa
## output: .txt files, structure of name will always be gbffname_hmmname.txt. hmm name can either be tpr or ank.
echo "script 2 complete"
mkdir hmmer_outputs
mv *.txt hmmer_outputs/
# Step 3 ## Hmmer output parsing, only accepts 1 input and output must be specified 
## Code should now take any number of input txt files and specify dom based on these
python3 "$script_dir"/3_parsing_hmmer.py hmmer_outputs/*.txt
## output: txt file
echo "script 3 complete"
## gbff modification
python3 "$script_dir"/4_overwriting_gbff.py -gen "$@" -txt hmmer_outputs/*_parsed.txt
echo "script 4 complete"
# To do:
## Add in the command line stuff for the overall script
## Test the script sequentially
## Determine the dependencies