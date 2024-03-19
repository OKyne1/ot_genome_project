#!/bin/bash

# Need to work out how to do the relative paths to these different scripts

# Step 1
## genbank to fasta format
python3 1_gen_2_faa.py [script arguments]
## The script will be called using: bash domain_annotations.sh <genbank files>
## outputs .faa files

mkdir faa_files
mv *.faa faa_files/

# Step 2
## Domain finding
bash [path]/2_hmmersearch.sh [path]/hmm_files/*.hmm faa_files/*.faa
## output: .txt files, structure of name will always be gbffname_hmmname.txt. hmm name can either be tpr or ank.

mkdir hmmer_outputs
mv *.txt hmmer_outputs/

# Step 3
## Hmmer output parsing, only accepts 1 input and output must be specified 
## Code should now take any number of input txt files and specify dom based on these
python3 [path]/3_parsing_hmmer.py hmmer_outputs/*.txt
## output: txt file

## Check where the outputs are written to. This is probably right, but I'm not certain.
mv *_parsed.txt hmmer_outputs/

# Step 4
## gbff parsing
## Should take all txt files present with _parsed at the end
## Should use all input gbff files
## This will be more efficient if gbff files are matched to txt files (optional - depends how lazy I am)
## Else can just do all vs all
## Modify the code so it determines the -dom from the file name
python3 [path]/4_overwriting_gbff.py [input arguments] hmmer_outputs/*_parsed.txt -dom [tpr/ank]
## outputs: overwrites the gbff file - may want it to make a new modified version?


