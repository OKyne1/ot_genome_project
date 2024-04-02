#!/bin/bash

# Usage: python gbff_2_txt.py <genbank_file1> <genbank_file2> ...

script_dir=$(dirname "$(realpath "$0")")

## genbank to txt table
python3 "$script_dir"/gbff_2_txt.py "$@"

## Changing the locations
mkdir gbff_tables
mv *.txt gbff_tables

## Define RAGE derived regions
python3 "$script_dir"/rage_derived_regions.py -tab gbff_tables/*.txt -list "$script_dir"/rage_proteins/lists_combined.txt -exclusion "$script_dir"/rage_proteins/exclusion_list.txt

mkdir bed_files
mv *_rage.bed bed_files 

# Script format
## 1. convert gbff to txt file
## 2. run the rage_id script
## 3. subtract 1 from the start values
## 4. then try scripts to identify complete RAGEs