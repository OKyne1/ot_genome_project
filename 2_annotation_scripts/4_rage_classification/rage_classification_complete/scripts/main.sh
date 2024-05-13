#!/bin/bash

# Usage: python gbff_2_txt.py <genbank_file1> <genbank_file2> ...

script_dir=$(dirname "$BASH_SOURCE")

## genbank to txt table
python3 "$script_dir"/1_gbff_2_txt_complement.py "$@"

## Changing the locations
mkdir gbff_tables
mv *.txt gbff_tables

for file in gbff_tables/*.txt; do sed -i '1d' "$file"; done

## Define RAGE derived regions
python3 "$script_dir"/2_rage_derived_regions.py -tab gbff_tables/*.txt -list "$script_dir"/rage_proteins/lists_combined.txt -exclusion "$script_dir"/rage_proteins/exclusion_list.txt

mkdir bed_files
#mv *_rage_derived.bed bed_files 

# Define potential RAGE boundaries
python3 "$script_dir"/2_rage_boundaries.py gbff_tables/*.txt

#mv *_range_boundaries.bed bed_files

# RAGE boundary validation

# Directory containing the BED files
#bed_dir="bed_files"

# Iterate through the directory
for file1 in *_range_boundaries.bed; do
    # Extract the base filename (without extension)
    base_name="${file1%_range_boundaries.bed}"
    
    # Form the corresponding file2 name
    file2="${base_name}_rage_derived.bed"

    # Check if the corresponding file2 exists
    if [ -e "$file2" ]; then
        # Run bedtools intersect for the pair and save output to a file
        output_file="${base_name}_validated_boundaries.bed"
        bedtools intersect -a "$file1" -b "$file2" -f 1 -wa > "$output_file"
        # echo "Output saved to $output_file"
    # else
    #     echo "Corresponding file not found for $file1"
    fi
done

#!/bin/bash

# Directories containing the files
gbff_dir="gbff_tables"
bed_dir="bed_files"

# Iterate through the gbff_tables directory
for gbff_file in "$gbff_dir"/*.txt; do
    # Extract the base filename (without extension) from gbff_file
    base_name=$(basename "$gbff_file" .txt)
    
    # Form the corresponding filename in bed_files
    bed_file="$bed_dir/${base_name}_validated_boundaries.bed"

    # Check if the corresponding bed_file exists
    if [ -e "$bed_file" ]; then
        # Run the Python script for the pair
        python3 "$script_dir"/3_boundaries_to_complete.py "$bed_file" "$gbff_file"
    else
        echo "Corresponding bed file not found for $gbff_file"
    fi
done



# Script format
## 1. convert gbff to txt file
## 2. run the rage_id script
## 3. subtract 1 from the start values
## 4. then try scripts to identify complete RAGEs