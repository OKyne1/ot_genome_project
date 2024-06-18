#!/bin/bash

# Usage: python gbff_2_txt.py <genbank_file1> <genbank_file2> ...
# This should be run in an empty directory (+gbff files, tests had these in a directory called gbff)

script_dir2=$(dirname "$BASH_SOURCE")
echo "Script directory: $script_dir2"

## genbank to txt table
echo "Running 1_gbff_2_txt_complement.py..."
python3 "$script_dir2"/1_gbff_2_txt_complement.py "$@"

## Changing the locations
#mkdir gbff_tables
#mv *.txt gbff_tables

for file in *.txt; do
    echo "Removing first line of $file..."
    sed -i '1d' "$file"
done

## Define RAGE derived regions
echo "Running 2_rage_derived_regions.py..."
python3 "$script_dir2"/2_rage_derived_regions.py -tab *.txt -list "$script_dir2"/rage_proteins/lists_combined.txt -exclusion "$script_dir2"/rage_proteins/exclusion_list.txt
echo "RAGE derived regions done"

#mkdir bed_files
#mv *_rage_derived.bed bed_files 

# Define potential RAGE boundaries
echo "Running 2_rage_boundaries.py..."
python3 "$script_dir2"/2_rage_boundaries.py *.txt
echo "RAGE boundaries done"
#mv *_range_boundaries.bed bed_files

# RAGE boundary validation
echo "Starting RAGE boundary validation..."

# Iterate through the directory
for file1 in *_rage_boundaries.bed; do
    # Extract the base filename (without extension)
    base_name=$(basename "$file1" _rage_boundaries.bed)
    
    # Form the corresponding file2 name
    file2="${base_name}_rage_derived.bed"
    echo "File1: $file1"
    echo "File2: $file2"
    
    # Check if the corresponding file2 exists
    if [ -e "$file2" ]; then
        # Run bedtools intersect for the pair and save output to a file
        output_file="${base_name}_validated_boundaries.bed"
        echo "Running bedtools intersect..."
        bedtools intersect -a "$file1" -b "$file2" -f 1 -wa > "$output_file"
        echo "Output saved to $output_file"
    else
        echo "Corresponding file not found for $file1"
    fi
done
echo "Combined beds done"

# Directories containing the files
gbff_dir="gbff_tables"
#bed_dir="bed_files"

echo "Processing GBFF tables..."

# Iterate through the gbff_tables directory
for gbff_file in *.txt; do
    # Extract the base filename (without extension) from gbff_file
    base_name=$(basename "$gbff_file" .txt)
    
    # Form the corresponding filename in bed_files
    bed_file="${base_name}_validated_boundaries.bed"
    echo "Processing GBFF file: $gbff_file"
    echo "Expected BED file: $bed_file"

    # Check if the corresponding bed_file exists
    if [ -e "$bed_file" ]; then
        echo "Found BED file: $bed_file"
        # Run the Python script for the pair
        echo "Running: python3 \"$script_dir2/3_boundaries_to_complete.py\" \"$bed_file\" \"$gbff_file\""
        python3 "$script_dir2/3_boundaries_to_complete.py" "$bed_file" "$gbff_file"
        echo "Processed GBFF file: $gbff_file"
    else
        echo "Corresponding BED file not found for $gbff_file"
    fi
done

echo "Script execution completed."
mkdir rage_derived complete_rage

mv *_rage_derived.bed rage_derived
mv *_complete_RAGEs.bed complete_rage