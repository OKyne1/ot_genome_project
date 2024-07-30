#!/bin/bash

############################################################################################################################################################
# This script needs to be run from a directory containing only the gbff files within the 3_annotation_scripts directory (dropbox fucked things up)
#
# The ank and tpr stuff is still a bit messy as it doesn't try to match files. This isn't an issue, it just adds ~1 min to the run time (with 8 gbff files).
# This also results in some error messages in the annotation.out file (file xxxx not found) this generatlly isn't an issue in this specific case.
# **This script does not work if files have names with "-" in** due to the way it handles contigs.
#
############################################################################################################################################################

# Get the directory of the current script
SCRIPT_DIR=$(dirname "$(realpath "$0")")

echo "$SCRIPT_DIR"

# run cog removal
echo "################################################################ COG removal messages ########################################################" >> annotation.out 2>&1
python "$SCRIPT_DIR/2_cog_removal/1_cog_removal.py" *.gbff >> annotation.out

echo "script 1 of 6 complete (cog removal)"

# Create directory for input files
mkdir -p input_gbff processing_outputs

# Move the original input gbff files to input_gbff directory
for file in *.gbff; do
    if [[ ! $file == modified_* ]]; then
        mv "$file" input_gbff/
    fi
done

# run domain annotation
echo "######################################################## domain annotation ########################################################" >> annotation.out 2>&1

bash "$SCRIPT_DIR/3_domain_annotation/scripts/domain_annotation.sh" modified*.gbff >> annotation.out 2>&1

echo "script 2 of 6 complete (domain_annotation.sh)"

echo "#################################################### ank protein annotation ########################################################" >> annotation.out 2>&1

bash "$SCRIPT_DIR/4_ank_annotation/ank_completeness_master_script.sh" modified*.gbff >> annotation.out 2>&1

echo "script 3 of 6 complete (ank annotation script)"

# Remove files matching modified_*.gbff but not modified_*_anked.gbff
found=false
for file in modified_*.gbff; do
    if [[ $file == modified_*_anked.gbff ]]; then
        found=true
    else
        rm "$file"
    fi
done

if ! $found; then
    echo "No files matching the pattern modified_*_anked.gbff were found. This means the script failed at or before the 4_ank_annotation/scripts/gene_completeness_master_script.sh script."
fi

# Run completeness checks
echo "################################################ RAGE protein completeness checks #################################################" >> annotation.out 2>&1

bash "$SCRIPT_DIR/5_completeness_checks/gene_completeness_master_script.sh" modified*anked.gbff >> annotation.out 2>&1

echo "script 4 of 6 complete (rage protein completeness checks)"

rm *.txt

for file in *.*; do
    if [[ ! $file == modified_*_anked_completeness_checked.gbff && $file != annotation.out ]]; then
        rm "$file"
    fi
done

echo "######################################################## spoT annotation ########################################################" >> annotation.out 2>&1

bash "$SCRIPT_DIR/6_spoT/spot_script.sh" modified_*_anked_completeness_checked.gbff >> annotation.out 2>&1

echo "script 5 of 6 complete (spot_script.sh)"

# Create directory for complete RAGEs and move files
mkdir -p processing_outputs

# Check if files exist before moving
if ls modified_*_anked_completeness_checked_spotted.gbff 1> /dev/null 2>&1; then
    mv modified_*_anked_completeness_checked_spotted.gbff processing_outputs/
    rm *_checked.faa *_checked.gbff *_complete.txt *_checked.txt
else
    echo "No files matching modified_*_anked_completeness_checked_spotted.gbff found."
    exit 1
fi

# Check if files were moved successfully
echo "######################################################## RAGE classification ####################################################" >> annotation.out 2>&1

if ls processing_outputs/modified_*_anked_completeness_checked_spotted.gbff 1> /dev/null 2>&1; then
    # Remove unnecessary files and directories
    # rm modified_*
    rm -r hmmer_outputs 
    rm -r faa_files

    # Change directory
    cd processing_outputs || { echo "Failed to change directory to processing_outputs."; exit 1; }

    # Ensure SCRIPT_DIR is set
    if [ -z "$SCRIPT_DIR" ]; then
        echo "SCRIPT_DIR is not set." 
        exit 1
    fi

    # Check if main.sh is executable
    if [ ! -x "$SCRIPT_DIR/7_rage_classification/main.sh" ]; then
        echo "main.sh script is not executable or not found."
        echo "$SCRIPT_DIR/7_rage_classification/main.sh"
        exit 1
    fi

    # Run RAGE classification
    for file in modified_*_anked_completeness_checked_spotted.gbff; do
        bash "$SCRIPT_DIR/7_rage_classification/main.sh" $file >> ../annotation.out 2>&1
        rm *.txt *.bed
    done
    echo "script 6 of 6 complete (7_rage_classification/main.sh)"
else
    echo "Failed to move files to processing_outputs."
    exit 1
fi