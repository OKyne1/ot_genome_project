#!/bin/bash

################################## hmmer search ##############################################
# This script is used to run faa files against hmmer files to identify tpr and ank domains
# The rage environment is used for this script
##############################################################################################

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <hmm_files...> <faa_files...>"  # Display usage message
    exit 1  # Exit with error code
fi

# Extract the HMM files
hmm_files=()
for file in "$@"; do
    if [[ "$file" == *.hmm ]]; then  # Check if file has .hmm extension
        hmm_files+=("$file")  # Add file to hmm_files array
    fi
done

# Extract the FAA files
faa_files=()
for file in "$@"; do
    if [[ "$file" == *.faa ]]; then  # Check if file has .faa extension
        faa_files+=("$file")  # Add file to faa_files array
    fi
done

# Loop through each HMM file
for hmm_file in "${hmm_files[@]}"; do
    # Extract the filename without extension
    hmm_name=$(basename "$hmm_file" .hmm)

    # Loop through each FAA file
    for faa_file in "${faa_files[@]}"; do
        # Extract the filename without extension
        faa_name=$(basename "$faa_file" .faa)

        # Run hmmsearch command
        # The E value given here is the one used in the bakta scripts and represents a fairly stringent threshold.
        hmmsearch -E 1E-10 --tblout "${faa_name}_${hmm_name}.txt" "$hmm_file" "$faa_file"
    done
done