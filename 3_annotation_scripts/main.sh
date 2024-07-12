#!/bin/bash

# rage environment needs to be activated

# This script needs to be run from a directory containing only the gbff files.
# The ank and tpr stuff is still a bit messy as it doesnt try to match files. This isn't an issue, it just adds ~3 min to the run time (with 8 gbff files).  

# Get the directory of the current script
SCRIPT_DIR=$(dirname "$0")

# Activate conda environment and run cog removal
#conda activate rage
python "$SCRIPT_DIR/2_cog_removal/scripts/1_cog_removal.py" *.gbff

# Activate conda environment and run domain annotation
conda activate hmmer
bash "$SCRIPT_DIR/3_domain_annotation/domain_annotation_package/scripts/domain_annotation.sh" modified*.gbff
# Clean up hmmer and faa files if necessary
# rm *.hmmer *.faa

# Reactivate rage environment and run Ankyrin ID
#conda activate rage
bash "$SCRIPT_DIR/4_ank_annotation/scripts/gene_completeness_master_script.sh" modified*.gbff

# Run completeness checks
bash "$SCRIPT_DIR/5_completeness_checks/scripts/gene_completeness_master_script.sh" modified*anked.gbff

# Create directory for complete RAGEs and move files
mkdir -p processing_outputs
mv modified_*_anked_completeness_checked.gbff processing_outputs/
#rm modified_*
#rm -r hmmer_outputs 
#rm -r faa_files
# Change directory and run RAGE classification
cd processing_outputs
bash "$SCRIPT_DIR/../6_rage_classification/scripts/main.sh" modified_*_anked_completeness_checked.gbff
#rm *.txt *.bed