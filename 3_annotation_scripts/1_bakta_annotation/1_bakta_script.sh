#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J bakta
#SBATCH --output=bakta_%j.out
#SBATCH --error=bakta_%j.err
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu=15G
#SBATCH -p short

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

######################################################### Bakta genome annotation #########################################################################
# Usage: $0 <prefix> <bakta_input_combined.faa> <genome>
# bakta_input_combined.faa is gene annotations from the Salje group 2023 genome paper and certain ncbi sequences
# Some of the sequences used above use modified confidence levels, these are; scaA, scaC, tsa22 and tsa56 due to their variability
# Additional sequences taken from ncbi were; cinA, scaA, scaB, scaD, scaE, secA, traA, traG, tsa22, tsa47 and tsa56 (excluding incomplete or partial genes)
###########################################################################################################################################################


# Activate bakta environment
source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate bakta

# Ensure the correct number of command line entries
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <prefix> <bakta_input_combined> <genome>"
    exit 1
fi

# Prefix for naming the files and also used in the locus tag
prefix="$1"
# Trusted proteins which we are testing/adding these are found in the same directory as the script and are called; bakta_input_combined.faa
input_proteins="$2"
# The genome to be annotated
genome="$3"
# Path to bakta database
database="/well/moru-batty/projects/ot/annotation/bakta_db/db"

# Bakta command line usage
bakta --db "$database" --prefix "$prefix" --complete --proteins "$input_proteins" --gram - --locus-tag "$prefix" --threads 8 --output "$prefix" "$genome" --keep-contig-headers

echo "Script finished"
exit 0