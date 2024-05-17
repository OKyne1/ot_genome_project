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

# Activate bakta environment
source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate bakta

# May want to modify the code so the protein identified based on the location of the script?

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <prefix> <jeanne_proteins> <genome>"
    exit 1
fi

# Prefix for naming the files and also used in the locus tag
prefix="$1"
# Trusted proteins which we are testing/adding
jeanne_proteins="$2"
# The genome to be annotated
genome="$3"
# Path to bakta database
database="/well/moru-batty/projects/ot/annotation/bakta_db/db"

bakta --db "$database" --prefix "$prefix" --complete --proteins "$jeanne_proteins" --gram - --locus-tag "$prefix" --threads 8 --output "$prefix" "$genome" --keep-contig-headers

echo "Script finished"
exit 0
