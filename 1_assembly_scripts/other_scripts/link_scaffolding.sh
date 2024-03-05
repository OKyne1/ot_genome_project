#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J link
#SBATCH --output=link_%j.out
#SBATCH --error=link_%j.err
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=$(date)
echo "job=$SLURM_JOB_ID"
echo "hostname=$(hostname)"
echo "OS=$(uname -s)"
echo "username=$(whoami)"
echo "Usage: $0 $@"

source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
# This environment is just links and its dependencies
conda activate links

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <draft_assembly> <reads>"
    exit 1
fi

draft="$1"
reads="$2"

mkdir scaffolded

LINKS -f "$draft" -s "$reads" -o scaffolded_"$draft" -t 4

echo "well something may have happened..."
