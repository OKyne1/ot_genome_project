#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J flye-overlap
#SBATCH --output=flye_overlap_%j.out
#SBATCH --error=flye_overlap_%j.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`
# This line has been added in but my not work properly - still hasn't been tested.
echo "Usage: $0 $@"

source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate assembly

# This code filters out reads below a specified length and then assembles the remaining reads using flye.

# Activate the correct environment (assembly)

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reads><overlap>"
    exit 1
fi

reads="$1"
overlap="$2"

flye --nano-hq "$reads" --threads 6 --min-overlap "$overlap" --out-dir assemblies_flye_overlap_"$overlap"

echo "flye is complete"
