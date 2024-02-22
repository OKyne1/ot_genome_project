#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J nanosim
#SBATCH --output=nanosim_%j.out
#SBATCH --error=nanosim_%j.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

#activate the nanosim environment

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reference_genome> <reads>"
    exit 1
fi

ref="$1"
reads="$2"

read_analysis.py genome -i "$reads" -rg "$ref" -a minimap2 -t 6

echo "well the things seems to have stopped..."
