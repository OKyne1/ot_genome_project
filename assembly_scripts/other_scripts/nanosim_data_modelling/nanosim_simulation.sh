#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J nanosim2
#SBATCH --output=nanosim2_%j.out
#SBATCH --error=nanosim2_%j.err
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
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference_genome><model><number_of_reads>"
    exit 1
fi

ref="$1"
model="$2"
no="$3"

simulator.py genome -rg "$ref" -c "$model" -n "$no" -dna_type circular --fastq -t 6 -b guppy -s 0.5 -c /well/moru-batty/projects/ot/assembly/240116_karp_super/240220_nanosim/models/all_reads_model/training

echo "I guess this is over then. Probably worth checking the results."
