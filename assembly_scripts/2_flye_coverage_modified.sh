#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J flye
#SBATCH --output=flye.out
#SBATCH --error=flye.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p short

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# 1 select length threshold
# 2 select percentage limit
# 3 ensure assembly environment is activated

# Check for the correct number of command-line arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <read_name>"
    exit 1
fi

read="$1"

flye --nano-hq "$read" --threads 6 --out-dir assemblies_flye --asm-coverage 50 --genome-size 2.3m

echo "flye is complete"
