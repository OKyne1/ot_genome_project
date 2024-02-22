#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J flye
#SBATCH --output=flye_%j.out
#SBATCH --error=flye_%j.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# 1 select length threshold
# 2 select percentage limit
# 3 ensure assembly environment is activated

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <read_name><coverage>"
    exit 1
fi

read="$1"
coverage="$2"

mkdir assemblies_flye_"$coverage"

flye --nano-hq "$read" --threads 6 --out-dir assemblies_flye_"$coverage" --asm-coverage "$coverage" --genome-size 2.47m

echo "flye is complete"

