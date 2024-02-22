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

# 1 ensure assembly environment is activated
# 2 Decide on the threshold limits

# Note that this is not the same as the filtering done by filtlong as the filtered reads will be used in later assembly stages and in O. tsutsugamushi can result in increased contig number.

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

