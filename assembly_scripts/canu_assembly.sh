#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J canu
#SBATCH --output=canu_%j.out
#SBATCH --error=canu_%j.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`


# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <read_name><coverage>"
    exit 1
fi

read="$1"
coverage="$2"

mkdir assemblies_canu_"$coverage"

canu gridEngine=slurm -p karp -d assemblies_canu_"$coverage" genomeSize=2.47m -trimmed -nanopore "$1" maxInputCoverage="$2"

echo "Canu is complete or atleast the code is finished"
