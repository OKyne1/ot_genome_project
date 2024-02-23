#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J flye-length
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

# This code filters out reads below a specified length and then assembles the remaining reads using flye.

# Activate the correct environment (assembly)

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reads><minlength>"
    exit 1
fi

reads="$1"
cov="$2"

mkdir filtered_reads 

filtlong --min_length "$cov" "$reads" > ./filtered_reads/reads_"$cov".fastq

echo "filtlong filtering seems to be done"

flye --nano-hq ./filtered_reads/reads_"$cov".fastq --threads 6 --out-dir assemblies_flye_"$cov"

echo "flye is complete"
