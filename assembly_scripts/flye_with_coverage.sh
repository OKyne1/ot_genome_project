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

# Activate the correct environment (assembly)

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reads><coverage>"
    exit 1
fi

reads="$1"
cov="$2"
genome_size=2470000
bases=$(( $2 * genome_size ))

mkdir filtered_reads 

filtlong --target_bases "$bases" "$reads" > ./filtered_reads/reads_"$cov".fastq

echo "filtlong filtering seems to be done"

flye --nano-hq ./filtered_reads/reads_"$cov".fastq --threads 6 --out-dir assemblies_flye_"$cov"

echo "flye is complete"
