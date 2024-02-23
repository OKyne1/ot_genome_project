#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J polypolish
#SBATCH --output=illumina_%j.out
#SBATCH --error=illumina_%j.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# This code is for polishing the assembly with illumina reads. We may not use it in our workflow but here it is.
## 1. Activate assembly environment
## 2. Polish!

# Check for the correct number of command-line arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <ill_1><ill_2><assembly>"
    exit 1
fi

ill_1="$1"
ill_2="$2"
assembly="$3"

mkdir illumina_filtered

fastp --in1 "$ill_1" --in2 "$ill_2" --out1 illumina_filtered/ill_1.fastq.gz --out2 illumina_filtered/ill_2.fastq.gz --unpaired1 illumina_filtered/unpaired_1.fastq.gz --unpaired2 illumina_filtered/unpaired_$

bwa index "$assembly"

bwa mem -t 6 -a "$assembly" illumina_filtered/ill_1.fastq.gz > align_1.sam
bwa mem -t 6 -a "$assembly" illumina_filtered/ill_2.fastq.gz > align_2.sam

polypolish_insert_filter.py --in1 align_1.sam --in2 align_2.sam --out1 filtered_1.sam --out2 filtered_2.sam

# Choose the appropriate output format
polypolish "$assembly" filtered_1.sam filtered_2.sam > assembly_polypolished.fastq


echo "Well the script seems to have finished..."


