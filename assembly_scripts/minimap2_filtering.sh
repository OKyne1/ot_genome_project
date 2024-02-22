#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J filter_host
#SBATCH --output=filter_host.out
#SBATCH --error=filter_host.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

# samtools often has conflicts with other packages
# Ensure the correct environment is activated for this task!

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reference.fasta> <ONT_reads.fastq>"
    exit 1
fi

reference="$1"
reads="$2"


minimap2 -x map-ont -t 6 -a -o data_out.sam "$reference" "$reads"
samtools sort -@ 4 -o data_out_sorted.bam data_out.sam
samtools index data_out_sorted.bam
samtools view -f 4 data_out_sorted.bam | samtools bam2fq - > unmapped.fastq

echo "mapping complete"
