#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J filter_host
#SBATCH --output=filter_host.out
#SBATCH --error=filter_host.err
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=$(date)
echo "job=$SLURM_JOB_ID"
echo "hostname=$(hostname)"
echo "OS=$(uname -s)"
echo "username=$(whoami)"

# Ensure that minimap2 and samtools are available in the environment (env: minimap2)
# The purpose of this code is to remove reads mapped to the Karp genome or unmapped and combine these into a single fastq file.

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reference.fasta> <ONT_reads.fastq>"
    exit 1
fi

reference="$1"
reads="$2"

# minimap2 mapping
minimap2 -x map-ont -t 6 -a -o aligned_reads.sam "$reference" "$reads"

# Convert SAM to BAM
samtools view -@ 4 -bS aligned_reads.sam > aligned_reads.bam

echo "alignment complete"

# Sort BAM file
samtools sort -@ 4 -o aligned_reads_sorted.bam aligned_reads.bam

echo "sorting complete"

# Index the sorted BAM file
samtools index aligned_reads_sorted.bam

echo "indexing complete"

# Extract unmapped reads
samtools view -f 4 -b aligned_reads_sorted.bam | samtools bam2fq - > unmapped.fastq

# This should allow us to view the chromosomes present.
#samtools view -H aligned_reads_sorted.bam # enables identification of the different chromosomes present. These are the things following SN:

# This code produces the Ot aligned reads. Change the chromosome name as appropriate
#samtools view -b -@ 4 aligned_reads_sorted.bam LS398548.1 > mapped_ot_reads.bam

# Filter out unmapped reads before extracting mapped reads
samtools view -F 4 -b -@ 4 aligned_reads_sorted.bam LS398548.1 > mapped_ot_reads.bam

# Convert to a fastq format
samtools bam2fq mapped_ot_reads.bam > mapped_ot_reads.fastq

# merge the fastq files
cat unmapped.fastq mapped_ot_reads.fastq > ot_and_unmapped.fastq

echo "mapping complete"
