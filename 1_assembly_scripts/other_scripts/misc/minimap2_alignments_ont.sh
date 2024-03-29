#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J minimap
#SBATCH --output=minimap_%j.out
#SBATCH --error=minimap_%j.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=12G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# This code is for the visualisation (e.g. with IGV) of reads mapped to a genome. It was used to see whether there was an uneven distribution of reads across the O.tsutsugamushi genomes.

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <genome> <reads>"
    exit 1
fi

genome="$1"
reads="$2"

# Run minimap2 to align reads to the reference genome
minimap2 -x map-ont -t 6 -a -o mapped.sam "$genome" "$reads"

# Sort the SAM file and convert it to BAM
samtools sort -@ 4 -o sorted_mapped.bam mapped.sam

# Index the sorted BAM file
samtools index sorted_mapped.bam

# Extract primary mapped reads
samtools view -F 260 -q 1 -b sorted_mapped.bam > primary_mapped.bam

# Index primary mapped reads
samtools index primary_mapped.bam

echo "Well, it seems to have finished running..."
