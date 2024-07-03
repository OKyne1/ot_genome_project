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
echo "Usage: $0 $@"

#################################################### Mouse removal code #################################################################################################
# This code filters out reads that map to the mouse genome and not the karp genome. Future use of this code may want to change the reference doc to include more Ot diversity and also get reads mapped to all of these different genomes.
# Lines 2-15 are specific to the slurm computing cluster (BMRC cluster). If running locally, remove these, threading may need to be changed in some of the lines.
# The environment activation stuff will need to be modified for the user, minimap2 and samtools should be present in the environment.
# Genomes used for this competitive mapping is found in ./genomes_for_mapping/ relative to this script and contains karp, mouse genome and the dilute control sequence (positive control)
#########################################################################################################################################################################

# conda environment activation
source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate minimap


# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reference.fasta> <ONT_reads.fastq>"
    exit 1
fi

reference="$1"
reads="$2"

# minimap2 mapping of reads against the combined reference
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

# This code produces the Ot aligned reads (based on the name of the karp chromosome: LS398548.1).  
# If the file/genomes are modified then this may need to be changed.
samtools view -b -@ 4 aligned_reads_sorted.bam LS398548.1 > mapped_ot_reads.bam

# Filter out unmapped reads before extracting LS398548.1 mapped reads
samtools view -F 4 -b -@ 4 aligned_reads_sorted.bam LS398548.1 > mapped_ot_reads.bam

# Convert to a fastq format
samtools bam2fq mapped_ot_reads.bam > mapped_ot_reads.fastq

# merge the fastq files
cat unmapped.fastq mapped_ot_reads.fastq > ot_and_unmapped.fastq
echo "mapping complete"
