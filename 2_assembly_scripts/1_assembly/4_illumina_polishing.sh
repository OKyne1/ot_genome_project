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
echo "Usage: $0 $@"

######################################### illumina polishing ############################################################
# We never really decided if illumina polishing was required - had too many other issues to work through.
# But here's a script to do it using polypolish
# Usage: $0 <ill_1><ill_2><assembly>
# We also never tested the impact of filtering out mouse reads from illumina reads first - may want to be considered
#########################################################################################################################

# Modify for your conda environment
source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate assembly

# Check for the correct number of command-line arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <ill_1><ill_2><assembly>"
    exit 1
fi

ill_1="$1"
ill_2="$2"
assembly="$3"

mkdir illumina_filtered

# fastp is used to identify the paired reads which we use for polishing
fastp --in1 "$ill_1" --in2 "$ill_2" --out1 illumina_filtered/ill_1.fastq.gz --out2 illumina_filtered/ill_2.fastq.gz --unpaired1 illumina_filtered/unpaired_1.fastq.gz --unpaired2 illumina_filtered/unpaired_$

# indexing assembly
bwa index "$assembly"

# alinging illumina reads to the assembly
bwa mem -t 6 -a "$assembly" illumina_filtered/ill_1.fastq.gz > align_1.sam
bwa mem -t 6 -a "$assembly" illumina_filtered/ill_2.fastq.gz > align_2.sam

# filtering only the higher quality read alignments
polypolish_insert_filter.py --in1 align_1.sam --in2 align_2.sam --out1 filtered_1.sam --out2 filtered_2.sam

# correcting the assembly based on the alignments. Change the output file type as desired
polypolish "$assembly" filtered_1.sam filtered_2.sam > assembly_polypolished.fasta

echo "Fini"