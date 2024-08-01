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
echo "Usage: $0 $@"

############################################# Script to assembly with a certain coverage level #########################################
# Usage: $0 <reads><coverage>
# Lines 2-15 are BMRC cluster specific
# Experiments found that >5000 bp reads generally resulted in the best assemblies; but coverage may want to be tested
########################################################################################################################################

# Conda activation will need to be modified for the user (by changing the source and environment name)
source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate assembly

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reads><coverage>"
    exit 1
fi

reads="$1"
cov="$2"

# this is the size of the karp genome, this should be fine for Ot, but may be modified for genomes of very different sizes
genome_size=2470000

# calculate base requirement for desired coverage
bases=$(( $2 * genome_size ))

mkdir filtered_reads 

# get the required reads for the desired coverage
filtlong --target_bases "$bases" "$reads" > ./filtered_reads/reads_coverage_"$cov".fastq
echo "filtlong filtering done"

# Assemble with the filtered reads
flye --nano-hq ./filtered_reads/reads_coverage_"$cov".fastq --threads 6 --out-dir assemblies_flye_coverage_"$cov"
echo "flye coverage code is complete"