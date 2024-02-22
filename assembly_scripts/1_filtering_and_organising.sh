#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J filter
#SBATCH --output=filter.out
#SBATCH --error=filter.err
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=12G
#SBATCH -p short

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# 1 select length threshold
# 2 select percentage limit
# 3 ensure assembly environment is activated

# Check for the correct number of command-line arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <read_name>"
    exit 1
fi

read="$1"

# Making new files for good organisation
mkdir reads reads_qc assemblies

# move the reads into reads
mv "$read" reads/"$read"

# Remove reads less than a specified length
## change this based on your data
filtlong --min_length 5000 reads/"$read" > reads_qc/ont_1k.fastq

# Remove the low quality reads
## adjust dependent on your data
filtlong --keep_percent 95 reads_qc/ont_1k.fastq > reads_qc/ont.fastq
rm reads_qc/ont_1k.fastq

echo "low quality and short reads have been removed"

# trycycler subsample --reads reads_qc/ont.fastq --out_dir read_subsets --genome_size 2.3m
# I commented this line as we aren't generating enough reads to subsample.
