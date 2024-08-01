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
echo "Usage: $0 $@"

############################################# Script to extract reads >x and run flye ###############################################
# Usage: $0 <reads><coverage>
# Lines 2-15 are BMRC cluster specific
# Experiments found that >5000 bp reads generally resulted in the best assemblies; but it may be worth testing a few thresholds
#####################################################################################################################################

# Modify for your conda system
source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate assembly

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reads><minlength>"
    exit 1
fi

reads="$1"
minlength="$2"

mkdir filtered_reads 

# filter out reads <minlength value
filtlong --min_length "$minlength" "$reads" > ./filtered_reads/reads_minlength_"$minlength".fastq
echo "filtlong filtering done"

# flye assembly using default parameters and filtered reads
flye --nano-hq ./filtered_reads/reads_minlength_"$minlength".fastq --threads 6 --out-dir assemblies_flye_minlength_"$minlength"
echo "flye is complete"