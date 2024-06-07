#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J kraken2
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=8G
#SBATCH -p short

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate qc

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome><reads><output>"
    echo "Please specify the input FASTQ file."
    exit 1
fi

reads="$2"
genome="$1"
output="$3"

snippy --cpus 4 --reference genome --se reads --output output