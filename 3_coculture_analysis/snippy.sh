#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J snippy
#SBATCH --output=snippy_%j.out
#SBATCH --error=snippy_%j.err
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=8G
#SBATCH -p short

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate snippy

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome><reads><output>"
    echo "Please specify the input FASTQ file."
    exit 1
fi

genome="$1"
reads="$2"
output="$3"

singularity run snippy.sif snippy --cpus 4 --ref genome --se reads --outdir output

# Need to add in the singularity run stuff to this command
# singularity run snippy.sif snippy ..........
