#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J kraken2
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu=15G
#SBATCH -p short

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# Need to remove mouse reads first! otherwise mis-hits will occur. Should ideally also remove lambda phage.

# This code takes a single fastq.gz file and runs kraken2 on it, the specific db needs to be specified

source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate kraken

rootdir="/well/moru-batty/projects/ot"
#db="/well/moru-batty/projects/krakendb/standard-plus-mouse"
db="/well/moru-batty/projects/krakendb/standard"
sample="sample"

# Check if the input FASTQ file is provided as a command-line argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <fastqfile>"
    echo "Please specify the input FASTQ file."
    exit 1
fi

# Get the input FASTQ file from the command-line argument
fastqfile="$1"

# Extract the base filename (without path and extension)
input_filename=$(basename "$fastqfile" .fastq.gz)

# Set the sample name to the input filename
sample="$input_filename"

kraken2 --db ${db} --threads 6 --output ${sample}.out --report ${sample}.report --gzip-compressed ${fastqfile}

echo "finished="`date`
exit 0

