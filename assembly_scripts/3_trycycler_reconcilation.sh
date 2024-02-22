#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J cluster_reconciliation
#SBATCH --output=cluster_rec.out
#SBATCH --error=cluster_rec.err
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=15G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# Check for the correct number of command-line arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <assembly>"
    exit 1
fi

assembly="$1"


trycycler reconcile --reads reads_qc/ont.fastq --cluster_dir "$assembly"
