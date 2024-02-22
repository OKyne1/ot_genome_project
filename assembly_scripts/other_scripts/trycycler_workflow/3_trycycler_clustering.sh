#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J cluster
#SBATCH --output=cluster.out
#SBATCH --error=cluster.err
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=15G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

#--------------------------------------------------------------
# New job
# Remove unnecessary data
rm -r read_subsets

trycycler cluster --assemblies assemblies/*.fasta --reads reads_qc/ont.fastq --out_dir trycycler

# At this point it is important to decide what clusters are important to keep and which to discard
## Here is an example of what to do with bad clusters: mv trycycler/cluster_003 trycycler/bad_cluster_003
echo "trycycler clustering complete"
