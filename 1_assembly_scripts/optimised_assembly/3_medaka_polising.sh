#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J polish
#SBATCH --output=polish_%j.out
#SBATCH --error=polish_%j.err
#SBATCH --cpus-per-task 2
#SBATCH --mem-per-cpu=15G
#SBATCH -p long

echo started=`date`
echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`

# Used for polishing assemblies with nanopore reads
## 1. Activate assembly environment
## 2. Modify the -m for the correct flow cell model

echo "script is running"
# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <assembly> <reads>"
    exit 1
fi

assembly="$1"
reads="$2"
# usage suggests not to thread above 2 as it scales inefficiently past this, but says that it may use an additional 4 threads for the reading and preparing of data (hence 6 threads)
medaka_consensus -t 2 -i "$reads" -d "$assembly" -o medaka_polished -m r1041_e82_400bps_hac_v4.2.0

echo "medaka polishing code is finished"

