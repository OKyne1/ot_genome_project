#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J fastq
#SBATCH --output=fastq_%j.out
#SBATCH --error=fastq_%j.err
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=12G
#SBATCH -p short

echo started=$(date)
echo "job=$SLURM_JOB_ID"
echo "hostname=$(hostname)"
echo "OS=$(uname -s)"
echo "username=$(whoami)"

# Activate the correct environment (assembly)

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fastq> <percentage>"
    exit 1
fi

fastq="$1"
percentage="$2"

cat "$fastq" | awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' |
awk 'BEGIN { srand(systime() + PROCINFO["pid"]); } NR % 4 == 0 { count++; } END { target = count * (percentage / 100); } { s = ++x <= target ? x : int(rand() * x); if (s <= target) { R[s] = $0; } } END { for (i in R) print R[i]; }' percentage="$percentage" |
awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > "subsampled_"$percentage".fastq"}'

echo "well an end has come..."
