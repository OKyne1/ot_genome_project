#!/bin/bash
#SBATCH -A moru-batty.prj
#SBATCH -J trycycler
#SBATCH --output=trycycler_%A_%a.out
#SBATCH --error=trycycler_%A_%a.err
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=15G
#SBATCH --array 1-12:3
#SBATCH -p long

# 1 make sure the environemnt is assembly 
# Changed this to long at some point - may want to consider this (was killed when short at one point)
# For some reason gfa files 2,5,8,11 are not being removed from the location it was run
# Also 01 failed to be produced on the first run

# Define the number of threads
threads=4
# Could change this so that there are 12 threads (not exactly sure how...)

# Calculate the formatted_task_id with leading zeros
formatted_task_id=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

# Calculate the formatted_task_id_plus1 and formatted_task_id_plus2 without leading zeros
task_id_plus1=$((SLURM_ARRAY_TASK_ID + 1))
task_id_plus2=$((SLURM_ARRAY_TASK_ID + 2))

# adding leading zeros
formatted_task_id_plus1=$(printf "%02d" $task_id_plus1)
formatted_task_id_plus2=$(printf "%02d" $task_id_plus2)

# Run the desired commands using the formatted_task_id
flye --nano-hq read_subsets/sample_"$formatted_task_id".fastq --threads "$threads" --out-dir assembly_"$formatted_task_id" && cp assembly_"$formatted_task_id"/assembly.fasta assemblies/assembly_"$formatted_task_id".fasta && rm -r assembly_"$formatted_task_id"
/well/moru-batty/projects/ot/assembly/code_assembly/miniasm_and_minipolish.sh read_subsets/sample_"$formatted_task_id_plus1".fastq "$threads" > assembly_"$formatted_task_id_plus1".gfa && any2fasta assembly_"$formatted_task_id_plus1".gfa > assemblies/assembly_"$formatted_task_id_plus1".fasta && rm assembly_"$formatted_task_id_plus1".gfa
raven --threads "$threads" --disable-checkpoints read_subsets/sample_"$formatted_task_id_plus2".fastq > assemblies/assembly_"$formatted_task_id_plus2".fasta

echo "formatted_task_id: $formatted_task_id"
echo "formatted_task_id_plus1: $formatted_task_id_plus1"
echo "formatted_task_id_plus2: $formatted_task_id_plus2"

