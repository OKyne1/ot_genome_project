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
echo "Usage: $0 $@"

######################################################## Medaka polishing ###############################################################################
# Usage: $0 <assembly> <reads>
# This should use all orientia reads
# Ensure the flow cell technology is correct for the flow cell used
#########################################################################################################################################################

# Change for your conda environemnt/location
source /well/moru-batty/users/vhs789/miniforge/etc/profile.d/conda.sh
conda activate assembly

echo "script is running"
# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <assembly> <reads>"
    exit 1
fi

assembly="$1"
reads="$2"

# medaka command line argument, modify the flow cell technology as required
# Change the model as appropriate
medaka_consensus -t 2 -i "$reads" -d "$assembly" -o medaka_polished -m r1041_e82_400bps_hac_v4.2.0
echo "medaka polishing finished"