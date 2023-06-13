#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --exclude=compute-0-11

# Load Anaconda environment
source /home/tsatler/anaconda3/etc/profile.d/conda.sh

input=$1
input_pdbs=($(cat "$input"))
pdb=${input_pdbs[$SLURM_ARRAY_TASK_ID]}

target_chain=$2
binder_chain=$3
out_file=$4

########################
# Binder optimization
########################
echo "Running binder analysis"

conda activate colabthread

python helper_scripts/binder_analysis.py $pdb $target_chain $binder_chain $out_file