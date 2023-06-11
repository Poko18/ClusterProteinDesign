#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --exclude=compute-0-11

# Load Anaconda environment
source /home/tsatler/anaconda3/etc/profile.d/conda.sh

target_chain="A"
binder_chain="B"

input=$1
input_pdbs=($input/*.pdb)
pdb=${input_pdbs[$SLURM_ARRAY_TASK_ID]}

output=$2
iterations=$3
initial_proteinmpnns=$4

designs_per_iteration=$5
mpnn_per_design=$6
sort=$7

num_recycles=$8
sampling_temp=$9 # ProteinMPNN sampling temperature

########################
# Binder optimization
########################
echo "Running binder optimization"

conda activate colabthread

python helper_scripts/binder_opt.py $pdb $output \
    $target_chain $binder_chain $iterations --initial_proteinmpnns $initial_proteinmpnns --designs_per_iteration $designs_per_iteration\
    --proteinmpnn_per_input $mpnn_per_design --sort_by $sort \
    --num_recycles $num_recycles --sampling_temp $sampling_temp
