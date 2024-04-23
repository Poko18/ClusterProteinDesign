#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2


# Load Anaconda environment
source /home/tsatler/anaconda3/etc/profile.d/conda.sh
conda activate colabthread

target_pdb=$1
target_chain=$2

scaf_dir=$3
input_pdbs=($scaf_dir/*.pdb)
sorted_input_pdbs=($(printf "%s\n" "${input_pdbs[@]}" | sort))
scaff_pdb=${sorted_input_pdbs[$SLURM_ARRAY_TASK_ID]}
echo $scaff_pdb

output=$4
iterations=$5

designs_per_iteration=$6
mpnn_per_design=$7
num_recycles=$8
sampling_temp=$9 # ProteinMPNN sampling temperature
mpnn_batch=$10


python helper_scripts/predict_one_binder_random_loc.py $target_pdb $target_chain $scaff_pdb \
        $output $iterations --designs_per_iteration $designs_per_iteration --proteinmpnn_per_input $mpnn_per_design \
        --num_recycles $num_recycles --sampling_temp $sampling_temp --mpnn_batch $mpnn_batch
