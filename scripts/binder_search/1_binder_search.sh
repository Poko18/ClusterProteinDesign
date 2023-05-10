#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1-111%20
#SBATCH --exclude=compute-6-2,compute-6-0,compute-6-4,compute-6-1,compute-6-3

# Load Anaconda environment
source /home/tsatler/anaconda3/etc/profile.d/conda.sh

########################
# Binder search
########################

total_mpnn=2000
mpnn_per_design=2000 # Filtered mpnn sequences for AF2
num_recycles=16
sampling_temp=0.15 # ProteinMPNN sampling temperature

prefix="cd5_binders_v1"
input_dir="../../targets/cd5/outputs/fold_threading/filtered/"
output="../../targets/cd5/outputs/binders" # output gets created by rfdiffusion

########################
# ProteinMPNN and AF2
########################
echo Running AFDesign with MPNN sampling...

conda activate colabdesign

input_files=($input_dir*.pdb)
pdb_file=${input_files[$SLURM_ARRAY_TASK_ID]}

echo $pdb_file AF2 design
af_out=$output/mpnn_af2/${prefix}_${SLURM_ARRAY_TASK_ID}
#echo "python ../helper_scripts/colabinder.py $pdb_file $af_out A B --sampling_temp $sampling_temp --num_recycles $num_recycles --num_seqs $total_mpnn --num_filt_seq $mpnn_per_design --results_dataframe $output --save_best_only"
python ../helper_scripts/colabinder.py $pdb_file $af_out A B --sampling_temp $sampling_temp --num_recycles $num_recycles --num_seqs $total_mpnn --num_filt_seq $mpnn_per_design --results_dataframe $output --save_best_only


