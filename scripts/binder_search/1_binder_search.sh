#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100
#SBATCH --exclude=compute-6-2,compute-6-0,compute-6-4,compute-6-1,compute-6-3

# Load Anaconda environment
source /home/tsatler/anaconda3/etc/profile.d/conda.sh

########################
# Binder search
########################

total_mpnn=1000
mpnn_per_design=100 # Filtered mpnn sequences for AF2
num_recycles=16
sampling_temp=0.15 # ProteinMPNN sampling temperature

prefix="cd5_lcb3_binder_v4"
target_pdb="../../targets/cd5/inputs/2ja4.pdb"
output="../../targets/cd5/outputs/fold_docking" # output gets created by rfdiffusion
#scaff_dir="../../scaffolds/rfdiff_filtered_scaffolds"
scaff_dir="../../testing/testing_scaffolds/lcb3"


########################
# ProteinMPNN and AF2
########################
echo Running AFDesign with MPNN sampling...

conda activate colabdesign

input_files=($rf_out*.pdb)
echo $input_files

# or run another bash script from here.. could be much faster.. but I dont know when it will end and would need to implement something for waiting

for ((i=0; i<${#input_files[@]}; i++)); do
  pdb_file=${input_files[$i]}
  echo $pdb_file AF2 design
  af_out=$output/mpnn_af2/${prefix}_${SLURM_ARRAY_TASK_ID}
  python ../helper_scripts/colabinder.py $pdb_file $af_out B A --sampling_temp $sampling_temp --num_recycles $num_recycles --num_seqs $total_mpnn --num_filt_seq $mpnn_per_design --results_dataframe $output --break_design --save_best_only
done

