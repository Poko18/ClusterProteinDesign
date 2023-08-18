#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

set -e
# Load Anaconda environment
source /home/tsatler/anaconda3/etc/profile.d/conda.sh


########################
# RFdocking inputs
########################

# Define input parameters for RFdocking
num_of_diffusions=$1 # Number of RF diffusions per script
total_mpnn=$2
mpnn_per_design=$3 # Filtered mpnn sequences for AF2
num_recycles=$4
sampling_temp=$5 # ProteinMPNN sampling temperature


prefix=$6
target_pdb=$7
output="output/$prefix" # output gets created by rfdiffusion
hotspots=$8
echo $hotspots
contigs=$9
echo "$contigs"
thread_sequence=${10}
echo "$thread_sequence"


########################
# RF diffusion
########################
echo "Doing $((num_of_diffusions * mpnn_per_design)) designs!"

# Generate ss, adj, and out_dir variables based on target_pdb
target_dir="$(dirname "$target_pdb")"
target_base="$(basename "${target_pdb%.*}")"
out_dir="${target_dir}/"
subfolder_name="random_scaff"


echo Running RFdiffusion for docking scaffolds to the target...

rf_out=$output/rf_dock/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}

# Do something if folder/files already exist

conda activate SE3nv

python /home/tsatler/RFdif/RFdiffusion/scripts/run_inference.py \
inference.input_pdb=/home/tsatler/RFdif/pcsk_binders/pcsk9_target_short.pdb \
inference.ckpt_override_path=/home/tsatler/RFdif/RFdiffusion/models/Complex_beta_ckpt.pt \
"$contigs" \
inference.output_prefix=${rf_out} \
$hotspots \
'potentials.guiding_potentials=["type:binder_ROG,weight:3","type:interface_ncontacts","type:binder_distance_ReLU"]' \
potentials.guide_scale=3 \
potentials.guide_decay="quadratic" \
inference.num_designs=$num_of_diffusions \
denoiser.noise_scale_ca=0 \
denoiser.noise_scale_frame=0


# Remove some files
current_dir=$(pwd)
if [ -d "$current_dir/outputs" ]; then
    rm -r "$current_dir/outputs" # remove hydra stuff
fi
if [ -e "$rf_out*.trb" ]; then
    rm "$rf_out*.trb" # delete trb files
fi
if [ -d "$output/rf_dock/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}/traj" ]; then
    rm -r "$output/rf_dock/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}/traj" # delete trajectories
fi


########################
# ProteinMPNN and AF2
########################
echo "Running AFDesign with MPNN sampling... (and threading)"

conda activate colabthread

input_files=($rf_out*.pdb)
echo $input_files

# or run another bash script from here.. could be much faster.. but I dont know when it will end and would need to implement something for waiting

if [ -f "colabinder.py" ]; then
    script="colabinder.py"
else
    script="helper_scripts/colabinder.py"
fi

for ((i=0; i<${#input_files[@]}; i++)); do
  pdb_file=${input_files[$i]}
  echo $pdb_file AF2 design - random scaffolds
  af_out=$output/mpnn_af2/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}
  python $script $pdb_file $af_out B A --sampling_temp $sampling_temp \
  --num_recycles $num_recycles --num_seqs $total_mpnn --num_filt_seq $mpnn_per_design \
  --results_dataframe $output --break_design --save_best_only $thread_sequence
done