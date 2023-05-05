#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1%5
#SBATCH --exclude=compute-6-2,compute-6-0,compute-6-4,compute-6-1,compute-6-3

# Load Anaconda environment
source /home/tsatler/anaconda3/etc/profile.d/conda.sh

########################
# RFdocking inputs
########################

# Define input parameters for RFdocking
num_of_diffusions=30 # Number of RF diffusions per script
total_mpnn=100
mpnn_per_design=10 # Filtered mpnn sequences for AF2
num_recycles=12
sampling_temp=0.2 # ProteinMPNN sampling temperature

prefix="cd5_helix_frag"
target_pdb="../../targets/cd5/inputs/2ja4.pdb"
output="../../targets/cd5/outputs/fold_docking"
hotspots='ppi.hotspot_res=[A288,A290,A298,A299,A328,A330,A362,A364]'
#scaff_dir="../../scaffolds/rfdiff_filtered_scaffolds"
#scaff_dir="../../testing/testing_scaffolds/lcb3"
scaff_dir="../../scaffolds/helix_seed_scaffolds"
########################
# RF diffusion
########################
echo "Doing $((num_of_diffusions * mpnn_per_design)) designs!"

# Generate ss, adj, and out_dir variables based on target_pdb
target_dir="$(dirname "$target_pdb")"
target_base="$(basename "${target_pdb%.*}")"
ss="${target_dir}/${target_base}_ss.pt"
adj="${target_dir}/${target_base}_adj.pt"
out_dir="${target_dir}/"

# Check if ss and adj files exist
if [ ! -f "$ss" ] || [ ! -f "$adj" ]; then
  echo "One or both of the SS and ADJ files do not exist. Running Python script to generate them..."
  
  conda activate pyro
  # Run Python script to generate ss and adj files
  python ../../scripts/helper_scripts/make_secstruc_adj.py --input_pdb "$target_pdb" --out_dir "$out_dir"
  
else
  echo "SS and ADJ files already exist."
fi

echo Running RFdiffusion for docking scaffolds to the target...

rf_out=$output/rf_dock/${prefix}_${SLURM_ARRAY_TASK_ID}/${prefix}_${SLURM_ARRAY_TASK_ID}

# Do something if folder/files already exist

conda activate SE3nv

python /home/tsatler/RFdif/RFdiffusion/scripts/run_inference.py \
scaffoldguided.target_path=$target_pdb \
inference.output_prefix=${rf_out} \
scaffoldguided.scaffoldguided=True \
$hotspots \
scaffoldguided.target_pdb=True \
scaffoldguided.target_ss=$ss \
scaffoldguided.target_adj=$adj \
scaffoldguided.scaffold_dir=$scaff_dir \
'potentials.guiding_potentials=["type:binder_distance_ReLU","type:binder_ROG","type:interface_ncontacts"]' \
potentials.guide_scale=2 \
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
if [ -d "$output/rf_dock/${prefix}_${SLURM_ARRAY_TASK_ID}/traj" ]; then
    rm -r "$output/rf_dock/${prefix}_${SLURM_ARRAY_TASK_ID}/traj" # delete trajectories
fi


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
  output="/home/tsatler/RFdif/ClusterProteinDesign/testing/helix_extraction"
  python ../helper_scripts/colabinder.py $pdb_file $af_out B A --sampling_temp $sampling_temp --num_recycles $num_recycles --num_seqs $total_mpnn --num_filt_seq $mpnn_per_design --results_dataframe $output --break_design --save_best_only
done

