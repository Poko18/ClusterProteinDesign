#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

set -e
# Load Anaconda environment
#source /home/tsatler/anaconda3/etc/profile.d/conda.sh

########################
# RFdocking inputs
########################

# Define input parameters for RFdocking
num_of_diffusions=$1 # Number of RF diffusions per script
num_seqs=$2
num_filtered_seqs=$3 # Filtered mpnn sequences for AF2
num_recycles=$4
sampling_temp=$5 # ProteinMPNN sampling temperature

scaff_string=$6 #python list as a string
IFS=',' read -r -a scaff_list <<< "$scaff_string" #convert string into a list
echo $scaff_string
prefix=$7
target_pdb=$8
output="output/$prefix" # output gets created by rfdiffusion
hotspots=$9
for scaff_dir in "${scaff_list[@]}"
do
  subfolder_name=$(basename "$scaff_dir")
  echo $subfolder_name


  ########################
  # RF diffusion
  ########################
  echo "Doing $((num_of_diffusions * num_seqs)) designs!"

  # Generate ss, adj, and out_dir variables based on target_pdb
  target_dir="$(dirname "$target_pdb")"
  target_base="$(basename "${target_pdb%.*}")"
  ss="${target_dir}/${target_base}_ss.pt"
  adj="${target_dir}/${target_base}_adj.pt"
  out_dir="${target_dir}/"
  # Check if ss and adj files exist
  if [ ! -f "$ss" ] || [ ! -f "$adj" ]; then
    echo "One or both of the SS and ADJ files do not exist. Running Python script to generate them..."
    
    source /home/lhafner/anaconda3/bin/activate pyrosetta
    # Run Python script to generate ss and adj files
    python helper_scripts/make_secstruc_adj.py --input_pdb "$target_pdb" --out_dir "$out_dir"
    
  else
    echo "SS and ADJ files already exist."
  fi

  echo Running RFdiffusion for docking scaffolds to the target...

  rf_out=$output/rf_dock/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}

  # Do something if folder/files already exist
  source /home/lhafner/anaconda3/bin/activate SE3nv-auto

  python /home/lhafner/RFdiffusion/scripts/run_inference.py \
  diffuser.T=50 \
  scaffoldguided.target_path=$target_pdb \
  inference.output_prefix=${rf_out} \
  scaffoldguided.scaffoldguided=True \
  $hotspots \
  scaffoldguided.target_pdb=True \
  scaffoldguided.target_ss=$ss \
  scaffoldguided.target_adj=$adj \
  scaffoldguided.scaffold_dir=$scaff_dir \
  scaffoldguided.mask_loops=False \
  'potentials.guiding_potentials=["type:binder_ROG,weight:3"]' \
  potentials.guide_scale=2 \
  potentials.guide_decay="quadratic" \
  inference.num_designs=$num_of_diffusions \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0

  #"type:binder_distance_ReLU",

  # Remove some files
  current_dir=$(pwd)
  #if [ -d "$current_dir/outputs" ]; then
  #    rm -r "$current_dir/outputs" # remove hydra stuff
  #fi
  for file in "$output/rf_dock/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}/"*.trb; do
    if [ -f "$file" ]; then
      rm "$file"
    fi
  done
  if [ -d "$output/rf_dock/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}/traj" ]; then
      rm -r "$output/rf_dock/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}/traj" # delete trajectories
  fi

  ########################
  # ProteinMPNN and AF2
  ########################
  echo "Running AFDesign with MPNN sampling... (and threading)"

  source /home/lhafner/anaconda3/bin/activate colabthread-tsatler

  input_files=(${rf_out}*.pdb)
  echo $input_files

  # or run another bash script from here.. could be much faster.. but I dont know when it will end and would need to implement something for waiting

  if [ -f "colabinder_v2.py" ]; then
      script="colabinder.py"
  else
      script="helper_scripts/colabinder.py"
  fi

  #for ((i=0; i<${#input_files[@]}; i++)); do
  #  pdb_file=${input_files[$i]}
  #  echo $pdb_file AF2 design
  #  af_out=$output/mpnn_af2/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}
  #  python $script $pdb_file $af_out B A --sampling_temp $sampling_temp \
  #  --num_recycles $num_recycles --num_seqs $total_mpnn --num_filt_seq $mpnn_per_design \
  #  --results_dataframe $output --break_design --thread_sequence --save_best_only \
  #  --scaffold_folder $scaff_dir
  #done

  for ((i=0; i<${#input_files[@]}; i++)); do
    pdb_file=${input_files[$i]}
    echo $pdb_file AF2 design
    af_out=$output/mpnn_af2/${prefix}_${subfolder_name}_${SLURM_ARRAY_TASK_ID}
    python $script $pdb_file $af_out B A --sampling_temp $sampling_temp \
    --num_recycles $num_recycles --num_seqs $num_seqs --num_filt_seq $num_filtered_seqs \
    --results_dataframe $output --save_best_only
  done
done
