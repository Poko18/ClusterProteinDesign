#!/bin/bash
#SBATCH -J mpnn
#SBATCH -p gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --nodelist=compute-0-12

source /home/tsatler/anaconda3/etc/profile.d/conda.sh
conda activate mlfold

# Check if input is provided
if [ -z "$1" ]; then
  echo "Error: Input file not provided. Input should be a folder with pdb files or a signle .pdb file."
  echo "Usage: script.sh <input_file>"
  exit 1
fi

###### FILE PATHS ######

input_file=$1
num_of_seq=1000
output_folder="$(pwd)/output/$(basename "$input_file" | sed 's/\.[^.]*$//')"

########################

if [ ! -d "$output_folder/tmp" ]
then
    mkdir -p "$output_folder/tmp"
fi

path_for_parsed_chains="$output_folder/tmp/parsed_pdbs.jsonl"
path_for_assigned_chains="$output_folder/tmp/assigned_pdbs.jsonl"
path_for_fixed_positions="$output_folder/tmp/fixed_pdbs.jsonl"
path_for_tied_positions="$output_folder/tmp/tied_pdbs.jsonl"

###### DESIGN PARAMETERS ######

chains_to_design="A"

###### HELPER SCRIPTS ######

# Check if input is a folder
if [ -d "$input_file" ]; then
  python /home/tsatler/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path="$input_file" --output_path="$path_for_parsed_chains"
# Check if input is a .pdb file
elif [ -f "$input_file" ]; then
  python /home/tsatler/ProteinMPNN/helper_scripts/parse_one_pdb.py --input_path="$input_file" --output_path="$path_for_parsed_chains"
else
  echo "Error: Invalid input. Input should be a folder with pdb files or a signle .pdb file."
  exit 1
fi

python /home/tsatler/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

#python /home/tsatler/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

#python /home/tsatler/ProteinMPNN/helper_scripts/make_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --chain_list "$chains_to_design" --position_list "$tied_positions"

###### RUNNING PARAMETERS ######

python /home/tsatler/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --tied_positions_jsonl $path_for_tied_positions \
        --out_folder $output_folder \
        --num_seq_per_target $num_of_seq \
        --sampling_temp "0.20" \
        --seed 37 \
        --batch_size 1
