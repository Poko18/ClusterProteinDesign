#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-2 #range for example 3HB_folder

### Runs ColabFold on multiple PDB files in a folder
# Note: Adjust the array range accordingly

# Arguments:
# 1) PDB file or Folder containing PDB files
# 2) Specify interaction analysis: "binder" or "binder-second"
#    (ignore if interaction analysis is not required)
# 3) Additional ColabFold arguments

# Example command: sbatch AF2_pdb_prediciton.sh examples/3HB_folder/ binder-second --num-recycle 16

# Activate the AlphaFold2 environment
. /home/aljubetic/bin/set_up_AF2.3.sh

# Prepare pdb file, fasta file and output folder
input=$1
if [ -d $input ]; then
    # Folder input
    input_files=($input/*.pdb)
    file=${input_files[$SLURM_ARRAY_TASK_ID]}
else
    # Single file input
    file=$input
fi

# Prepare pdb file, fasta file and output folder
base_name="$(basename "$file" | sed 's/\.[^.]*$//')"
current_path=$(pwd)
folder_base_name="$(basename "$input" | sed 's/\.[^.]*$//')"
out_dir="$current_path/output/$folder_base_name/tmp"
fasta_file=$out_dir/$base_name.fasta

# Create the output directory
if [ ! -d $out_dir ]
then
    mkdir -p $out_dir
fi

# Generate fasta file
if [ ! -f $fasta_file ]
then
    python helper_scripts/pdb2fasta.py $file $fasta_file
fi

# Check if interaction analysis is required
if [ "$2" == "binder" ]; then
    flags="--calculate_interaction"
elif [ "$2" == "binder-second" ]; then
    flags="--calculate_interaction --binder_sequence_second"
fi

shift 2
additional_args="$@"
echo analysis flags: $flags
echo colabfold args: $additional_args

# Run AF2
CMD="/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/python /home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch $additional_args $fasta_file $out_dir"
echo running AF2 command:
echo $CMD
# Run the prediction command
$CMD

# Copy the best-ranked PDB and JSON files to the output directory
final_pdb=$current_path/output/$folder_base_name/$base_name-AF2.pdb
final_json=$current_path/output/$folder_base_name/$base_name-AF2.json
cp $out_dir/$base_name*rank_001*.pdb $final_pdb
cp $out_dir/$base_name*rank_001*.json $final_json

### Analysis
echo running analysis
python helper_scripts/pdb_to_pdb_analysis.py \
        $file $final_pdb $final_json $current_path/output/$folder_base_name $flags
echo done!