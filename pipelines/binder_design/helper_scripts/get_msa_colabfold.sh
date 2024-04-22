#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Activate the AlphaFold2 environment
. /home/aljubetic/bin/set_up_AF2.3.sh

# Store the arguments in variables
input_pdb="$1"
output_folder="$2"
additional_af2_args="--num_recycles 1"

fasta_file=$output_folder/fasta.fasta
# Generate fasta file
if [ ! -f $fasta_file ]
then
    python helper_scripts/pdb2fasta.py $input_pdb $fasta_file
fi

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Run ColabFold
CMD="/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/python /home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch $additional_args $fasta_file $output_folder"
$CMD
