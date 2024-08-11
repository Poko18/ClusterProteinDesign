#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

# Setup RF2 environment
source /home/tsatler/RFdif/RoseTTAFold2/set_up_RF2.sh

rf2_path="/home/tsatler/RFdif/RoseTTAFold2"

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_pdb> <output_folder>"
    exit 1
fi

# Store the arguments in variables
input_pdb="$1"
output_folder="$2"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Run the Python script to generate the MSA
python helper_scripts/get_msa_rf.py "$input_pdb" "$output_folder"
