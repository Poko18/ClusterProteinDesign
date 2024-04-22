#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Activate the AlphaFold2 environment
. /home/aljubetic/bin/set_up_AF2.3.sh

# Arguments
msa_folder=$1
filtered_binders_path=$2
binder_analysis=$3
prediction_tools="${@:4}"  # Store all arguments from the 4th position onward as prediction tools

# Prepare a3m folder
if [[ $prediction_tools == *"colabfold"* ]]; then
    # ColabFold mode
    msa_folder_fullpath="$msa_folder/msa_inputs_af2"  # Point to the af2 folder
elif [[ $prediction_tools == *"rosettafold2"* ]]; then
    # RoseTTAFold2 mode
    msa_folder_fullpath="$msa_folder/msa_inputs_rf2"  # Point to the rf2 folder
else
    # Handle the case if neither tool is specified
    echo "No prediction tools specified or unsupported tools."
    exit 1
fi

# Prepare a3m file
if [ -d $msa_folder_fullpath ]; then
    # Folder input
    input_files=($msa_folder_fullpath/*.a3m)
    file=${input_files[$SLURM_ARRAY_TASK_ID]}
else
    # Single file input
    file=$msa_folder
fi

# Get basename from file
basename=$(basename "$file" .a3m)
echo $basename

# A3M files
file_af2=$msa_folder/msa_inputs_af2/$basename.a3m
file_rf2=$msa_folder/msa_inputs_rf2/$basename.a3m
echo $file_af2

# Create the output directory
scores_path=$filtered_binders_path/scores
if [ ! -d $scores_path ]
then
    mkdir -p $scores_path
fi

# Check if interaction analysis is required
if [ "$binder_analysis" == "binder" ]; then
    flags="--calculate_interaction"
elif [ "$binder_analysis" == "binder-second" ]; then
    flags="--calculate_interaction --binder_sequence_second"
fi
echo binder analysis - $binder_analysis - $flags

# Prediction tools to use
echo using prediction tools: $prediction_tools


if [[ $prediction_tools == *"colabfold"* ]]; then

    # Create AF2 out dir if it doesnt exist
    af2_out_dir=$filtered_binders_path/af2
    if [ ! -d $af2_out_dir ]
    then
        mkdir -p $af2_out_dir
    fi

    additional_af2_args="--num_recycles 12"

    # Run ColabFold
    CMD="/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/python /home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch $additional_args $file_af2 $af2_out_dir"
    echo running AF2 command:
    echo $CMD
    # Run the prediction command
    $CMD
    echo AF2 prediction done!

    # Prepare for analysis
    ref_pdb=$filtered_binders_path/$basename.pdb
    final_pdb=$filtered_binders_path/$basename-AF2.pdb
    final_json=$filtered_binders_path/$basename-AF2.json

    # Copy files analysis
    cp $af2_out_dir/$basename*rank_001*.pdb $final_pdb
    cp $af2_out_dir/$basename*rank_001*.json $final_json

    # Run analysis
    echo running analysis
    python helper_scripts/analyze_pdbs.py \
            $ref_pdb $final_pdb af2 $scores_path --out_name scores_af2 --sample_json $final_json $flags
    echo AF2 analysis done!

fi

if [[ $prediction_tools == *"rosettafold2"* ]]; then

    #filtered_binders_path_rf2=$filtered_binders_path/msa_inputs_rf2
    # Create AF2 out dir if it doesnt exist
    rf2_out_dir=$filtered_binders_path/rf2
    if [ ! -d $rf2_out_dir ]
    then
        mkdir -p $rf2_out_dir
    fi

    additional_rf2_args="--num_recycles 12"

    # Run Rosettafold2
    source /home/tsatler/RFdif/RoseTTAFold2/set_up_RF2.sh
    rf2_path="/home/tsatler/RFdif/RoseTTAFold2"

    # Set model_params if not provided
    if [ -z "$model_params" ]; then
        model_params="$rf2_path/network/weights/RF2_apr23.pt"
    fi

    # Print the python command
    cmd="python $rf2_path/network/predict_mmseq.py $file_rf2 $rf2_out_dir $additional_rf2_args --model_params $model_params"
    echo running RF2 command:
    echo $cmd

    # Run the Python script
    $cmd
    echo RF2 prediction is done

    # Prepare for analysis
    ref_pdb=$filtered_binders_path/$basename.pdb
    final_pdb=$filtered_binders_path/$basename-RF2.pdb
    final_json=$filtered_binders_path/$basename-RF2.json

    # Copy files analysis
    cp $rf2_out_dir/$basename*.pdb $final_pdb
    cp $rf2_out_dir/$basename*.json $final_json

    # Run analysis

    source /home/tsatler/anaconda3/etc/profile.d/conda.sh
    conda activate pyro

    echo running analysis
    python helper_scripts/analyze_pdbs.py \
            $ref_pdb $final_pdb rf2 $scores_path --out_name scores_rf2 --sample_json $final_json $flags
    echo RF2 analysis done!
fi

echo done!
