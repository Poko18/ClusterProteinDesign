#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Activate the AlphaFold2 environment
#. /home/aljubetic/bin/set_up_AF2.sh
. /home/aljubetic/bin/set_up_AF2.3.sh

input_file=$1

current_path=$(pwd)
base_name="$(basename "$input_file" | sed 's/\.[^.]*$//')"
out_dir="$current_path/output/$base_name"

# Create the output directory
if [ ! -d $out_dir ]
then
    mkdir -p $out_dir
fi

shift
additional_args="$@"

#Set the AlphaFold2 prediction command
CMD="/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch $additional_args $input_file $out_dir"
#CMD="/home/aljubetic/AF2/CF2/bin/colabfold_batch --num-recycle 12 --model-type AlphaFold2-multimer-v2 $input_file $out_dir"

echo running AF2 command:
echo $CMD
# Run the prediction command
$CMD
