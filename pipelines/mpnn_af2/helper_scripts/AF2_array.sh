#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Activate the AlphaFold2 environment
. /home/aljubetic/bin/set_up_AF2.3.sh

input=$1
input_files=($input/*.fasta)
file=${input_files[$SLURM_ARRAY_TASK_ID]}
out_dir=$input/AF2out_$SLURM_ARRAY_TASK_ID

# Create the output directory
if [ ! -d $out_dir ]
then
    mkdir -p $out_dir
fi

shift
additional_args="$@"

#Set the AlphaFold2 prediction command
CMD="/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/python /home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch $additional_args $file $out_dir"

echo running AF2 command:
echo $CMD

# Run the prediction command
$CMD
