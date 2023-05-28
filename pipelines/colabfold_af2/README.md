# Colabfold scripts
Basic scripts to run AF2 with ColabFold.

## Predict single fasta/a3m file
### Usage
``` 
SBATCH predict_single_AF2.sh <input_file (fasta,a3m)> <additional_args (--msa-mode single sequence --num-recycle 16)>
``` 

### Example
``` 
SBATCH predict_single_AF2.sh examples/3HB.fasta --num-recycle 16
``` 

## Predict pdb files and compare
Adjust the array range accordingly
### Usage
Note: ignore binder args if interaction analysis is not required
``` 
sbatch predict_pdbs_AF2.sh <input_file (pdb,folder with pdbs)> <binder args (binder,binder-second)> <additional_args (--msa-mode single sequence --num-recycle 16)>
``` 

### Example
``` 
sbatch predict_pdbs_AF2.sh examples/3HB_folder/ binder-second --num-recycle 16
``` 