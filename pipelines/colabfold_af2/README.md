# Colabfold script
Basic script to run AF2 with ColabFold. It saves results in output folder

### Usage
``` 
SBATCH 01_AF2.sh <input_file (fasta,a3m)> <additional_args (--msa-mode single sequence --num-recycle 16)>
``` 

### Example
``` 
SBATCH 01_AF2.sh examples/3HB.fasta --num-recycle 16
``` 

### TO DO
- [ ] add a script that would parse folder of pdbs, predict and do rmsd (notebook? or both.. so it is not running)
