# ProteinMPNN-AF2 pipeline
Basic pipeline for protein design with ProteinMPNN and AF2 (ColabFold). It saves results in output folder

### Usage
#### 1. ProteinMPNN
Before running the script, prepare input file/folder and setup the script parameters.
```
SBATCH 01_mpnn.sh <input_file (pdb file/ folder with pdbs)>
```
#### 2. AF2 prediction (ColabFold)
Follow the `02_prep_af2_inputs.ipynb` notebook to prepare the inputs for AF2 prediction. Make sure to specify all the neccessary parameters. For instance, for prediction without MSA add `--msa-mode single_sequence` flag.

To do:
- [ ] notebook that uses only target MSA (and not for de novo designed binder)

#### 3. Basic analysis of predicted structures
Example of basic analysis in `03_analyze_results.ipynb` notebook.

To do:
- [ ] organize the code for metric calculation
