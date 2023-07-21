# ClusterProteinDesign

This is a personal repository for running protein design scripts on Slurm. The scripts are designed to use various tools such as RFdiffusion, ProteinMPNN, AlphaFold2, and Rosetta.

### Usage:
To use the scripts from this repository, one must install several different things:
- RFdiffusion
- ColabDesign
- PyRosetta/Rosetta
- ColabFold
- RoseTTAfold 2
- ESMfold

### Pipelines
- [ColabFold AF2 predictions](/pipelines/colabfold_af2)
    - [Predict single fasta/a3m file](/pipelines/colabfold_af2/README.md#Predict-single-fasta/a3m-file)
    - [Predict pdb files and compare](/pipelines/colabfold_af2)
- [Protein design with ProteinMPNN and AF2](/pipelines/mpnn_af2)
    - MPNN-AF2 pipeline:
        - [ProteinMPNN](/pipelines/mpnn_af2)
        - [AF2 prediction (ColabFold)](/pipelines/mpnn_af2)
        - [Structure analysis](/pipelines/mpnn_af2)
- [Binder design with RFdiffusion, ProteinMPNN and AF2](/pipelines/binder_design)
    - RF-MPNN-AF2 pipeline:
        - Binder docking:
            - [Random binder docking](/pipelines/binder_design)
            - [RFdiffusion scaffold docking](/pipelines/binder_design)
            - [RFdiffusion random scaffold docking (WIP)]()
        - [Binder optimization](/pipelines/binder_design)
        - [Binder analysis](/pipelines/binder_design)
        - [Binder filter best (WIP)]()