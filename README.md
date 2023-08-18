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
    - [Predict single fasta/a3m file](/pipelines/colabfold_af2/#Predict-single-fasta/a3m-file)
    - [Predict pdb files and compare](/pipelines/colabfold_af2/#Predict-pdb-files-and-compare)
- [Protein design with ProteinMPNN and AF2](/pipelines/mpnn_af2)
    - MPNN-AF2 pipeline:
        - [ProteinMPNN](/pipelines/mpnn_af2#1.-ProteinMPNN)
        - [AF2 prediction (ColabFold)](/pipelines/mpnn_af2#2.-AF2-prediction-(ColabFold))
        - [Structure analysis](/pipelines/mpnn_af2#3.-Basic-analysis-of-predicted-structures)
- [Binder design with RFdiffusion, ProteinMPNN and AF2](/pipelines/binder_design)
    - RF-MPNN-AF2 pipeline:
        - [Binder docking](/pipelines/binder_design#Round-1---binder-scaffold-docking):
            - [Random binder docking](/pipelines/binder_design#1a-Random-binder-docking)
            - [RFdiffusion scaffold docking](/pipelines/binder_design#1b-RFdiffusion-selected-scaffold-docking)
            - [RFdiffusion random binder scaffold docking](/pipelines/binder_design#1c-RFdiffusion-random-binder-scaffold-docking)
        - [Binder optimization](/pipelines/binder_design#Round-2---binder-optimization)
        - [Binder analysis](/pipelines/binder_design#Round-3---binder-analysis)
        - [Binder filter best (WIP)]()