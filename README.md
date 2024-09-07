# ClusterProteinDesign

This is a personal repository for running protein design scripts on Slurm. The scripts are designed to use various tools such as RFdiffusion, ProteinMPNN, AlphaFold2, Rosetta, etc..

### Usage:
To use the scripts from this repository, one must install several different things:
- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)
- [ColabDesign](https://github.com/sokrypton/ColabDesign)

What worked for me was installing cuda before colabdesign:

```
conda create -n colabdesign python=3.9 pip
conda activate colabdesign
pip install --upgrade "jax[cuda11_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1

#Downgrade to match cluster version of cuda
pip install nvidia-cuda-nvcc-cu11==11.7.99
pip install nvidia-cuda-runtime-cu11==11.7.99
```

- PyRosetta/Rosetta
- [ColabFold](https://github.com/sokrypton/ColabFold/tree/main)
- [RoseTTAfold2](https://github.com/Poko18/RoseTTAFold2)
> **Note:** I linked my fork that has some additional scripts
- [ESMfold](https://github.com/facebookresearch/esm)

### Pipelines
- [ColabFold AF2 predictions](/pipelines/colabfold_af2)
    - [Predict single fasta/a3m file](/pipelines/colabfold_af2/#Predict-single-fasta/a3m-file)
    - [Predict pdb files and compare](/pipelines/colabfold_af2/#Predict-pdb-files-and-compare)
- [Protein design with ProteinMPNN and AF2](/pipelines/mpnn_af2)
    - MPNN-AF2 pipeline:
        - [ProteinMPNN](/pipelines/mpnn_af2#1.-ProteinMPNN)
        - [AF2 prediction (ColabFold)](/pipelines/mpnn_af2#2.-AF2-prediction-(ColabFold))
        - [Structure analysis](/pipelines/mpnn_af2#3.-Basic-analysis-of-predicted-structures)

    - [Fix conserved positions MPNN-AF2 pipeline](/pipelines/msa_mpnn_af2):
        - based on article: [Sumida2024 - Improving Protein Expression, Stability, and Function with ProteinMPNN](https://doi.org/10.1021/jacs.3c10941)

- [Binder design with RFdiffusion, ProteinMPNN and AF2](/pipelines/binder_design)
    - RF-MPNN-AF2 pipeline:
        - [Binder docking](/pipelines/binder_design#Round-1---binder-scaffold-docking):
            - [Random binder docking](/pipelines/binder_design#1a-Random-binder-docking)
            - [RFdiffusion scaffold docking](/pipelines/binder_design#1b-RFdiffusion-selected-scaffold-docking)
            - [RFdiffusion random binder scaffold docking](/pipelines/binder_design#1c-RFdiffusion-random-binder-scaffold-docking)
        - [Binder optimization](/pipelines/binder_design#Round-2---binder-optimization)
        - [Binder analysis](/pipelines/binder_design#Round-3---binder-analysis)
        - [Binder filtering](/pipelines/binder_design#Round-4---binder-filtering-and-sequence-clustering)

- [Binder design with partial diffusion](/pipelines/binder_design_partial_diff)

- [Motif scaffolding](/pipelines/motif_scaffolding)
