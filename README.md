# ClusterProteinDesign

This is a personal repository for running protein design scripts on Slurm. The scripts are designed to use various tools such as RFdiffusion, ProteinMPNN, AlphaFold2, and Rosetta.

### Usage:
To use these scripts, one must install several different things:
- RFdiffusion
- ColabDesign
- PyRosetta/Rosetta

### TO DO:
- [x] Add proteinMPNN/FastRelax step between Folddocking and ProteinMPNN/AF2 to get better protein backbone
- [x] Thread and FastRelax proteinMPNN sequences to RFdiffused binder backbones before ProteinMPNN/AF2 metrics
    - currently best approach was to use original binder sequence (LCB3 for example)
    - [ ] Prepare sequences that fold into a nice binder for all scaffolds!
    - [ ] Redesign only interacting residues and fix other for proteinMPNN
    - idea - whole design pipeline --> RFdif with scaffolding, ProteinMPNN, thread (10?) and relax, proteinMPNN newly relaxed backbones, af2 predict, best rmsd / plddt binders --> proteinMPNN (much more), af2 predict and calculate metrics for filtering
- [ ] write cleaner metric scripts and more of them
- [ ] implement saving all generated pdbs in one file (silent?)
