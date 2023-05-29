# ClusterProteinDesign

This is a personal repository for running protein design scripts on Slurm. The scripts are designed to use various tools such as RFdiffusion, ProteinMPNN, AlphaFold2, and Rosetta.

### Usage:
To use the scripts from this repository, one must install several different things:
- RFdiffusion
- ColabDesign
- PyRosetta/Rosetta
- ColabFold

### Pipelines
[ColabFold AF2 predictions](/pipelines/colabfold_af2)

[Protein design with ProteinMPNN and AF2](/pipelines/mpnn_af2)

[Binder design (WIP)]

#### TO DO:
- [x] Add proteinMPNN/FastRelax step between Folddocking and ProteinMPNN/AF2 to get better protein backbone
- [x] Thread and FastRelax proteinMPNN sequences to RFdiffused binder backbones before ProteinMPNN/AF2 metrics
    - currently best approach was to use original binder sequence (LCB3 for example)
    - [ ] Prepare sequences that fold into a nice binder for all scaffolds!
    - [ ] Redesign only interacting residues and fix other for proteinMPNN
    - idea - whole design pipeline --> RFdif with scaffolding, ProteinMPNN, thread (10?) and relax, proteinMPNN newly relaxed backbones, af2 predict, best rmsd / plddt binders --> proteinMPNN (much more), af2 predict and calculate metrics for filtering
- [ ] write cleaner metric scripts and more of them
- [ ] implement saving all generated pdbs in one file (silent?)
- [ ] get sequences for all scaffolds
- [ ] make metric functions more general and group them in common file
- [ ] BINDER_DESIGN: binder search + filter notebook + final binder with metrics
- [x] QUICK_ProteinMPNN: (notebook for fixing/tieing residues) proteinMPNN + AF2 (with/without MSA)
- [x] QUICK_AF2_test: folder with pdbs/file with sequences --> AF2 (with/without MSA)(+ESMfold?) --> RMSD to original/ plddt metrics