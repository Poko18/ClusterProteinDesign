# Binder design
Binder design pipeline, that worked well for me, also for harder targets.

### Round 1 - binder scaffold docking
- RFdiffusion for docking of 3HB scaffolds near hotspot residues
- threading the poly-glycine backbone with scaffold sequence and FastRelax (yields better ProteinMPNN sequences)
- ProteinMPNN and AF2 prediction of binders

### Filtering
- Mostly RMSD filtering..?

### Round 2 - binder optimization (iterating?)
- ProteinMPNN and AF2 prediction of binders
- Metrics calculation
- Filtering?

