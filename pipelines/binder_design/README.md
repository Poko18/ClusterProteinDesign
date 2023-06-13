# Binder design
Binder design pipeline, that worked well for me, also for harder targets.

### Round 1 - binder scaffold docking
- [x] RFdiffusion for docking of 3HB scaffolds near hotspot residues
- [x] threading the poly-glycine backbone with scaffold sequence and FastRelax (yields better ProteinMPNN sequences)
- [x] ProteinMPNN and AF2 prediction of binders

### Round 2 - binder optimization
- Filtering of docked scaffolds
- Multiple iterations of ProteinMPNN and AF2 prediction of binders

### Round 3 - binder analysis
- calculate metrics