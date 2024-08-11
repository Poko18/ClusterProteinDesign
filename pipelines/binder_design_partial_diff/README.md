# Binder design
Binder design using partial diffusion.

![binders](../binder_design/imgs/binder_design.png)

## Round 1 - binder scaffold docking
Idea is to get get initial binder sequence, that binds to the target protein. Start with `01_binder_scaffold_dock.ipynb` for docking scaffolds (usually 3HB)

## Round 1.5 - analyze diffused binders and optimize them with ColabDesign
Get stats about diffused binders from the second part of `01_binder_scaffold_dock.ipynb` notebook and then run ColabDesign optimization on the best ones with `015_cd_binder_opt.ipynb`. Optimized binders can serve as great complementary hits for partial diffusion in notebook `02_partial_diff.ipynb`

## Round 2 - use partial diffusion to generate more diverse (designable) backbones
Use `02_partial_diff.ipynb` notebook to run partial diffusion on best hits from the first notebook

## Round 3 - in parallel calculate all of the metrics with PyRosetta
Use `03_binder_analysis.ipynb` notebook.

## Round 4 - Filter designs and repeat AF2 prediction with ColabFold without Initial Guess
Use `04_binder_filter.ipynb` notebook.
