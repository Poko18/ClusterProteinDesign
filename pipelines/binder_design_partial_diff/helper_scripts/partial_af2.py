import os
import sys
import glob
import fcntl
import time
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af import mk_af_model
from colabdesign.mpnn.model import aa_order
from Bio import PDB
import pandas as pd
import numpy as np
import argparse

import pyrosetta
from pyrosetta import init, pose_from_file, pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import LayerSelector, ChainSelector, AndResidueSelector
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover

init()

parser = argparse.ArgumentParser(description='Run AFDesign with MPNN sampling')
parser.add_argument('pdb', type=str, help='input pdb file path')
parser.add_argument('output_folder', type=str, help='output folder path')
parser.add_argument('target_chains', type=str, help='target chain IDs (comma-separated)')
parser.add_argument('binder_chains', type=str, help='binder chain IDs (comma-separated)')
parser.add_argument('--relax', action='store_true', help='relax the structure before designing)')
parser.add_argument('--num_seqs', type=int, default=1000, help='number of sequences to generate (default: 1000)')
parser.add_argument('--sampling_temp', type=float, default=0.1, help='mpnn sampling temperature (default: 0.1)')
parser.add_argument('--mpnn_batch', type=int, default=8, help='mpnn batch size (default: 8)')
parser.add_argument('--num_recycles', type=int, default=12, help='number of repacking cycles (default: 12)')
parser.add_argument('--results_dataframe', type=str, help='save results')
parser.add_argument('--save_best_only', action='store_true', help='save only the best structures')
args = parser.parse_args()


pdb = args.pdb
output_folder = args.output_folder
target_chains = args.target_chains.split(',')
binder_chains = args.binder_chains.split(',')
num_seqs = args.num_seqs
sampling_temp = args.sampling_temp
mpnn_batch = args.mpnn_batch
num_recycles = args.num_recycles
results_dataframe = args.results_dataframe

assert num_seqs % mpnn_batch == 0, "num_seqs must be divisible by mpnn_batch"

### Functions ###

def copy_b_factors(source_pose, target_pose):
    """Copy B factors from source_pose to target_pose."""
    for resid in range(1, source_pose.total_residue() + 1):
        if source_pose.residue(resid).is_protein():
            # Get the B factor of the first heavy atom in the residue
            bfactor = source_pose.pdb_info().bfactor(resid, 1)
            for atom_id in range(1, source_pose.residue(resid).natoms() + 1):
                target_pose.pdb_info().bfactor(resid, atom_id, bfactor)

def relax_structure(pdb_file, relaxed_pdb_path):
    if not os.path.exists(relaxed_pdb_path):
        try:
            # Generate pose from PDB file
            pose = pose_from_pdb(pdb_file)
            start_pose = pose.clone()

            # Generate MoveMap
            movemap = MoveMap()
            movemap.set_chi(True) # Enable sidechain movement
            movemap.set_bb(True) # Enable backbone movement
            movemap.set_jump(False) # Disable whole chain movement

            # Run FastRelax
            fastrelax = FastRelax()
            scorefxn = get_fa_scorefxn()
            fastrelax.set_scorefxn(scorefxn)
            fastrelax.set_movemap(movemap)
            fastrelax.max_iter(200) # Reduced iterations
            fastrelax.min_type("lbfgs_armijo_nonmonotone")
            fastrelax.constrain_relax_to_start_coords(True)
            fastrelax.apply(pose)

            # Align relaxed structure to original trajectory
            align_mover = AlignChainMover()
            align_mover.source_chain(0)
            align_mover.target_chain(0)
            align_mover.pose(start_pose)
            align_mover.apply(pose)

            # Copy B factors from start_pose to relaxed pose
            copy_b_factors(start_pose, pose)

            # Output relaxed and aligned PDB
            pose.dump_pdb(relaxed_pdb_path)
        except Exception as e:
            print(f"Error during relaxation: {e}")

def select_surface_binder_residues(pdb, binder_chains):
    pose = pose_from_file(pdb)

    layer_selector = LayerSelector()
    layer_selector.set_layers(False, False, True)

    chain_selector = ChainSelector(",".join(binder_chains))

    combined_selector = AndResidueSelector(layer_selector, chain_selector)
    selected_residues = np.array(combined_selector.apply(pose))

    return selected_residues[list(chain_selector.apply(pose))]

def write_df_to_csv(df, output_path):
    mode = 'a' if os.path.isfile(output_path) else 'w'
    with open(output_path, mode) as f:
        fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
        time.sleep(np.random.uniform(0, 0.05))  # introduce random delay
        df.to_csv(f, header=(mode == 'w'), index=False)  # write header only if in write mode
        fcntl.flock(f, fcntl.LOCK_UN)  # release the lock

### Initialize AF2 and ProteinMPNN ###
mpnn_terms = ["score", "seq"]
af_terms = ["plddt","i_ptm","i_pae","rmsd"]
other_terms = ["model_path", "input_pdb"]

af_model = mk_af_model(protocol="binder",
                       initial_guess = True,
                       best_metric = "rmsd",
                       use_multimer = True,
                       data_dir="/home/tsatler/projects/AFdesign_playground",
                       model_names = ["model_1_multimer_v3"])

af_model.prep_inputs(pdb,
                    target_chain=",".join(target_chains),
                    binder_chain=",".join(binder_chains),
                    rm_aa="C")

mpnn_model = mk_mpnn_model()
mpnn_model.get_af_inputs(af_model)

### MPNN makes binders always on the second place.. thats why we modify the bias at the end ###

### Add alanine anti bias
binder_len = af_model._len
for residue_index in range(binder_len):
        mpnn_model._inputs["bias"][(residue_index - binder_len), aa_order["A"]] += -0.5


### Add surface residue bias to MPNN model ###
surface_aa = "DEHKNQRSTY"
polar_bias_value = 0.8
binder_sel = select_surface_binder_residues(pdb, binder_chains)

print(f"Surface binder residues: {[i+1 for i, x in enumerate(binder_sel) if x]}")
print(f"Adding bias ({polar_bias_value}) to surface residues: {surface_aa}")

assert len(binder_sel) == af_model._len, "Binder selection length does not match AF model length"

# Add bias to mpnn model
binder_len = len(binder_sel)
for residue_index in range(binder_len):
    if binder_sel[residue_index]:
        for aa in surface_aa:
            mpnn_model._inputs["bias"][(residue_index - binder_len), aa_order[aa]] += polar_bias_value


### Apply relaxation
if args.relax:
    relaxed_pdb_path = pdb.replace('.pdb', '_relaxed.pdb')
    relax_structure(pdb, relaxed_pdb_path)
    pdb = relaxed_pdb_path


### Sample sequences with MPNN ###
mpnn_out = mpnn_model.sample(num=num_seqs//mpnn_batch, batch=mpnn_batch, temperature=sampling_temp)

# Initialize mpnn_out dictionary with keys from af_terms
for k in af_terms + other_terms:
    if k not in mpnn_out:
        mpnn_out[k] = []



### Run AF2 on MPNN sequences ###
print("Starting AF2 structure prediction on MPNN sequences")

# Get basename and create output folder
pdb_name = os.path.splitext(os.path.basename(pdb))[0]
output_path = os.path.join(output_folder, "af2")
os.makedirs(output_path, exist_ok=True)

# Iterate over MPNN sequences
for n in range(num_seqs):
    seq = mpnn_out["seq"][n][-af_model._len:]
    print(f"Designing sequence {n+1}/{num_seqs} with sequence {seq}")

    # Predict structure with AF2
    af_model.predict(seq=seq, num_recycles=num_recycles, num_models=1, verbose=False)
    print(af_model.aux["log"])

    # Append logs to mpnn_out
    for term in af_terms:
        mpnn_out[term].append(af_model.aux["log"][term])

    # Adjust PAE values if they exist in the output
    if "i_pae" in mpnn_out:
        mpnn_out["i_pae"][-1] *= 31
    if "pae" in mpnn_out:
        mpnn_out["pae"][-1] *= 31

    # Save the current model
    current_model_path = f"{output_path}/{pdb_name}_{n}.pdb"
    mpnn_out["model_path"].append(current_model_path)
    mpnn_out["input_pdb"].append(pdb)

    # Save PDB file based on save_best_only flag
    save_pdb = not args.save_best_only or (mpnn_out["plddt"][n] > 0.7 and mpnn_out["rmsd"][n] < 3)
    if save_pdb:
        af_model.save_current_pdb(current_model_path)
    af_model._save_results(save_best=args.save_best_only, verbose=False)

    # Increment model counter
    af_model._k += 1

# Generate model paths for all sequences
model_paths = [f"{output_path}/{pdb_name}_{n}.pdb" for n in range(num_seqs)]

# Combine labels for data extraction
all_labels = mpnn_terms + af_terms + other_terms

# Extract data for all sequences
data = [[mpnn_out[label][n] for label in all_labels] for n in range(num_seqs)]
all_labels[0] = "mpnn"



### Write results to CSV ###
df = pd.DataFrame(data, columns=all_labels)

# Write df to CSV
output_path_all = f'{results_dataframe}/af2_results_all.csv' if results_dataframe else f'{output_folder}/af2_results_all.csv'
write_df_to_csv(df, output_path_all)

# Filter df for best results and write to CSV
df_best = df[(df["rmsd"] < 3) & (df["plddt"] > 0.7)]
output_path_best = f'{results_dataframe}/af2_best.csv' if results_dataframe else f'{output_folder}/af2_best.csv'
write_df_to_csv(df_best, output_path_best)
