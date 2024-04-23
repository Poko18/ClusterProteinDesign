import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from Bio import PDB
from colabdesign.af import mk_af_model
from colabdesign.mpnn import mk_mpnn_model
from pyrosetta import *
from pyrosetta.rosetta import *

parser = argparse.ArgumentParser(description="Run AFDesign with MPNN sampling")
parser.add_argument("pdb", type=str, help="input pdb file path")
parser.add_argument("output_folder", type=str, help="output folder path")
parser.add_argument("target_chains", type=str, help="target chain IDs (comma-separated)")
parser.add_argument("binder_chains", type=str, help="binder chain IDs (comma-separated)")
parser.add_argument(
    "--num_recycles",
    type=int,
    default=12,
    help="number of repacking cycles (default: 12)",
)
parser.add_argument(
    "--num_seqs",
    type=int,
    default=1000,
    help="number of sequences to generate (default: 1000)",
)
parser.add_argument(
    "--num_filt_seq",
    type=int,
    default=50,
    help="number of sequences to keep after filtering (default: 50)",
)
parser.add_argument(
    "--sampling_temp",
    type=float,
    default=0.15,
    help="sampling temperature (default: 0.15)",
)
parser.add_argument("--results_dataframe", type=str, help="save results")
parser.add_argument(
    "--break_design",
    action="store_true",
    help="stop making new sequences if rmsd < 2 and plddt > 0.9",
)
parser.add_argument("--save_best_only", action="store_true", help="save only the best structures")
parser.add_argument(
    "--thread_sequence",
    action="store_true",
    help="Threads RFdiffusion designs with ProteinMPNN sequence",
)
parser.add_argument("--scaffold_folder", type=str, help="Scaffold folder to do rmsd")
args = parser.parse_args()

# Initialize PyRosetta
# init("-beta_nov16 -in:file:silent_struct_type binary")
init("-beta_nov16 -holes:dalphaball")
xml = "/home/tsatler/RFdif/ClusterProteinDesign/testing/mpnnthreading/RosettaFastRelaxUtil.xml"
objs = protocols.rosetta_scripts.XmlObjects.create_from_file(xml)

# Load the movers we will need

FastRelax = objs.get_mover("FastRelax")

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
states = len(alpha_1)
alpha_3 = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "GAP",
]
aa_1_3 = {a: b for a, b in zip(alpha_1, alpha_3)}


def thread_mpnn_seq(pose, binder_seq):
    rsd_set = pose.residue_type_set_for_pose(core.chemical.FULL_ATOM_t)

    for resi, mut_to in enumerate(binder_seq):
        resi += 1  # 1 indexing
        name3 = aa_1_3[mut_to]
        new_res = core.conformation.ResidueFactory.create_residue(rsd_set.name_map(name3))
        pose.replace_residue(resi, new_res, True)

    return pose


def relax_pose(pose):
    FastRelax.apply(pose)
    return pose


def get_sequence_from_pdb(pdb_file_path):
    # Create a PDB parser object
    parser = PDB.PDBParser()
    structure = parser.get_structure("pdb_structure", pdb_file_path)
    model = structure[0]
    sequence = ""

    # Iterate over all chains in the model
    for chain in model:
        # Iterate over all residues in the chain
        for residue in chain:
            # Get the residue name (3-letter code)
            residue_name = residue.get_resname()

            # Ignore non-standard residues
            if PDB.is_aa(residue_name):
                # Append the residue name to the sequence
                sequence += PDB.Polypeptide.three_to_one(residue_name)

    return sequence


def align_structures(pdb1, pdb2, save_aligned=False):
    """Take two structure and superimpose pdb1 on pdb2"""
    import Bio.PDB

    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    # Get the structures
    ref_structure = pdb_parser.get_structure("ref", pdb1)
    sample_structure = pdb_parser.get_structure("sample", pdb2)

    aligner = Bio.PDB.cealign.CEAligner()
    aligner.set_reference(ref_structure)
    aligner.align(sample_structure)

    # Save aligned coordinates
    if save_aligned:
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure)
        io.save(pdb2)

    return aligner.rms


pdb = args.pdb
output_folder = args.output_folder
target_chains = args.target_chains.split(",")
binder_chains = args.binder_chains.split(",")
num_recycles = args.num_recycles
num_seqs = args.num_seqs
num_filt_seq = args.num_filt_seq
sampling_temp = args.sampling_temp
results_dataframe = args.results_dataframe
scaffold_folder = args.scaffold_folder

# Threading and relaxing diffused backbone!
if args.thread_sequence:

    # make thread folder
    pdb_name = pdb.split("/")[-1].split(".pdb")[0]
    pdb_dir_path = os.path.dirname(pdb)
    thread_dir = os.path.join(pdb_dir_path, "threaded")
    os.makedirs(thread_dir, exist_ok=True)

    # path to save

    print("threading RF design")
    threaded_pose = pose_from_pdb(pdb)

    if scaffold_folder:
        scaff_pdb_path = glob.glob(f"{scaffold_folder}/*pdb")[0]
        thread_seq = get_sequence_from_pdb(scaff_pdb_path)
    else:
        mpnn_model = mk_mpnn_model(weights="soluble")
        mpnn_model.prep_inputs(pdb_filename=pdb, chain="A", verbose=False)
        out = mpnn_model.sample(num=10 // 10, batch=10, temperature=0.15, verbose=True)
        score_list = out["score"]
        sorted_index = sorted(range(len(score_list)), key=lambda k: score_list[k], reverse=False)
        for key, value in out.items():
            sorted_value = [value[i] for i in sorted_index][:10]  # slice the sorted list to only keep the first x elements
            out[key] = sorted_value
        thread_seq = out["seq"][0]

    threaded_pose = thread_mpnn_seq(threaded_pose, thread_seq)
    threaded_pose = relax_pose(threaded_pose)

    thread_dir = os.path.join(pdb_dir_path, "threaded")
    thread_pdb_name = f"{pdb_name}_threaded.pdb"
    save_thread = os.path.join(thread_dir, thread_pdb_name)
    threaded_pose.dump_pdb(save_thread)

    pdb = save_thread
    print(f"Threaded protein: {pdb}")


### Initialize ###
protocol = "binder"
if protocol == "binder":
    af_terms = ["plddt", "i_ptm", "i_pae", "rmsd"]

af_model = mk_af_model(
    protocol="binder",
    initial_guess=True,
    best_metric="rmsd",
    use_multimer=True,
    data_dir="/home/tsatler/projects/AFdesign_playground",
    model_names=["model_1_multimer_v3"],
)

af_model.prep_inputs(
    pdb,
    target_chain=",".join(target_chains),
    binder_chain=",".join(binder_chains),
    rm_aa="C",
)


# Generate ProteinMPNN (100-1000)
print("starting ProteinMPNN")
mpnn_model = mk_mpnn_model()
mpnn_model.get_af_inputs(af_model)
out = mpnn_model.sample(num=num_seqs // 10, batch=10, temperature=sampling_temp)

# Filter top X sequences based on score
# Extract the score list from the dictionary and create a sorted index
score_list = out["score"]
sorted_index = sorted(range(len(score_list)), key=lambda k: score_list[k], reverse=True)  # False --> lower score is first

# Iterate through the dictionary and sort each value (list) using the sorted index
for key, value in out.items():
    sorted_value = [value[i] for i in sorted_index][:num_filt_seq]  # slice the sorted list to only keep the first x elements
    out[key] = sorted_value

for k in af_terms:
    out[k] = []  # ?

if "model_path" not in out:
    out["model_path"] = []

if "input_pdb" not in out:
    out["input_pdb"] = []

if "binder-rmsd" not in out:
    out["binder-rmsd"] = []

print(out.keys())

pdb_name = pdb.split("/")[-1].split(".pdb")[0]
loc = f"{output_folder}/af2"
os.system(f"mkdir -p {loc}")

# AF2
print("starting AF2")
# with open(f"{loc}/design.fasta","w") as fasta:
for n in range(num_filt_seq):
    seq = out["seq"][n][-af_model._len :]
    af_model.predict(seq=seq, num_recycles=num_recycles, num_models=1, verbose=False)
    print(af_model.aux["log"])
    for t in af_terms:
        out[t].append(af_model.aux["log"][t])
    if "i_pae" in out:
        out["i_pae"][-1] = out["i_pae"][-1] * 31
    if "pae" in out:
        out["pae"][-1] = out["pae"][-1] * 31

    current_model = f"{loc}/{pdb_name}_{n}.pdb"
    out["model_path"].append(current_model)
    out["input_pdb"].append(pdb)

    if scaffold_folder:
        af_model.save_current_pdb(f"{current_model}")
        binder_rmsd = align_structures(scaff_pdb_path, current_model)
        out["binder-rmsd"].append(binder_rmsd)
        os.remove(current_model)
    else:
        out["binder-rmsd"].append(None)

    if args.save_best_only:
        if out["plddt"][n] > 0.7 or out["rmsd"][n] < 3:
            af_model.save_current_pdb(f"{current_model}")
        af_model._save_results(save_best=True, verbose=False)
    else:
        af_model.save_current_pdb(f"{current_model}")
        af_model._save_results(save_best=True, verbose=False)
    af_model._k += 1

    if args.break_design:
        if out["plddt"][n] > 0.9 and out["rmsd"][n] < 1:
            print("Breaking out of loop due to --break_design flag")
            num_filt_seq = n  # update the number to save later
            break

        # fasta.write(line+"\n")

# af_model.save_pdb(f"{loc}/best.pdb")

model_paths = [f"{loc}/{pdb_name}_{n}.pdb" for n in range(num_filt_seq)]

labels = ["score"] + af_terms + ["seq"] + ["model_path"] + ["binder-rmsd"] + ["input_pdb"]

# data = [[out[k][n] for k in labels] + [model_paths[n], pdb] for n in range(num_filt_seq)]

data = [[out[k][n] for k in labels] for n in range(num_filt_seq)]


# add additional labels
# labels = ["score"] + af_terms + ["seq"] + ["model_path"] + ["input_pdb"]
labels[0] = "mpnn"


import fcntl
import time

df = pd.DataFrame(data, columns=labels)

if results_dataframe:
    output_path = f"{results_dataframe}/af2_results.csv"
else:
    output_path = f"{output_folder}/af2_results.csv"

if os.path.isfile(output_path):
    with open(output_path, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
        time.sleep(np.random.uniform(0, 0.05))
        df.to_csv(f, header=False, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)  # release the lock
else:
    with open(output_path, "w") as f:
        fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
        time.sleep(np.random.uniform(0, 0.05))
        df.to_csv(f, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)  # release the lock

# Save also the best
df_best = df[df["rmsd"] < 3]

if results_dataframe:
    output_path = f"{results_dataframe}/af2_best.csv"
else:
    output_path = f"{output_folder}/af2_best.csv"

if os.path.isfile(output_path):
    with open(output_path, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
        time.sleep(np.random.uniform(0, 0.05))
        df_best.to_csv(f, header=False, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)  # release the lock
else:
    with open(output_path, "w") as f:
        fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
        time.sleep(np.random.uniform(0, 0.05))
        df_best.to_csv(f, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)  # release the lock
