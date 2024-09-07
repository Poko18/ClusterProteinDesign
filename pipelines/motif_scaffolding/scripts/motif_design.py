"""
Use ColabDesign library to validate diffusied protein scaffolds
"""

import argparse
import fcntl
import logging
import math
import os
import time
from string import ascii_uppercase, ascii_lowercase

import numpy as np
import pandas as pd
from utils import parse_pdb

from colabdesign.af import mk_af_model
from colabdesign.mpnn import mk_mpnn_model


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("binder_cycling.py")
logging.getLogger("jax").setLevel(logging.ERROR)
logging.getLogger("tensorflow").setLevel(logging.ERROR)
logging.getLogger("matplotlib").setLevel(logging.ERROR)


parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument("pdb", type=str, help="input pdb file path")
parser.add_argument("output_folder", type=str, help="output folder path")
parser.add_argument("contigs", type=str, help="RFdiffusion contigs")
parser.add_argument(
    "--input_pdb",
    type=str,
    help="Input pdb that was used with RFdiffusion",
)
parser.add_argument(
    "--input_pdb_as",
    type=str,
    help="Input pdb active site",
)
parser.add_argument(
    "--num_seqs",
    type=int,
    default=80,
    help="number of sequences to generate (default: 80)",
)
parser.add_argument(
    "--sampling_temp",
    type=float,
    default=0.1,
    help="mpnn sampling temperature (default: 0.1)",
)
parser.add_argument(
    "--num_recycles",
    type=int,
    default=12,
    help="number of repacking cycles (default: 12)",
)
parser.add_argument(
    "--rm_aa", type=str, default="C", help="residue to remove from the design"
)
parser.add_argument(
    "--mpnn_batch", type=int, default=8, help="mpnn batch size (default: 8)"
)
parser.add_argument("--results_dataframe", type=str, help="save results")
parser.add_argument(
    "--save_best_only", action="store_true", help="save only the best structures"
)
parser.add_argument(
    "--initial_guess",
    action="store_true",
    help="use initial guess for alphafold2 validation",
)
parser.add_argument(
    "--use_multimer",
    action="store_true",
    help="use multimer weights for alphafold2 validation",
)
parser.add_argument(
    "--use_soluble", action="store_true", help="use soluble weights for mpnn"
)
args = parser.parse_args()

pdb = args.pdb
output_folder = args.output_folder
rf_contigs = args.contigs
input_pdb = args.input_pdb
input_pdb_as = args.input_pdb_as
num_seqs = args.num_seqs
sampling_temp = args.sampling_temp
mpnn_batch = args.mpnn_batch
num_recycles = args.num_recycles
rm_aa = args.rm_aa
results_dataframe = args.results_dataframe
save_best_only = args.save_best_only
initial_guess = args.initial_guess
use_multimer = args.use_multimer
use_soluble = args.use_soluble


### Functions ###
def write_df_to_csv(df, output_path):
    """Save/append dataframe to CSV."""
    mode, header = ("a", False) if os.path.isfile(output_path) else ("w", True)
    with open(output_path, mode) as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        time.sleep(np.random.uniform(0, 0.05))
        df.to_csv(f, header=header, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)


def parse_range(_range) -> tuple:
    """Return the start and end residue in a range str."""
    if "-" in _range:
        s, e = _range.split("-")
    else:
        s, e = _range, _range
    return int(s), int(e)


def parse_contig(contig) -> tuple:
    """
    Return the chain, start and end residue in a contig or gap str.
    Ex:
    'A4-8' --> 'A', 4, 8
    'A5'   --> 'A', 5, 5
    '4-8'  --> None, 4, 8
    'A'    --> 'A', None, None
    """
    if contig[0].isalpha():  # contig
        ch = contig[0]
        if len(contig) > 1:
            s, e = parse_range(contig[1:])
        else:
            s, e = None, None
    else:  # gap
        ch = None
        s, e = parse_range(contig)
    return ch, s, e


def expand(mask_str) -> list:
    """
    Ex: '2,A3-5,3' --> [None, None, (A,3), (A,4), (A,5), None, None, None]
    """
    expanded = []
    for l in mask_str.split(","):
        ch, s, e = parse_contig(l)
        if ch:
            expanded += [(ch, res) for res in range(s, e + 1)]
        else:
            expanded += [None for _ in range(s)]
    return expanded


def get_designed_motif_offset(inptu_pdb_as: str, contig: str) -> str:
    """Calculate the offset values based on the given test_as and contig strings."""
    contig_split = contig.split("/")
    expanded_split = expand(",".join(contig_split))
    expanded_dict = {i: aa for i, aa in enumerate(expanded_split) if aa}

    offset_as = []
    parsed_ranges = [parse_contig(l) for l in inptu_pdb_as.split(",")]

    for ch, s, e in parsed_ranges:
        for i, aa in expanded_dict.items():
            if aa[0] == ch and s <= aa[1] <= e:
                if s == e:
                    offset_as.append(f"{aa[0]}{i}")
                else:
                    add = e - s
                    offset_as.append(f"{aa[0]}{i}-{i+add}")
                break

    return ",".join(offset_as)


def calc_rmsd(ref_coords, targ_coords, epsilon=1e-6):
    """Calculate the RMSD between two sets of coordinates."""
    # Center the coordinates by subtracting their means
    ref_centered = ref_coords - ref_coords.mean(axis=0)
    targ_centered = targ_coords - targ_coords.mean(axis=0)

    # Compute covariance matrix and perform SVD
    cov_matrix = targ_centered.T @ ref_centered
    U, _, Vt = np.linalg.svd(cov_matrix)

    # Adjust rotation matrix for correct handedness
    correction = np.ones((3, 3))
    correction[:, -1] = np.sign(np.linalg.det(U) * np.linalg.det(Vt))

    # Calculate rotation matrix and apply to target coordinates
    rot_matrix = (correction * U) @ Vt
    rot_targ_coords = targ_centered @ rot_matrix

    # Compute RMSD
    rmsd = np.sqrt(
        np.sum((rot_targ_coords - ref_centered) ** 2) / rot_targ_coords.shape[0]
        + epsilon
    )
    return rmsd, rot_matrix


def calc_motif_rmsd(input_path, af_path, input_indices, af_indices):
    """
    Calculate the RMSD between motifs of two PDB structures.

    params: input_path (str): Path to the input PDB file.
    params: af_path (str): Path to the AlphaFold PDB file.
    params: input_indices (list): List of indices for the motif in the input PDB.
    params: af_indices (list): List of indices for the motif in the AlphaFold PDB.

    returns: float: The RMSD value for the motif.
    returns: np.ndarray: The rotation matrix used for aligning the motifs.
    """
    # Parse PDB files to extract coordinates
    input_pdb = parse_pdb(input_path)
    af_pdb = parse_pdb(af_path)

    # Expand motif indices
    input_idx = expand(input_indices)
    af_idx = expand(af_indices)

    # Map PDB indices to internal indices
    idx_map_input = dict(zip(input_pdb["pdb_idx"], range(len(input_pdb["pdb_idx"]))))
    idx_map_af = dict(zip(af_pdb["pdb_idx"], range(len(af_pdb["pdb_idx"]))))
    input_motif_idx = [idx_map_input[idx] for idx in input_idx]
    af_motif_idx = [idx_map_af[idx] for idx in af_idx]

    # Get motif coordinates
    input_motif_xyz = input_pdb["xyz"][input_motif_idx, :3].reshape(-1, 3)
    af_motif_xyz = af_pdb["xyz"][af_motif_idx, :3].reshape(-1, 3)

    # Calculate RMSD and rotation matrix
    rmsd, rot_matrix = calc_rmsd(input_motif_xyz, af_motif_xyz)

    return rmsd, rot_matrix


def get_info(contig):
    F = []
    free_chain = False
    fixed_chain = False
    sub_contigs = [x.split("-") for x in contig.split("/")]
    for n, (a, b) in enumerate(sub_contigs):
        if a[0].isalpha():
            L = int(b) - int(a[1:]) + 1
            F += [1] * L
            fixed_chain = True
        else:
            L = int(b)
            F += [0] * L
            free_chain = True
    return F, [fixed_chain, free_chain]


def run_mpnn_sampling(af_model, mpnn_model, num_samples, batch_size, temperature):
    """Run MPNN sampling to generate protein sequences."""
    mpnn_model.get_af_inputs(af_model)
    mpnn_out = mpnn_model.sample(
        num=num_samples // batch_size,
        batch=batch_size,
        temperature=temperature,
    )
    for k in af_terms:
        mpnn_out[k] = []
    for term in other_terms:
        if term not in mpnn_out:
            mpnn_out[term] = []
    return mpnn_out


def run_af2_predictions(
    mpnn_out, af_model, num_seqs, num_recycles, af_terms, output_path, pdb_name, pdb
):
    """Run AF2 predictions on MPNN sampled sequences."""
    for n in range(num_seqs):
        seq = mpnn_out["seq"][n][-af_model._len :]
        logger.info(f"Running AF2 predictions for sequence {n+1}/{num_seqs}...")
        af_model.predict(
            seq=seq, num_recycles=num_recycles, num_models=1, verbose=False
        )

        for t in af_terms:
            mpnn_out[t].append(af_model.aux["log"][t])
        if "i_pae" in mpnn_out:
            mpnn_out["i_pae"][-1] *= 31
        if "pae" in mpnn_out:
            mpnn_out["pae"][-1] *= 31

        current_model_path = f"{output_path}/{pdb_name}_{n}.pdb"
        mpnn_out["model_path"].append(current_model_path)
        mpnn_out["input_pdb"].append(pdb)

        if input_pdb and input_pdb_as:
            motif_rmsd = calc_motif_rmsd(
                input_pdb,
                current_model_path,
                input_pdb_as,
                current_pdb_as,
            )
        mpnn_out["motif_rmsd"].append(motif_rmsd[0])

        if not save_best_only or (
            mpnn_out["plddt"][n] > 0.7 and mpnn_out["rmsd"][n] < 3
        ):
            af_model.save_current_pdb(current_model_path)

        af_model._save_results(save_best=save_best_only, verbose=False)
        af_model._k += 1
    return mpnn_out


### Checks ###
if num_seqs < mpnn_batch:
    mpnn_batch = num_seqs
    logger.warning(
        f"num_seqs must be greater than or equal to mpnn_batch. Setting mpnn_batch to {mpnn_batch}"
    )
elif num_seqs % mpnn_batch != 0:
    mpnn_batch = math.gcd(num_seqs, mpnn_batch)
    logger.warning(
        f"num_seqs must be divisible by mpnn_batch. Setting mpnn_batch to {mpnn_batch}"
    )

if rm_aa == "":
    rm_aa = None

### Parse contigs
contigs = []
for contig_str in rf_contigs.replace(" ", ":").replace(",", ":").split(":"):
    if len(contig_str) > 0:
        contig = []
        for x in contig_str.split("/"):
            if x != "0":
                contig.append(x)
        contigs.append("/".join(contig))

alphabet_list = list(ascii_uppercase + ascii_lowercase)
chains = alphabet_list[: len(contigs)]
info = [get_info(x) for x in contigs]
fixed_pos = []
fixed_chains = []
free_chains = []
both_chains = []
for pos, (fixed_chain, free_chain) in info:
    fixed_pos += pos
    fixed_chains += [fixed_chain and not free_chain]
    free_chains += [free_chain and not fixed_chain]
    both_chains += [fixed_chain and free_chain]
rm_template = np.array(fixed_pos) == 0

### Get pdb as
if input_pdb and input_pdb_as:
    current_pdb_as = get_designed_motif_offset(input_pdb_as, rf_contigs)

### Initializations ###
logger.info("Initializing symmetrical protein design...")
pdb_basename = pdb.split("/")[-1].split(".pdb")[0]

mpnn_terms = ["score", "seq"]
af_terms = ["plddt", "i_ptm", "i_pae", "rmsd"]
other_terms = ["model_path", "input_pdb"]
if input_pdb_as and input_pdb:
    other_terms.append("motif_rmsd")
labels = ["score"] + af_terms + other_terms + ["seq"]

af_model = mk_af_model(
    protocol="fixbb",
    use_templates=True,
    initial_guess=initial_guess,
    best_metric="rmsd",
    use_multimer=use_multimer,
    data_dir="/home/tsatler/projects/AFdesign_playground",
    model_names=["model_1_multimer_v3" if use_multimer else "model_1_ptm"],
)

mpnn_model = mk_mpnn_model(weights="soluble" if use_soluble else "original")


### Run MPNN sampling and AF2 predictions
logger.info("Running MPNN sampling and AF2 predictions...")
af_model.prep_inputs(
    pdb,
    chain=",".join(chains),
    rm_template=rm_template,
    rm_template_seq=rm_template,
    copies=1,
    homooligomer=False,
    rm_aa=rm_aa,
)
mpnn_out = run_mpnn_sampling(af_model, mpnn_model, num_seqs, mpnn_batch, sampling_temp)

logger.info("Running AF2 predictions...")
af2_out = run_af2_predictions(
    mpnn_out,
    af_model,
    num_seqs,
    num_recycles,
    af_terms,
    output_folder,
    pdb_basename,
    pdb,
)

# Generate model paths for all sequences
model_paths = [f"{output_folder}/{pdb_basename}_{n}.pdb" for n in range(num_seqs)]
all_labels = mpnn_terms + af_terms + other_terms
data = [[af2_out[label][n] for label in all_labels] for n in range(num_seqs)]
all_labels[0] = "mpnn"


### Save data to CSV ###
logger.info("Saving results to CSV...")
df = pd.DataFrame(data, columns=all_labels)

# Write df to CSV
output_path_all = (
    f"{results_dataframe}/af2_results_all.csv"
    if results_dataframe
    else f"{output_folder}/af2_results_all.csv"
)
write_df_to_csv(df, output_path_all)

# Filter df for best results and write to CSV
df_best = df[(df["rmsd"] < 3) & (df["plddt"] > 0.7)]
output_path_best = (
    f"{results_dataframe}/af2_best.csv"
    if results_dataframe
    else f"{output_folder}/af2_best.csv"
)
write_df_to_csv(df_best, output_path_best)
