import argparse
import fcntl
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Script to calculate and save scores in a CSV file.")
parser.add_argument("ref_pdb", type=str, help="Path to the reference PDB file")
parser.add_argument("sample_pdb", type=str, help="Path to the sample PDB file")
parser.add_argument("sample_json", type=str, help="Path to the sample JSON file")
parser.add_argument("output_df", type=str, help="Output folder")

parser.add_argument("--out_name", type=str, default="af2_scores", help="Path to the sample JSON file")


args = parser.parse_args()

ref_pdb = args.ref_pdb
sample_pdb = args.sample_pdb
sample_json = args.sample_json
output_df = args.output_df
out_name = args.out_name

### Define functions


def generate_fasta_sequence_from_pdb(pdb_file, output_folder=None):
    ca_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
    aa3to1 = {
        "ALA": "A",
        "VAL": "V",
        "PHE": "F",
        "PRO": "P",
        "MET": "M",
        "ILE": "I",
        "LEU": "L",
        "ASP": "D",
        "GLU": "E",
        "LYS": "K",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "TYR": "Y",
        "HIS": "H",
        "CYS": "C",
        "ASN": "N",
        "GLN": "Q",
        "TRP": "W",
        "GLY": "G",
        "MSE": "M",
    }
    filename = os.path.basename(pdb_file).split(".")[0]
    chain_dict = dict()
    chain_list = []

    with open(pdb_file, "r") as fp:
        for line in fp:
            if line.startswith("ENDMDL"):
                break
            match_list = ca_pattern.findall(line)
            if match_list:
                resn = match_list[0][0] + match_list[0][2]
                chain = match_list[0][1] + match_list[0][3]
                if chain in chain_dict:
                    chain_dict[chain] += aa3to1[resn]
                else:
                    chain_dict[chain] = aa3to1[resn]
                    chain_list.append(chain)

    fasta_sequence = f">{filename}\n"
    for i, chain in enumerate(chain_list):
        fasta_sequence += chain_dict[chain]
        if i < len(chain_list) - 1:
            fasta_sequence += ":"

    if output_folder:
        output_file = os.path.join(output_folder, f"{filename}.fasta")
        with open(output_file, "w") as fp:
            fp.write(fasta_sequence)

    return chain_dict


def calculate_scores(scores_path, binder_len=None, is_binder_second=False):
    scores = json.loads(Path(scores_path).read_text())

    plddt = np.mean(scores["plddt"])
    pae = np.array(scores["pae"])

    if is_binder_second:
        pae_binder = np.mean(pae[binder_len:, binder_len:]) if binder_len else None
        pae_target = np.mean(pae[:binder_len, :binder_len]) if binder_len else None
        plddt_binder = np.mean(scores["plddt"][binder_len:]) if binder_len else None
        plddt_target = np.mean(scores["plddt"][:binder_len]) if binder_len else None
        pae_int1 = np.mean(pae[binder_len:, :binder_len]) if binder_len else None
        pae_int2 = np.mean(pae[:binder_len, binder_len:]) if binder_len else None
    else:
        pae_binder = np.mean(pae[:binder_len, :binder_len]) if binder_len else None
        pae_target = np.mean(pae[binder_len:, binder_len:]) if binder_len else None
        plddt_binder = np.mean(scores["plddt"][:binder_len]) if binder_len else None
        plddt_target = np.mean(scores["plddt"][binder_len:]) if binder_len else None
        pae_int1 = np.mean(pae[:binder_len, binder_len:]) if binder_len else None
        pae_int2 = np.mean(pae[binder_len:, :binder_len]) if binder_len else None

    pae_int_tot = (pae_int1 + pae_int2) / 2 if binder_len else None

    results = {"plddt": plddt, "pae": np.mean(pae)}
    if binder_len:
        results.update(
            {
                "binder_plddt": plddt_binder,
                "target_plddt": plddt_target,
                "pae_binder": pae_binder,
                "pae_target": pae_target,
                "pae_int_tot": pae_int_tot,
            }
        )

    return results


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


def get_plddt_values(input_pdb):
    plddt_values = []

    with open(input_pdb, "r") as f:
        for line in f:
            if line[0:4] == "ATOM":
                plddt_value = float(line[60:66])
                plddt_values.append(plddt_value)

    return plddt_values


### Run script

# Get plddt from pdb
plddt_values = get_plddt_values(sample_pdb)
mean_plddt_pdb = np.mean(plddt_values)

scores = calculate_scores(sample_json)
rmsd = align_structures(ref_pdb, sample_pdb, save_aligned=True)
scores["rmsd"] = rmsd
scores["name"] = sample_pdb
scores["model_id"] = sample_pdb.split("/")[-1].split(".")[0]

# Save scores in a DataFrame
scores_df = pd.DataFrame(scores, index=[0])
scores_filepath = f"{output_df}/{out_name}.csv"

if os.path.isfile(scores_filepath):
    # Lock the file for exclusive access
    with open(scores_filepath, "r+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)

        existing_scores_df = pd.read_csv(scores_filepath)
        scores_df = pd.concat([existing_scores_df, scores_df], ignore_index=True)
        scores_df = scores_df.drop_duplicates("name")

        # Save the combined DataFrame back to the CSV file
        scores_df.to_csv(scores_filepath, index=False)

        fcntl.flock(f, fcntl.LOCK_UN)  # Release the lock
else:
    # If the file doesn't exist, just save the DataFrame without locking
    scores_df.to_csv(scores_filepath, index=False)
