import argparse
import os
import re
import shutil
from collections import Counter, OrderedDict
from string import ascii_lowercase, ascii_uppercase

import numpy as np
from mmseq_api import get_hash, get_msa, get_unique_sequences, parse_fasta, run_mmseqs2

parser = argparse.ArgumentParser(description="Argument Parser")

parser.add_argument("input_file", type=str, help="input pdb or fasta file")
parser.add_argument("output_folder", type=str, help="output folder path")
parser.add_argument("--prefix", type=str, default="sample_msa", help="out pdb prefix")

parser.add_argument("--msa_concat_mode", type=str, default="diag", help="MSA concatenation mode")
parser.add_argument(
    "--msa_mode",
    type=str,
    default="mmseqs2",
    help="msa mode: mmseqs2 or single_sequence",
)
parser.add_argument(
    "--pair_mode",
    type=str,
    default="unpaired_paired",
    help="msa pair mode: unpaired_paired, paired or unpaired",
)
parser.add_argument("--random_seed", type=int, default=0, help="Random seed")
parser.add_argument("--collapse_identical", action="store_true", help="collapse identical sequences")
parser.add_argument("--max_msa", type=int, default=256, help="Maximum number of MSA sequences")

args = parser.parse_args()

input_file = args.input_file
output_folder = args.output_folder
prefix = args.prefix

msa_concat_mode = args.msa_concat_mode
msa_mode = args.msa_mode
pair_mode = args.pair_mode

random_seed = args.random_seed
collapse_identical = args.collapse_identical
max_msa = args.max_msa
max_extra_msa = max_msa * 8
copies = 1

### Prepare inputs ###

os.makedirs(output_folder, exist_ok=True)


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


# Generate the sample.fasta file from the input PDB or copy the input FASTA to sample.fasta.
fasta_output_path = os.path.join(output_folder, "sample.fasta")

if input_file.endswith(".fasta") or input_file.endswith(".fa"):
    # If the input file is already in FASTA format, simply copy it to the sample.fasta file.
    shutil.copyfile(input_file, fasta_output_path)
else:
    # If the input file is a PDB file, generate the sample.fasta file from it.
    generate_fasta_sequence_from_pdb(input_file, output_folder)
    created_fasta_path = f"{output_folder}/{os.path.basename(input_file).split('.')[0]}.fasta"
    shutil.copyfile(created_fasta_path, fasta_output_path)


# Read the content of sample.fasta to get the sequences for MSA.
with open(fasta_output_path, "r") as file:
    fasta_string = file.read()

if fasta_string.strip():  # remove any leading or trailing whitespace
    sequences, descriptions = parse_fasta(fasta_string)
else:
    print("The FASTA file is empty.")

# Generate MSA and predict all sequences
for sequence, description in zip(sequences, descriptions):
    description = description.split("|")[0]

    sequence = re.sub("[^A-Z:]", "", sequence.replace("/", ":").upper())
    sequence = re.sub(":+", ":", sequence)
    sequence = re.sub("^[:]+", "", sequence)
    sequence = re.sub("[:]+$", "", sequence)

    sequences = sequence.replace(":", "/").split("/")

    # u_sequences = sequences
    # sequences = sum(sequences, [])
    sequence = "/".join(sequences)
    print(f"sequences: {sequence}")

    # Perform MSA for each sequence
    print("Starting MSA...")
    # get_msa(sequence, output_folder, mode=pair_mode, max_msa=max_extra_msa, msa_name=f"{description}")
    if msa_mode == "mmseqs2":
        get_msa(
            sequences,
            output_folder,
            mode=pair_mode,
            max_msa=max_extra_msa,
            msa_name=f"{description}",
        )
        # os.rename(f"{output_folder}/msa.a3m", f"{output_folder}/{description}.a3m")

    elif msa_mode == "single_sequence":
        u_sequence = "/".join(sequences)
        with open(f"{output_folder}/{description}.a3m", "w") as a3m:
            a3m.write(f">{description}\n{u_sequence}\n")

a3m_output_path = os.path.join(output_folder, "sample.a3m")
created_a3m_path = f"{output_folder}/{os.path.basename(input_file).split('.')[0]}.a3m"
shutil.copyfile(created_a3m_path, a3m_output_path)
