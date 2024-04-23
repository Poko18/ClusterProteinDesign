import argparse
import os
import re


def generate_fasta_sequence_from_pdb(pdb_file, output_file):
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

    with open(output_file, "w") as fp:
        fp.write(fasta_sequence)

    print(f"FASTA sequence saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate FASTA sequence from PDB file")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("output_file", help="Output file for the FASTA sequence")
    args = parser.parse_args()

    generate_fasta_sequence_from_pdb(args.pdb_file, args.output_file)
