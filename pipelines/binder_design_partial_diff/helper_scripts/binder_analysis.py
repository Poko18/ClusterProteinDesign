import argparse
import pandas as pd
import argparse
import os
import sys
import fcntl
import time
import numpy as np
from pyrosetta.rosetta import *
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.protocols import *
from Bio.PDB import *
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.SeqUtils.ProtParam import ProteinAnalysis

os.getcwd()

from pyrosetta import *

init()
from pyrosetta.rosetta import *

init("-beta_nov16 -holes:dalphaball")
# init("-beta_nov16 -holes:dalphaball")


parser = argparse.ArgumentParser(description="Binder analysis")

parser.add_argument("input_pdb", type=str, help="Path to input PDB files")
parser.add_argument("target_chain", type=str, help="name of target chain (C,D..)")
parser.add_argument("binder_chain", type=str, help="name of binder chain (C,D..)")
# parser.add_argument('partners', default="B_C", type=str, help='Partners to unbind, separated by "_". Example: "A_B"')
parser.add_argument("output_df", type=str, help="Output csv file for adding metrics")
parser.add_argument("xml_file", type=str, help="Path to rosetta XML file")

args = parser.parse_args()
input_pdb = args.input_pdb
target_chain = args.target_chain
binder_chain = args.binder_chain
# partners = args.partners
output_df = args.output_df

xml = args.xml_file
objs = protocols.rosetta_scripts.XmlObjects.create_from_file(xml)

# Check if output_df is a file
if not os.path.isfile(output_df):
    print(f"{output_df} is not a file. Exiting the program.")
    sys.exit()

# Check if input_pdb is in dataframe
with open(output_df, "r") as f:
    fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
    time.sleep(np.random.uniform(0, 0.05))

    df = pd.read_csv(output_df)

    # Check if input_pdb is in df['model_path']
    if input_pdb not in df["model_path"].values:
        raise ValueError(f"Input pdb '{input_pdb}' not found in the DataFrame.")

    fcntl.flock(f, fcntl.LOCK_UN)  # release the lock

print(f"Input PDB is in metrics dataframe!")

num2aa = [
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
    "UNK",
    "MAS",
]

one_letter = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
    "?",
    "-",
]

aa2long = [
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "3HB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # ala
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD ",
        " NE ",
        " CZ ",
        " NH1",
        " NH2",
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "1HG ",
        "2HG ",
        "1HD ",
        "2HD ",
        " HE ",
        "1HH1",
        "2HH1",
        "1HH2",
        "2HH2",
    ),  # arg
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " OD1",
        " ND2",
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "1HD2",
        "2HD2",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # asn
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " OD1",
        " OD2",
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # asp
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " SG ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        " HG ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # cys
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD ",
        " OE1",
        " NE2",
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "1HG ",
        "2HG ",
        "1HE2",
        "2HE2",
        None,
        None,
        None,
        None,
        None,
    ),  # gln
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD ",
        " OE1",
        " OE2",
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "1HG ",
        "2HG ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # glu
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        "1HA ",
        "2HA ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # gly
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " ND1",
        " CD2",
        " CE1",
        " NE2",
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        " HD2",
        " HE1",
        " HE2",
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # his
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG1",
        " CG2",
        " CD1",
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        " HB ",
        "1HG2",
        "2HG2",
        "3HG2",
        "1HG1",
        "2HG1",
        "1HD1",
        "2HD1",
        "3HD1",
        None,
        None,
    ),  # ile
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD1",
        " CD2",
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        " HG ",
        "1HD1",
        "2HD1",
        "3HD1",
        "1HD2",
        "2HD2",
        "3HD2",
        None,
        None,
    ),  # leu
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD ",
        " CE ",
        " NZ ",
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "1HG ",
        "2HG ",
        "1HD ",
        "2HD ",
        "1HE ",
        "2HE ",
        "1HZ ",
        "2HZ ",
        "3HZ ",
    ),  # lys
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " SD ",
        " CE ",
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "1HG ",
        "2HG ",
        "1HE ",
        "2HE ",
        "3HE ",
        None,
        None,
        None,
        None,
    ),  # met
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD1",
        " CD2",
        " CE1",
        " CE2",
        " CZ ",
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        " HD1",
        " HD2",
        " HE1",
        " HE2",
        " HZ ",
        None,
        None,
        None,
        None,
    ),  # phe
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " HA ",
        "1HB ",
        "2HB ",
        "1HG ",
        "2HG ",
        "1HD ",
        "2HD ",
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # pro
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " OG ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HG ",
        " HA ",
        "1HB ",
        "2HB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # ser
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " OG1",
        " CG2",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HG1",
        " HA ",
        " HB ",
        "1HG2",
        "2HG2",
        "3HG2",
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # thr
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD1",
        " CD2",
        " NE1",
        " CE2",
        " CE3",
        " CZ2",
        " CZ3",
        " CH2",
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        " HD1",
        " HE1",
        " HZ2",
        " HH2",
        " HZ3",
        " HE3",
        None,
        None,
        None,
    ),  # trp
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG ",
        " CD1",
        " CD2",
        " CE1",
        " CE2",
        " CZ ",
        " OH ",
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        " HD1",
        " HE1",
        " HE2",
        " HD2",
        " HH ",
        None,
        None,
        None,
        None,
    ),  # tyr
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        " CG1",
        " CG2",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        " HB ",
        "1HG1",
        "2HG1",
        "3HG1",
        "1HG2",
        "2HG2",
        "3HG2",
        None,
        None,
        None,
        None,
    ),  # val
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "3HB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # unk
    (
        " N  ",
        " CA ",
        " C  ",
        " O  ",
        " CB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        " H  ",
        " HA ",
        "1HB ",
        "2HB ",
        "3HB ",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ),  # mask
]

aa2num = {x: i for i, x in enumerate(num2aa)}

aa_321 = {a: b for a, b in zip(num2aa, one_letter)}
aa_123 = {val: key for key, val in aa_321.items()}

# Sokrypton's code
MODRES = {
    "MSE": "MET",
    "MLY": "LYS",
    "FME": "MET",
    "HYP": "PRO",
    "TPO": "THR",
    "CSO": "CYS",
    "SEP": "SER",
    "M3L": "LYS",
    "HSK": "HIS",
    "SAC": "SER",
    "PCA": "GLU",
    "DAL": "ALA",
    "CME": "CYS",
    "CSD": "CYS",
    "OCS": "CYS",
    "DPR": "PRO",
    "B3K": "LYS",
    "ALY": "LYS",
    "YCM": "CYS",
    "MLZ": "LYS",
    "4BF": "TYR",
    "KCX": "LYS",
    "B3E": "GLU",
    "B3D": "ASP",
    "HZP": "PRO",
    "CSX": "CYS",
    "BAL": "ALA",
    "HIC": "HIS",
    "DBZ": "ALA",
    "DCY": "CYS",
    "DVA": "VAL",
    "NLE": "LEU",
    "SMC": "CYS",
    "AGM": "ARG",
    "B3A": "ALA",
    "DAS": "ASP",
    "DLY": "LYS",
    "DSN": "SER",
    "DTH": "THR",
    "GL3": "GLY",
    "HY3": "PRO",
    "LLP": "LYS",
    "MGN": "GLN",
    "MHS": "HIS",
    "TRQ": "TRP",
    "B3Y": "TYR",
    "PHI": "PHE",
    "PTR": "TYR",
    "TYS": "TYR",
    "IAS": "ASP",
    "GPL": "LYS",
    "KYN": "TRP",
    "CSD": "CYS",
    "SEC": "CYS",
}

restype_1to3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

restype_3to1 = {v: k for k, v in restype_1to3.items()}


def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res, pdb_idx = [], []
    for l in lines:
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            res.append((l[22:26], l[17:20]))
            # chain letter, res num
            pdb_idx.append((l[21:22].strip(), int(l[22:26].strip())))
    seq = [aa2num[r[1]] if r[1] in aa2num.keys() else 20 for r in res]

    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(res), 14, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = (
            l[21:22],
            int(l[22:26]),
            " " + l[12:16].strip().ljust(3),
            l[17:20],
        )
        if (chain, resNo) in pdb_idx:
            idx = pdb_idx.index((chain, resNo))
            # for i_atm, tgtatm in enumerate(util.aa2long[util.aa2num[aa]]):
            for i_atm, tgtatm in enumerate(
                aa2long[aa2num[aa]][:14]
            ):  # Nate's proposed change
                if (
                    tgtatm is not None and tgtatm.strip() == atom.strip()
                ):  # ignore whitespace
                    xyz[idx, i_atm, :] = [
                        float(l[30:38]),
                        float(l[38:46]),
                        float(l[46:54]),
                    ]
                    break

    # save atom mask
    mask = np.logical_not(np.isnan(xyz[..., 0]))
    xyz[np.isnan(xyz[..., 0])] = 0.0

    # remove duplicated (chain, resi)
    new_idx = []
    i_unique = []
    for i, idx in enumerate(pdb_idx):
        if idx not in new_idx:
            new_idx.append(idx)
            i_unique.append(i)

    pdb_idx = new_idx
    xyz = xyz[i_unique]
    mask = mask[i_unique]

    seq = np.array(seq)[i_unique]

    out = {
        "xyz": xyz,  # cartesian coordinates, [Lx14]
        "mask": mask,  # mask showing which atoms are present in the PDB file, [Lx14]
        "idx": np.array(
            [i[1] for i in pdb_idx]
        ),  # residue numbers in the PDB file, [L]
        "seq": np.array(seq),  # amino acid sequence, [L]
        "pdb_idx": pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
    }

    # heteroatoms (ligands, etc)
    if parse_hetatom:
        xyz_het, info_het = [], []
        for l in lines:
            if l[:6] == "HETATM" and not (ignore_het_h and l[77] == "H"):
                info_het.append(
                    dict(
                        idx=int(l[7:11]),
                        atom_id=l[12:16],
                        atom_type=l[77],
                        name=l[16:20],
                    )
                )
                xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

        out["xyz_het"] = np.array(xyz_het)
        out["info_het"] = info_het

    return out


def parse_pdb(filename, **kwargs):
    """extract xyz coords for all heavy atoms"""
    lines = open(filename, "r").readlines()
    return parse_pdb_lines(lines, **kwargs)


### Functions ###


def fast_relax_pdb(pose, scorefxn=get_fa_scorefxn()):
    scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_rep, 0)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.constrain_relax_to_start_coords(True)
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(True)
    relax.set_movemap(mm)
    relax.apply(pose)


def unbind(pose, partners):
    STEP_SIZE = 100
    JUMP = 1  # was set to 2 but not working
    docking.setup_foldtree(pose, partners, Vector1([-1, -1, -1]))
    trans_mover = rigid.RigidBodyTransMover(pose, JUMP)
    trans_mover.step_size(STEP_SIZE)
    trans_mover.apply(pose)


def calculate_ddg(pose, partners, scorefxn=get_fa_scorefxn(), relax=True):
    # Load the PDB file as a PyRosetta Pose
    start_pose = pose_from_file(f"{input_pdb}")
    pose = start_pose.clone()

    # Relax the pose
    if relax:
        relax_pose(pose)

    # Save the relaxed structure
    relaxPose = pose.clone()

    # Calculate the bound score
    scorefxn = get_fa_scorefxn()
    scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_rep, 0)
    bound_score = scorefxn(relaxPose)
    # scorefxn.show(relaxPose)

    # Unbind chains
    unbind(relaxPose, partners)

    # Calculate the unbound score
    unbound_score = scorefxn(relaxPose)
    unboundPose = relaxPose.clone()

    # Calculate ddG
    return round((bound_score - unbound_score), 3)


def align_structures(pdb1, pdb2):
    """Take two structure and superimpose pdb1 on pdb2"""
    import Bio.PDB

    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    # Get the structures
    ref_structure = pdb_parser.get_structure("ref", pdb1)
    sample_structure = pdb_parser.get_structure("sample", pdb2)

    aligner = Bio.PDB.cealign.CEAligner()
    aligner.set_reference(ref_structure)
    aligner.align(sample_structure)

    return aligner.rms


def get_sasa(pose):
    """Calculate the total and hydrophobic sasa"""
    rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
    rsd_hydrophobic_sasa = pyrosetta.rosetta.utility.vector1_double()
    rosetta.core.scoring.calc_per_res_hydrophobic_sasa(
        pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4
    )  # The last arguement is the probe radius

    return sum(rsd_sasa), sum(rsd_hydrophobic_sasa)


def calculate_charge(chain, ph=7.4):
    """
    Calculates the charge of the protein chain at a specific pH.

    Parameters:
        chain (Bio.PDB.Chain.Chain): The protein chain to calculate the charge of.
        ph (float): The pH value to calculate the charge at. Default is 7.4.

    Returns:
        The charge of the protein at the specified pH.
    """
    # Extract the sequence of amino acid residues in the chain
    sequence = ""
    for residue in chain.get_residues():
        resname = residue.get_resname()
        if resname in aa3:
            sequence += aa1[aa3.index(resname)]
        else:
            print(
                f"Skipping residue {resname} because it is not a standard amino acid."
            )
            continue

    # Create a ProteinAnalysis object from the sequence
    protein_analysis = ProteinAnalysis(sequence)

    # Calculate the charge of the protein at a specific pH
    charge = protein_analysis.charge_at_pH(ph)

    return charge


def calculate_sap_score(pose, chain="B"):

    # Select only chain B using a SwitchChainOrder mover
    select_chain = XmlObjects.static_get_mover(
        f'<SwitchChainOrder name="so" chain_order="{chain}"/>'
    )
    chain = pose.clone()
    select_chain.apply(chain)

    # Calculate the SAP score for chain B
    sap_score_metric = XmlObjects.static_get_simple_metric(
        '<SapScoreMetric name="sap_metric"/>'
    )
    sap_score_value = sap_score_metric.calculate(chain)

    # Return the SAP score value
    # sap_score_value = sap_score_metric.get(1)
    return sap_score_value


def calculate_rg(chain):
    """
    Calculates the radius of gyration (Rg) of a protein chain using only the alpha carbons (CA).

    Parameters:
        chain (Bio.PDB.Chain.Chain): The protein chain to calculate the Rg of.

    Returns:
        The Rg of the protein chain.
    """
    # Get the coordinates of all alpha carbons in the chain
    ca_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == "CA"]
    ca_coords = np.array([atom.get_coord() for atom in ca_atoms])

    # Calculate the Rg of the protein
    Rg = np.sqrt(
        np.mean(np.sum((ca_coords - np.mean(ca_coords, axis=0)) ** 2, axis=-1)) + 1e-8
    )

    return Rg


def get_binder_ca_xyz(pdb_file, binder_chain="A"):
    pdb = parse_pdb(pdb_file)
    ca_xyz = pdb["xyz"][:, 1, :]
    chain_mask = np.array(
        [chain_id == binder_chain for chain_id, res_num in pdb["pdb_idx"]]
    )
    binder_xyz = ca_xyz[chain_mask]
    return binder_xyz


def max_distance_between_atoms(coords):
    dists = np.linalg.norm(coords[:, np.newaxis] - coords[np.newaxis, :], axis=-1)
    return np.max(dists)


def interface_terms(pdb):
    # Returns:
    # - The dG (delta G) energy value of the interface
    # - The dSASA (delta Solvent Accessible Surface Area) value of the interface
    # - The dG_dSASA_ratio (delta G to delta SASA ratio) multiplied by 100
    # - The number of delta unsatisfied hydrogen bonds at the interface
    # - The number of hydrogen bonds formed at the interface

    pose = pose_from_pdb(pdb)
    interface_analyzer = analysis.InterfaceAnalyzerMover()
    interface_analyzer.apply(pose)
    data = interface_analyzer.get_all_data()
    return (
        data.dG[1],
        data.dSASA[1],
        (data.dG_dSASA_ratio * 100),
        data.delta_unsat_hbonds,
        data.interface_hbonds,
    )


def relax_pose(pose, binder_chain="A"):
    if binder_chain == "A":
        FastRelax = objs.get_mover("FastRelax")
        FastRelax.apply(pose)
    elif binder_chain == "B":
        FastRelax = objs.get_mover("FastRelax_ChainB_bbtrue")
        FastRelax.apply(pose)
    return pose


def get_ddg(pose, relax=True):
    if relax:
        relax_pose(pose)
    ddg = objs.get_filter("ddg")
    ddg.apply(pose)
    return ddg.score(pose)


def shape_complementarity(pose):
    shape_comp = objs.get_filter("interface_sc")
    shape_comp.apply(pose)
    return shape_comp.score(pose)


def interface_buried_sasa(pose):
    dsasa = objs.get_filter("interface_buried_sasa")
    dsasa.apply(pose)
    return dsasa.score(pose)


def hydrophobic_residue_contacts(pose):
    hyd_res = objs.get_filter("hydrophobic_residue_contacts")
    hyd_res.apply(pose)
    return hyd_res.score(pose)


def get_cms(pose):
    cms = objs.get_filter("cms")
    cms.apply(pose)
    return cms.score(pose)


def get_vbuns(pose):
    vbuns = objs.get_filter("vbuns")
    vbuns.apply(pose)
    return vbuns.score(pose)


def get_sbuns(pose):
    sbuns = objs.get_filter("sbuns")
    sbuns.apply(pose)
    return sbuns.score(pose)


def interface_vbuns(pose, partners):
    bound_vbuns = get_vbuns(pose)
    UnboundPose = pose.clone()
    unbind(UnboundPose, partners)
    unbound_vbuns = get_vbuns(UnboundPose)
    return round(bound_vbuns, 3), round(unbound_vbuns, 3)


def interface_sbuns(pose, partners):
    bound_sbuns = get_sbuns(pose)
    UnboundPose = pose.clone()
    unbind(UnboundPose, partners)
    unbound_sbuns = get_sbuns(UnboundPose)
    return round(bound_sbuns, 3), round(unbound_sbuns, 3)


##################

### Metric calculations ###

pdb = input_pdb
pose = pose_from_file(input_pdb)
parser = PDBParser()
structure = parser.get_structure("protein", pdb)
chain = structure[0][binder_chain]

# TO DO: call functions in dict, based on metric_columns defined

rg = calculate_rg(chain)
charge = calculate_charge(chain, ph=7.4)
sap = calculate_sap_score(pose, binder_chain)
dG, dSASA, dG_dSASA_ratio, int_unsat_hbonds, int_hbonds = interface_terms(pdb)
hyd_con = hydrophobic_residue_contacts(pose)
shape_comp = shape_complementarity(pose)
interface_sasa = interface_buried_sasa(pose)
rpose = relax_pose(pose, binder_chain)
ddg = calculate_ddg(
    rpose, partners=f"{binder_chain}_{''.join(target_chain.split(','))}", relax=False
)
ddg_score = get_ddg(rpose, relax=False)
ddg_dsasa_100 = (ddg / interface_sasa) * 100
ddgscore_dsasa_100 = (ddg_score / interface_sasa) * 100
cms = get_cms(pose)
ddg_cms_100 = (ddg / cms) * 100
vbuns_bound, vbuns_unbound = interface_vbuns(
    rpose, partners=f"{binder_chain}_{''.join(target_chain.split(','))}"
)
vbuns_int = vbuns_bound - vbuns_unbound
sbuns_bound, sbuns_unbound = interface_sbuns(
    rpose, partners=f"{binder_chain}_{''.join(target_chain.split(','))}"
)
sbuns_int = sbuns_bound - sbuns_unbound
max_binder_distance = max_distance_between_atoms(get_binder_ca_xyz(pdb, binder_chain))
################## ddg_dsasa_100 , ddgscore_dsasa_100

metric_columns = [
    "ddg",
    "rg",
    "charge",
    "sap",
    "dG",
    "dSASA",
    "dG_dSASA_ratio",
    "int_unsat_hbonds",
    "int_hbonds",
    "hyd_contacts",
    "shape_comp",
    "ddg_score",
    "ddg_dsasa_100",
    "ddgscore_dsasa_100",
    "cms",
    "vbuns_bound",
    "vbuns_unbound",
    "vbuns_int",
    "sbuns_bound",
    "sbuns_unbound",
    "sbuns_int",
    "max_binder_distance",
]

# Write the data to the file, acquiring a lock if necessary
with open(output_df, "r+") as f:
    fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
    time.sleep(np.random.uniform(0, 0.05))

    df = pd.read_csv(output_df)

    # Add missing metric columns if they don't exist
    for column in metric_columns:
        if column not in df.columns:
            df[column] = None

    df.loc[df["model_path"] == input_pdb, metric_columns] = [
        ddg,
        rg,
        charge,
        sap,
        dG,
        interface_sasa,
        dG_dSASA_ratio,
        int_unsat_hbonds,
        int_hbonds,
        hyd_con,
        shape_comp,
        ddg_score,
        ddg_dsasa_100,
        ddgscore_dsasa_100,
        cms,
        vbuns_bound,
        vbuns_unbound,
        vbuns_int,
        sbuns_bound,
        sbuns_unbound,
        sbuns_int,
        max_binder_distance,
    ]

    # Save the DataFrame back to the CSV file
    df.to_csv(output_df, index=False)

    fcntl.flock(f, fcntl.LOCK_UN)  # release the lock
