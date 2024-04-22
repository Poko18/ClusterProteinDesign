import argparse
import fcntl
import os
import sys
import time

import numpy as np
import pandas as pd
from Bio.PDB import *
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

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
    rosetta.core.scoring.calc_per_res_hydrophobic_sasa(pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4)  # The last arguement is the probe radius

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
            print(f"Skipping residue {resname} because it is not a standard amino acid.")
            continue

    # Create a ProteinAnalysis object from the sequence
    protein_analysis = ProteinAnalysis(sequence)

    # Calculate the charge of the protein at a specific pH
    charge = protein_analysis.charge_at_pH(ph)

    return charge


def calculate_sap_score(pose, chain="B"):

    # Select only chain B using a SwitchChainOrder mover
    select_chain = XmlObjects.static_get_mover(f'<SwitchChainOrder name="so" chain_order="{chain}"/>')
    chain = pose.clone()
    select_chain.apply(chain)

    # Calculate the SAP score for chain B
    sap_score_metric = XmlObjects.static_get_simple_metric('<SapScoreMetric name="sap_metric"/>')
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
    Rg = np.sqrt(np.mean(np.sum((ca_coords - np.mean(ca_coords, axis=0)) ** 2, axis=-1)) + 1e-8)

    return Rg


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
ddg = calculate_ddg(rpose, partners=f"{binder_chain}_{''.join(target_chain.split(','))}", relax=False)
ddg_score = get_ddg(rpose, relax=False)
ddg_dsasa_100 = (ddg / interface_sasa) * 100
ddgscore_dsasa_100 = (ddg_score / interface_sasa) * 100
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
    ]

    # Save the DataFrame back to the CSV file
    df.to_csv(output_df, index=False)

    fcntl.flock(f, fcntl.LOCK_UN)  # release the lock
