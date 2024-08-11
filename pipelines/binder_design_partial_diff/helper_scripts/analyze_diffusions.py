import argparse
import fcntl
import os
import time

import numpy as np
import pandas as pd

# ARGS PARSER
def parse_args():
    parser = argparse.ArgumentParser(description="Analyze PDB files and save results to CSV.")
    parser.add_argument('--input_txt', type=str, help='Text file with paths to PDB files.')
    parser.add_argument('--hotspots', type=str, help='Hotspots to analyze.')
    parser.add_argument('--target', type=str, help='Target PDB file to analyze.')
    parser.add_argument('--binder_chain', type=str, help='Chain identifier for the binder.')
    parser.add_argument('--target_chain', type=str, help='Chain identifier for the target.')
    parser.add_argument('--output_csv', type=str, help='Output CSV file to save results.')
    return parser.parse_args()


# dictionaries
num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    'UNK','MAS',
    ]

one_letter = ["A", "R", "N", "D", "C", \
             "Q", "E", "G", "H", "I", \
             "L", "K", "M", "F", "P", \
             "S", "T", "W", "Y", "V", "?", "-"]

aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD "," HE ","1HH1","2HH1","1HH2","2HH2"), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD2","2HD2",  None,  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE2","2HE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ",  None,  None,  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  ","1HA ","2HA ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HD2"," HE1"," HE2",  None,  None,  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG2","2HG2","3HG2","1HG1","2HG1","1HD1","2HD1","3HD1",  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2",  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ","1HE ","2HE ","1HZ ","2HZ ","3HZ "), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE ","2HE ","3HE ",  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HD2"," HE1"," HE2"," HZ ",  None,  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ",  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HG "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HG1"," HA "," HB ","1HG2","2HG2","3HG2",  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HZ2"," HH2"," HZ3"," HE3",  None,  None,  None), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HE2"," HD2"," HH ",  None,  None,  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2",  None,  None,  None,  None), # val
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # unk
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # mask
]

aa2num= {x:i for i,x in enumerate(num2aa)}

aa_321 = {a:b for a,b in zip(num2aa,one_letter)}
aa_123 = {val:key for key,val in aa_321.items()}

# Sokrypton's code
MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

restype_1to3 = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}

restype_3to1 = {v: k for k, v in restype_1to3.items()}

# RFdif
def parse_pdb(filename, **kwargs):
  """extract xyz coords for all heavy atoms"""
  lines = open(filename, "r").readlines()
  return parse_pdb_lines(lines, **kwargs)


def get_pdb_seq(lines) -> str:
    """get the sequence from a pdb file"""
    seq = ""
    for line in lines:
        if line[:4] == "ATOM":
            if line[13:15] == "CA":
                seq += line[17:20]
    return seq

def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res, pdb_idx = [],[]
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
        if (chain,resNo) in pdb_idx:
          idx = pdb_idx.index((chain, resNo))
          # for i_atm, tgtatm in enumerate(util.aa2long[util.aa2num[aa]]):
          for i_atm, tgtatm in enumerate(
              aa2long[aa2num[aa]][:14]
          ):  # Nate's proposed change
              if (
                  tgtatm is not None and tgtatm.strip() == atom.strip()
              ):  # ignore whitespace
                  xyz[idx, i_atm, :] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
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

def get_binder_ca_xyz(pdb_file, binder_chain="A"):
    pdb = parse_pdb(pdb_file)
    ca_xyz = pdb["xyz"][:, 1, :]
    chain_mask = np.array([chain_id == binder_chain for chain_id, res_num in pdb["pdb_idx"]])
    binder_xyz = ca_xyz[chain_mask]
    return binder_xyz

def get_hotspots_ca_xyz(pdb_file, hotspots, target_chain="B"):
    pdb = parse_pdb(pdb_file)
    ca_xyz = pdb["xyz"][:, 1, :]
    hotspots_mask = np.array([chain_id == target_chain and res_num in hotspots for chain_id, res_num in pdb["pdb_idx"]])
    hotspots_xyz = ca_xyz[hotspots_mask]
    return hotspots_xyz

def extract_hotspot_numbers(hotspots):
    return [value[1:] for value in hotspots.split('=')[1].strip('[]').split(',')]
def extract_target_chain(hotspots):
    return hotspots.split('=')[1][1]

def get_target_offset(pdb_file, target_chain="B"):
    pdb = parse_pdb(pdb_file)
    for chain_id, res_num in pdb["pdb_idx"]:
        if chain_id == target_chain:
            return res_num-1
    return None

# Rf diffusion have binder chain A and target chain B
# Predicted sequences have binder chain B and target chain A
def analyze_binder_diffusions(pdb_file, hotspots, hotspot_offset=0, binder_chain="A", target_chain="B"):
    """Analyze the binder diffusions"""
    def max_distance_between_atoms(coords):
        dists = np.linalg.norm(coords[:, np.newaxis] - coords[np.newaxis, :], axis=-1)
        return np.max(dists)
    
    def min_distances_between_hotspots_and_binder(hotspots_xyz, binder_xyz):
        hotspots_xyz = np.expand_dims(hotspots_xyz, axis=1)
        binder_xyz = np.expand_dims(binder_xyz, axis=0)
        squared_diffs = (hotspots_xyz - binder_xyz) ** 2
        squared_distances = np.sum(squared_diffs, axis=2)
        distances = np.sqrt(squared_distances)
        return np.min(distances, axis=1)

    binder_xyz = get_binder_ca_xyz(pdb_file, binder_chain)
    hotspots = np.array(hotspots, dtype=int)
    if binder_chain == "A":
        # print("Binder chain is A... RF diffusion generated")
        corrected_hotspots = hotspots + len(binder_xyz) - hotspot_offset
    elif binder_chain == "B":
        # print("Binder chain is B... AF2 generated")
        corrected_hotspots = hotspots - hotspot_offset
    hotspots_xyz = get_hotspots_ca_xyz(pdb_file, corrected_hotspots, target_chain)

    # Get an average distance from hotspot ca atoms to the closest atom of the binder
    distances = min_distances_between_hotspots_and_binder(hotspots_xyz, binder_xyz)

    # Get biggest max distance between Ca atoms of the binder and binder
    max_binders = max_distance_between_atoms(binder_xyz)

    # Calculate radius of gyration of binder
    center = np.mean(binder_xyz, axis=0)
    squared_distances = np.sum((binder_xyz - center)**2, axis=1)
    radius_of_gyration = np.sqrt(np.mean(squared_distances) + 1e-8)

    return np.mean(distances), max_binders, radius_of_gyration

# hotspots='ppi.hotspot_res=[A69,A70,A39,A146,A156,A157,A160]'


# rf_distance_dict = {}
# for rf_pdb in rf_pdbs:
#     distance = analyze_binder_diffusions(rf_pdb, extract_hotspot_numbers(hotspots), hotspot_offset=int(get_target_offset(target)), binder_chain="A", target_chain="B")
#     rf_distance_dict[rf_pdb] = distance

# rf_dist = []
# af_dist = []
# for i, row in all_predictions.iterrows():
#     rf_pdb = row["input_pdb"]
#     af_pdb = row["model_path"]
#     distance = analyze_binder_diffusions(af_pdb, extract_hotspot_numbers(hotspots), hotspot_offset=int(get_target_offset(target)), binder_chain="B", target_chain="A")
#     rf_dist.append(rf_distance_dict[rf_pdb])
#     af_dist.append(distance)

# all_predictions["rf_distance"] = rf_dist
# all_predictions["af_distance"] = af_dist
# all_predictions

def main():
    args = parse_args()
    input_txt = args.input_txt
    hotspots = args.hotspots
    target = args.target
    binder_chain = args.binder_chain
    target_chain = args.target_chain
    output_csv = args.output_csv
    
    # Read PDB file paths from the input text file
    with open(input_txt, 'r') as file:
        pdb_files = [line.strip() for line in file]

    # Initialize lists to collect results
    results = []

    for pdb_file in pdb_files:
        # Extract hotspots from a predefined source (can be customized as needed)
        hotspot_numbers = extract_hotspot_numbers(hotspots)
        hotspot_chain = extract_target_chain(hotspots)
        hotspot_offset = int(get_target_offset(target, target_chain=hotspot_chain))
        
        mean_distance, binder_max_distance, rg = analyze_binder_diffusions(
            pdb_file,
            hotspot_numbers,
            hotspot_offset=hotspot_offset,
            binder_chain=binder_chain,
            target_chain=target_chain
        )
        results.append({
            'model_path': pdb_file,
            'mean_hotspot_distance': mean_distance,
            'binder_max_distance': binder_max_distance,
            'rg': rg
        })

    results_df = pd.DataFrame(results)
    file_exists = os.path.isfile(output_csv)

    with open(output_csv, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        time.sleep(np.random.uniform(0, 0.05))
        results_df.to_csv(f, header=not file_exists, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)

if __name__ == "__main__":
    main()
