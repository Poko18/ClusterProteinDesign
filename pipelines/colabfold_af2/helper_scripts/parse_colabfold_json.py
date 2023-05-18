import argparse
import glob
import os
import json
import numpy as np
from pathlib import Path
import pandas as pd
import shutil

def calculate_scores(scores_path, binder_len=None, is_binder_second=False):
    scores = json.loads(Path(scores_path).read_text())

    plddt = np.mean(scores['plddt'])
    pae = np.array(scores['pae'])

    if is_binder_second:
        pae_binder = np.mean(pae[binder_len:, binder_len:]) if binder_len else None
        pae_target = np.mean(pae[:binder_len, :binder_len]) if binder_len else None
        plddt_binder = np.mean(scores['plddt'][binder_len:]) if binder_len else None
        plddt_target = np.mean(scores['plddt'][:binder_len]) if binder_len else None
        pae_int1 = np.mean(pae[binder_len:, :binder_len]) if binder_len else None
        pae_int2 = np.mean(pae[:binder_len, binder_len:]) if binder_len else None
    else:
        pae_binder = np.mean(pae[:binder_len, :binder_len]) if binder_len else None
        pae_target = np.mean(pae[binder_len:, binder_len:]) if binder_len else None
        plddt_binder = np.mean(scores['plddt'][:binder_len]) if binder_len else None
        plddt_target = np.mean(scores['plddt'][binder_len:]) if binder_len else None
        pae_int1 = np.mean(pae[:binder_len, binder_len:]) if binder_len else None
        pae_int2 = np.mean(pae[binder_len:, :binder_len]) if binder_len else None

    pae_int_tot = (pae_int1 + pae_int2) / 2 if binder_len else None

    results = {'plddt': plddt, 'pae': np.mean(pae)}
    if binder_len:
        results.update({
            'binder_plddt': plddt_binder,
            'target_plddt': plddt_target,
            'pae_binder': pae_binder,
            'pae_target': pae_target,
            'pae_int_tot': pae_int_tot
        })

    return results


parser = argparse.ArgumentParser(description='Calculate scores.')
parser.add_argument('scores_path', type=str, help='Path to score JSON file or folder with score file')
parser.add_argument('--binder_len', type=int, help='Length of binder protein')
parser.add_argument('--is_binder_second', action='store_true', help='Whether binder protein is the second protein')
parser.add_argument('--output', type=str, help='Path to the output file')

args = parser.parse_args()

if os.path.isdir(args.scores_path):
    score_file = glob.glob(os.path.join(args.scores_path, '*_rank_001_*.json'))[0] #we know is just 1 file
    print(score_file)
else:
    score_file = args.scores_path

results = calculate_scores(score_file, args.binder_len, args.is_binder_second)

output_path = args.output or os.path.join(os.path.dirname(score_file), 'score_results.json')
with open(output_path, 'w') as outfile:
    json.dump(results, outfile, indent=4)

print(f"Results saved to: {output_path}")
