import os,sys

from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af import mk_af_model
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Run AFDesign with MPNN sampling')
parser.add_argument('pdb', type=str, help='input pdb file path')
parser.add_argument('output_folder', type=str, help='output folder path')
parser.add_argument('target_chains', type=str, help='target chain IDs (comma-separated)')
parser.add_argument('binder_chains', type=str, help='binder chain IDs (comma-separated)')
parser.add_argument('--num_recycles', type=int, default=12, help='number of repacking cycles (default: 12)')
parser.add_argument('--num_seqs', type=int, default=1000, help='number of sequences to generate (default: 1000)')
parser.add_argument('--num_filt_seq', type=int, default=50, help='number of sequences to keep after filtering (default: 50)')
parser.add_argument('--sampling_temp', type=float, default=0.15, help='sampling temperature (default: 0.15)')
parser.add_argument('--results_dataframe', type=str, help='save results')
parser.add_argument('--break_design', action='store_true', help='stop making new sequences if rmsd < 2 or plddt > 0.85')
args = parser.parse_args()

pdb = args.pdb
output_folder = args.output_folder
target_chains = args.target_chains.split(',')
binder_chains = args.binder_chains.split(',')
num_recycles = args.num_recycles
num_seqs = args.num_seqs
num_filt_seq = args.num_filt_seq
sampling_temp = args.sampling_temp
results_dataframe = args.results_dataframe

### Initialize ###
protocol = "binder"
if protocol == "binder":
    af_terms = ["plddt","i_ptm","i_pae","rmsd"]

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


#Generate ProteinMPNN (100-1000)
print("starting ProteinMPNN")
mpnn_model = mk_mpnn_model()
mpnn_model.get_af_inputs(af_model)
out = mpnn_model.sample(num=num_seqs//10, batch=10, temperature=sampling_temp)

# Filter top X sequences based on score
# Extract the score list from the dictionary and create a sorted index
score_list = out['score']
sorted_index = sorted(range(len(score_list)), key=lambda k: score_list[k], reverse=True) #False --> lower score is first

# Iterate through the dictionary and sort each value (list) using the sorted index
for key, value in out.items():
    sorted_value = [value[i] for i in sorted_index][:num_filt_seq] # slice the sorted list to only keep the first x elements
    out[key] = sorted_value

for k in af_terms: out[k] = [] #?
print(out.keys())

pdb_name = pdb.split('/')[-1].split('.pdb')[0]
loc = f"{output_folder}/af2"
os.system(f"mkdir -p {loc}")

#AF2
print("starting AF2")
#with open(f"{loc}/design.fasta","w") as fasta:
for n in range(num_filt_seq):
  seq = out["seq"][n][-af_model._len:]
  af_model.predict(seq=seq, num_recycles=num_recycles, num_models=1, verbose=False)
  print(af_model.aux["log"])
  for t in af_terms: out[t].append(af_model.aux["log"][t])
  if "i_pae" in out:
    out["i_pae"][-1] = out["i_pae"][-1] * 31
  if "pae" in out:
    out["pae"][-1] = out["pae"][-1] * 31


  current_model=f"{loc}/{pdb_name}_{n}.pdb"
  af_model.save_current_pdb(f"{current_model}")
  af_model._save_results(save_best=True, verbose=False)
  af_model._k += 1
  score_line = [f'mpnn:{out["score"][n]:.3f}']
  for t in af_terms:
    score_line.append(f'{t}:{out[t][n]:.3f}')
  print(n, " ".join(score_line)+" "+seq)
  line = f'>{"|".join(score_line)}\n{seq}'


  if args.break_design:
    if out["plddt"][n] > 0.9 and out["rmsd"][n] < 2:
      print('Breaking out of loop due to --break_design flag')
      num_filt_seq=n #update the number to save later
      break

  
    #fasta.write(line+"\n")

#af_model.save_pdb(f"{loc}/best.pdb")

model_paths = [f"{loc}/{pdb_name}_{n}.pdb" for n in range(num_filt_seq)]

labels = ["score"] + af_terms + ["seq"] #+ ["model_path"]

data = [[out[k][n] for k in labels] + [model_paths[n], pdb] for n in range(num_filt_seq)]

# add additional labels
labels = ["score"] + af_terms + ["seq"] + ["model_path"] + ["input_pdb"]
labels[0] = "mpnn"


import fcntl
import time

df = pd.DataFrame(data, columns=labels)

if results_dataframe:
   output_path = f'{results_dataframe}/af2_results.csv'
else:
   output_path = f'{output_folder}/af2_results.csv'

if os.path.isfile(output_path):
    with open(output_path, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
        time.sleep(np.random.uniform(0, 0.05))
        df.to_csv(f, header=False, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)  # release the lock
else:
    with open(output_path, 'w') as f:
        fcntl.flock(f, fcntl.LOCK_EX)  # lock the file for exclusive access
        time.sleep(np.random.uniform(0, 0.05))
        df.to_csv(f, index=False)
        fcntl.flock(f, fcntl.LOCK_UN)  # release the lock
