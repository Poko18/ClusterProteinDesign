import os,sys

from colabdesign.mpnn import mk_mpnn_model, clear_mem
from colabdesign.af import mk_af_model
import pandas as pd
import numpy as np
import argparse
import fcntl
import time
import matplotlib.pyplot as plt
#from pyrosetta import *
#from pyrosetta.rosetta import *
#from pyrosetta.rosetta.protocols import *


parser = argparse.ArgumentParser(description='Run AFDesign with MPNN sampling')
parser.add_argument('pdb', type=str, help='input pdb file path')
parser.add_argument('output_folder', type=str, help='output folder path')
parser.add_argument('target_chains', type=str, help='target chain IDs (comma-separated)')
parser.add_argument('binder_chains', type=str, help='binder chain IDs (comma-separated)')
parser.add_argument('iterations', type=int, default=10, help='number of ProteinMPNN/AF2 iterations')
parser.add_argument('--initial_proteinmpnns', type=int, default=100, help='initial proteinmpnns (default: 4)')
parser.add_argument('--designs_per_iteration', type=int, default=10, help='number of designs for further optimization (default: 10)')
parser.add_argument('--proteinmpnn_per_input', type=int, default=4, help='proteinmpnn_per_input pdb in iteration (default: 4)')
parser.add_argument('--sort_by', type=str, default="i_pae", help='score to sort during iteration (default: i_pae)')
parser.add_argument('--num_recycles', type=int, default=12, help='number of af2 cycles (default: 12)')
parser.add_argument('--sampling_temp', type=float, default=0.15, help='sampling temperature (default: 0.15)')
parser.add_argument('--mpnn_batch', type=int, default=5, help='proteinMPNN batch (default: 5)')
args = parser.parse_args()

# Initialize PyRosetta
# init("-beta_nov16 -holes:dalphaball")

pdb = args.pdb
output_folder = args.output_folder
target_chains = args.target_chains.split(',')
binder_chains = args.binder_chains.split(',')
iterations = args.iterations

initial_proteinmpnns=args.initial_proteinmpnns
designs_per_iteration=args.designs_per_iteration
proteinmpnn_per_input=args.proteinmpnn_per_input
sort_by=args.sort_by

num_recycles = args.num_recycles
sampling_temp = args.sampling_temp
mpnn_batch = args.mpnn_batch

t0=time.time()

### Functions ###
def save_to_csv(df,output_path):
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

def create_dataframe(iterations, pdb_output):
    result = []
    labels = ["iterations"]
    iteration_folder = pdb_output
    existing_keys = set()  # Keep track of keys already added to labels

    for iteration in range(1, iterations + 1):
        try:
            data = pd.read_csv(f"{iteration_folder}/it_{iteration}.csv")
            keys = data.head(10).describe().keys().tolist()

            # Add keys to labels if they are not already present
            for key in keys:
                if key not in existing_keys:
                    labels.append(key)
                    existing_keys.add(key)

            row = [iteration]
            for col in keys:
                row.append(data.head(10).describe()[col]["mean"])

            result.append(row)
        except FileNotFoundError:
            continue

    df = pd.DataFrame(result, columns=labels)
    return df

def plot_columns(df, columns, save_path=None):
    # Calculate the number of rows and columns in the grid
    num_plots = len(columns)
    num_cols = math.ceil(math.sqrt(num_plots))
    num_rows = math.ceil(num_plots / num_cols)

    # Create the grid of subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12,6))
    
    # Flatten the axes array for easier indexing
    axes = axes.flatten()

    # Plot each specified column against the "iterations" column
    for i, column in enumerate(columns):
        ax = axes[i]  # Get the current subplot
        ax.plot(df['iterations'], df[column], label=column)

        # Add labels and a legend to each subplot
        ax.set_xlabel('Iterations')
        ax.set_ylabel('Value')
        ax.set_title(column)
        ax.legend()

    # Remove any unused subplots
    if num_plots < len(axes):
        for j in range(num_plots, len(axes)):
            fig.delaxes(axes[j])

    # Adjust the spacing between subplots
    fig.tight_layout()

    # Save the plot if a save path is provided
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

#################

print("running initial proteinMPNN")

### Initialize ###
# output_folder all output
pdb_basename=pdb.split("/")[-1].split(".pdb")[0]
pdb_output=f"{output_folder}/binders/{pdb_basename}"
os.system(f"mkdir -p {output_folder}")
os.system(f"mkdir -p {pdb_output}")

af_terms = ["plddt","i_ptm","i_pae","rmsd"]
other_terms= ["model_path","input_pdb"]
labels = ["score"] + af_terms + other_terms

af_model = mk_af_model(protocol="binder",
                       initial_guess = True,
                       best_metric = "rmsd",
                       use_multimer = True,
                       data_dir="/home/tsatler/projects/AFdesign_playground",
                       model_names = ["model_1_multimer_v3"])
mpnn_model = mk_mpnn_model(weights="soluble")


# Run initial mpnn/af2
initial_csv=f"{pdb_output}/it_0.csv"

if not os.path.exists(initial_csv):
   af_model.prep_inputs(pdb,target_chain=",".join(target_chains),binder_chain=",".join(binder_chains),rm_aa="C")

   mpnn_model.get_af_inputs(af_model)
   mpnn_out = mpnn_model.sample(num=initial_proteinmpnns//mpnn_batch, batch=mpnn_batch, temperature=sampling_temp)
   for k in af_terms: mpnn_out[k] = []
   for term in other_terms: 
      if term not in mpnn_out:
         mpnn_out[term] = []

   # Predict initial sequences
   for n in range(initial_proteinmpnns):
      seq = mpnn_out["seq"][n][-af_model._len:]
      af_model.predict(seq=seq, num_recycles=num_recycles, num_models=1, verbose=False)
      
      for t in af_terms: mpnn_out[t].append(af_model.aux["log"][t])
      if "i_pae" in mpnn_out:
         mpnn_out["i_pae"][-1] = mpnn_out["i_pae"][-1] * 31
      if "pae" in mpnn_out:
         mpnn_out["pae"][-1] = mpnn_out["pae"][-1] * 31
      
      iteration=0
      iteration_folder=f"{pdb_output}/iterations/it_{iteration}"
      os.system(f"mkdir -p {iteration_folder}")
      current_model=f"{iteration_folder}/{pdb_basename}_{iteration}_{n}.pdb"
      mpnn_out["model_path"].append(current_model)
      mpnn_out["input_pdb"].append(pdb)

      af_model.save_current_pdb(f"{current_model}")
      af_model._save_results(save_best=True, verbose=False)
      af_model._k += 1


   iteration_data = [[mpnn_out[k][n] for k in labels] for n in range(initial_proteinmpnns)]
   initial_dataframe=pd.DataFrame(iteration_data, columns=labels)

   if sort_by=="plddt":
      initial_dataframe = initial_dataframe.sort_values("rmsd", ascending=False).reset_index(drop=True)
   else:
      initial_dataframe = initial_dataframe.sort_values(by="rmsd").reset_index(drop=True)

   save_to_csv(initial_dataframe,initial_csv)
else:
    print("Initial dataframe (iteration0) already exists, loading dataframe")
    initial_dataframe = pd.read_csv(initial_csv)



def mpnn_seqs(af_model, pdb, target_chains, binder_chains, num_seqs, rm_aa="C", sampling_temp=sampling_temp, batch=mpnn_batch):
  af_model.prep_inputs(pdb,target_chain=",".join(target_chains),binder_chain=",".join(binder_chains),rm_aa=rm_aa)
  mpnn_model.get_af_inputs(af_model)
  mpnn_out = mpnn_model.sample(num=num_seqs//batch, batch=batch, temperature=sampling_temp)
  for k in af_terms: mpnn_out[k] = []
  for term in other_terms:
      if term not in mpnn_out:
         mpnn_out[term] = []

  return mpnn_out
                    
def af2_predict(af_model, pdb, mpnn_out, output_folder, prefix="", num_models=1):
  
  num_seq=len(mpnn_out["seq"])
  os.makedirs(output_folder, exist_ok=True)

  for n in range(num_seq):
      seq = mpnn_out["seq"][n][-af_model._len:]
      af_model.predict(seq=seq, num_recycles=num_recycles, num_models=num_models, verbose=False)
      
      for t in af_terms: mpnn_out[t].append(af_model.aux["log"][t])
      if "i_pae" in mpnn_out:
         mpnn_out["i_pae"][-1] = mpnn_out["i_pae"][-1] * 31
      if "pae" in mpnn_out:
         mpnn_out["pae"][-1] = mpnn_out["pae"][-1] * 31
      
      current_model=f"{output_folder}/{pdb_basename}_{prefix}_{n}.pdb"
      mpnn_out["model_path"].append(current_model)
      mpnn_out["input_pdb"].append(pdb)

      af_model.save_current_pdb(f"{current_model}")
      af_model._save_results(save_best=True, verbose=False)
      af_model._k += 1
  
  # Prepare data after prediction
  labels = ["score"] + af_terms + other_terms
  iteration_data = [[mpnn_out[k][n] for k in labels] for n in range(num_seq)]

  return iteration_data           

iteration_data=initial_dataframe.head(designs_per_iteration).reset_index(drop=True)

# Start iterations
for i in range(1,iterations+1):
   print(f"starting iteration {i}")
   t_it=time.time()
   
   # Setup iteration folder
   iteration_folder=f"{pdb_output}/iterations/it_{i}"
   os.makedirs(iteration_folder, exist_ok=True)

   # Check if the output CSV file for the current iteration exists
   output_csv_path = f"{pdb_output}/it_{i}.csv"

   if os.path.isfile(output_csv_path):
      iteration_data=pd.read_csv(output_csv_path)
      print(f"Skipping iteration {i} as output file already exists.")
      continue

   # Run MPNNs/AF2s on X number of input designs
   iteration_data = iteration_data.head(designs_per_iteration).reset_index(drop=True)
   best_10=iteration_data
   for z,row in best_10.iterrows():

      mpnn_out=mpnn_seqs(af_model, row["model_path"],"A","B",proteinmpnn_per_input)
      row_data = af2_predict(af_model, row["model_path"], mpnn_out, iteration_folder, prefix=f"{i}_{z}")
      row_dataframe = pd.DataFrame(row_data, columns=labels)

      # Add row data to iteration data
      iteration_data = pd.concat([iteration_data, row_dataframe])
   
   if sort_by=="plddt":
      iteration_data = iteration_data.sort_values(sort_by, ascending=False).reset_index(drop=True)
   else:
      iteration_data = iteration_data.sort_values(sort_by).reset_index(drop=True)
   # Save data
   save_to_csv(iteration_data,output_csv_path) # to iteration file
   save_to_csv(iteration_data,f"{pdb_output}/all.csv") # to combined file
   save_to_csv(iteration_data,f"{output_folder}/all.csv") # to all binders

   # Plot progress
   columns = ['plddt','i_pae']
   try:
      #iterations, iteration_folder, columns=None, save_path=None
      df_plot=create_dataframe(iterations,pdb_output)
      plot_columns(df_plot, columns, save_path=f"{pdb_output}/iteration_plot.png")
   except:
      print("could not plot")
   print(f"finished iteration in {time.time()-t_it} seconds")

print(f"done iterating in {time.time()-t0} seconds")