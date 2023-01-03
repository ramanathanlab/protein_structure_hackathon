import os
import pickle as pkl
import random
from itertools import combinations
import sys
sys.path.append(os.path.dirname(__file__))
from distance_function import rxn_distance
import pandas as pd
from multiprocessing import Pool
from mpi4py import MPI

MPI_RANK = MPI.COMM_WORLD.Get_rank()
MPI_SIZE = MPI.COMM_WORLD.Get_size()

jobs_per_rank = 16  # how many serial jobs to run per rank. 

# run_num = int(sys.argv[1])
model_name = "gcn_1024"
substrate_weighting = 0.25
chunk_size = 400000
output_file = '/vol/structure/ModelSEEDDatabase/reaction_distances/reaction_distances.csv'
output_dir = os.path.dirname(output_file)

# Load data
current_dir = os.path.dirname(__file__)
# TODO: if file doesn't exist..
rxn_embeddings_file = os.path.join(current_dir, f"reaction_embeddings/embeddings/{model_name}_rxn_embeddings.pk")  
with open(rxn_embeddings_file, 'rb') as f:
    rxn_embeddings = pkl.load(f)

# Construct a table of distances for the range of combinations that are assigned to this rank. 
combos_ = list(combinations(rxn_embeddings, 2))
len_combos = len(combos_)
range_combos = combos_[min(len_combos, MPI_RANK*chunk_size):min(len_combos, (MPI_RANK+1)*chunk_size)]

# 
chunk_size = max(len(range_combos) // jobs_per_rank, 1)
chunk_ids = [( i * chunk_size, (i + 1) * chunk_size) for i in range(jobs_per_rank)[:-1]] + [((jobs_per_rank-1) * chunk_size, len(range_combos))]

def f(ids):
    id0, id1 = ids
    distance_data = pd.DataFrame({'rxn1':[], 'rxn2':[], 'distance':[]})  # ID-columns : rxn1 rxn2, Distance-column
    # Iterate over the pairwise combinations of rxns
    for i, (rxn_id1, rxn_id2) in enumerate(range_combos[id0: id1]):
        r1 = rxn_embeddings[rxn_id1]
        r2 = rxn_embeddings[rxn_id2]
        r_dist = rxn_distance(r1, r2, substrate_weighting)
        distance_data.loc[id0+i] = [rxn_id1, rxn_id2, r_dist]
        output_file_wo_ext = os.path.splitext(output_file)[0]
        distance_data.to_csv(f'{output_file_wo_ext}_{MPI_RANK * chunk_size + id0}_{MPI_RANK * chunk_size + id1}.csv', header=True, index=False)

output_dir = os.path.dirname(output_file)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for chunk_id in chunk_ids: 
    f(chunk_id)

