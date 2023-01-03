import os
import modelseedpy
from molr.src.featurizer import MolEFeaturizer
import pandas as pd 
import numpy as np
import pickle
import random

embedding_model_name = 'gcn_1024'

modelseed = modelseedpy.biochem.from_local('../ModelSEED/database/ModelSEEDDatabase')

smiles = []
cp_ids = []
for id in modelseed.compounds:
    cp_ids.append(id)
    smiles.append(modelseed.compounds[id]['smiles'])


# Build (and save) an embedding of for the compounds in ModelSEEDDatabase
if os.path.exists(f'embeddings/{embedding_model_name}_compound_embeddings.pk'):
    compounds_embeddings_file = open(f'embeddings/{embedding_model_name}_compound_embeddings.pk', 'rb')
    compound_embeddings = pickle.load(compounds_embeddings_file)
else: 
    model = MolEFeaturizer(path_to_model=f'../MolR/molr/saved/{embedding_model_name}')
    embeddings, flags = model.transform(smiles[:20000])
    compound_embeddings = {cp_ids[i]: (embeddings[i], flags[i]) for i, _ in enumerate(embeddings)}
    del model

    # We run out of memory at around 20000 embeddings
    model = MolEFeaturizer(path_to_model=f'../MolR/molr/saved/{embedding_model_name}')
    embeddings, flags = model.transform(smiles[20000:])  
    compound_embeddings.update({cp_ids[20000+i]: (embeddings[i], flags[i]) for i, _ in enumerate(embeddings)})

    compounds_embeddings_file = open(f'embeddings/{embedding_model_name}_compound_embeddings.pk', 'wb')
    pickle.dump(compound_embeddings, compounds_embeddings_file)


# For the reactions (which are annotated by compound) calculate embeddings.
if os.path.exists(f'embeddings/{embedding_model_name}_rxn_embeddings.pk'):
    reaction_embeddings_file = open(f'embeddings/{embedding_model_name}_rxn_embeddings.pk', 'rb')
    reaction_embeddings = pickle.load(reaction_embeddings_file)
else: 
    reaction_embeddings = dict()
    for reaction in modelseed.reactions: 
        # Check rxn for nans
        rxn_compound_string = modelseed.reactions[reaction]['compound_ids']
        if not isinstance(rxn_compound_string, str):
            continue
        
        # Check if the rxn's compounds were successfully embedded. 
        rxn_compounds = rxn_compound_string.split(';')  
        if not all([compound_embeddings[rc][1] for rc in rxn_compounds]):
            continue
        
        stoichiometry_string = modelseed.reactions[reaction]['stoichiometry']
        if not isinstance(stoichiometry_string, str):
            continue

        stoichiometry = {stoich.split(':')[1]:float(stoich.split(':')[0]) for stoich in stoichiometry_string.split(';')}
        substrate_embedding = -1 * sum([min(0, stoichiometry[rc])*compound_embeddings[rc][0] for rc in rxn_compounds])
        transformation_embedding = sum([stoichiometry[rc]*compound_embeddings[rc][0] for rc in rxn_compounds])

        reaction_embeddings.update({reaction:(substrate_embedding, transformation_embedding)})

    reaction_embeddings_file = open(f'embeddings/{embedding_model_name}_rxn_embeddings.pk', 'wb')
    pickle.dump(reaction_embeddings, reaction_embeddings_file)

if __name__ == "__main__":
    # For pairs of reactions calculate distances
    random_keys = random.sample(reaction_embeddings.keys(), 100)
    random_keys = [rk for rk in random_keys if isinstance(modelseed.reactions[rk]['ec_numbers'], str)][:min(len(random_keys), 50)]

    substrate_embeddings_sample = np.vstack([reaction_embeddings[k][0] for k in random_keys])
    transformation_embeddings_sample = np.vstack([reaction_embeddings[k][1] for k in random_keys])

    w = 0.75
    embedded_reactions = w * substrate_embeddings_sample + (1-w) * transformation_embeddings_sample
    embedded_reactions_df = pd.DataFrame(embedded_reactions, index=[modelseed.reactions[k]['ec_numbers'] for k in random_keys])

    # Do distances match EC numbers?
    from scipy.cluster.hierarchy import dendrogram, linkage
    from matplotlib import pyplot as plt

    Z = linkage(embedded_reactions_df, 'single')

    fig, ax = plt.subplots(1, 1, figsize=(8, 6), dpi=100)
    ax.set_title(f'Hierarchical Clustering Dendrogram \n {w} substrate and {1-w} tranformation weights')
    ax.set_xlabel('')
    ax.set_ylabel('distance (Euclidean)')

    # Make the dendrogram
    dendrogram(Z, labels=embedded_reactions_df.index, leaf_rotation=0, orientation='left')
    plt.show()
