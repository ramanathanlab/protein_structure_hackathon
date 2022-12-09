from molr.src.featurizer import MolEFeaturizer
import numpy as np

vanillin = 'COc1cc(C=O)ccc1O'
vanillylamine = 'COc1cc(CN)ccc1O'
vanillyl_alcohol = 'COC1=CC(CO)=CC=C1O'
CoA = 'O=C(NCCS)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3OP(=O)(O)O'
capsaicin = 'CC(C)/C=C/CCCC/C(=N/Cc1ccc(c(c1)OC)O)/O'
capsiate = 'CC(C)C=CCCCCC(=O)OCC1=CC(=C(C=C1)O)OC'

reactions_smiles = [
    [(vanillin,), (vanillylamine,)],
    [(vanillin,), (vanillyl_alcohol,)],
    [(vanillylamine, CoA), (capsaicin,)],
    [(vanillyl_alcohol, CoA), (capsiate,)]
]

# For each reaction, give the index for every molecule in the reaction. 
# TODO: Can we make this faster?
reactions_smiles_idx = []
molecules = []
for reaction in reactions_smiles:
    reagents_idx = []  # reactants and products list
    for reagents in reaction:
        reagent_idx_tuple = [] # reactant tuple (or product tuple)
        for reagent in reagents:
            if reagent not in molecules:             
                molecules.append(reagent)
            reagent_idx_tuple.append(molecules.index(reagent))
        reagents_idx.append(tuple(reagent_idx_tuple))
    reactions_smiles_idx.append(reagents_idx)

model = MolEFeaturizer(path_to_model='../saved/gcn_1024')
embeddings, flags = model.transform(molecules)

# Embedding for each reaction
reaction_embeddings = []
for reaction in reactions_smiles_idx:
    reactants_ids = reaction[0]
    reactants_vector = sum([embeddings[idx] for idx in reactants_ids])

    products_ids = reaction[1]
    transformation_vector = sum([embeddings[idx] for idx in products_ids]) - reactants_vector
    reaction_embeddings.append((reactants_vector, transformation_vector))

# Same substrate different reaction
dist_0 = (np.linalg.norm(reaction_embeddings[0][0] - reaction_embeddings[1][0]), 
          np.linalg.norm(reaction_embeddings[0][1] - reaction_embeddings[1][1]))

# Different substrates same reaction
dist_1 = (np.linalg.norm(reaction_embeddings[2][0] - reaction_embeddings[3][0]), 
          np.linalg.norm(reaction_embeddings[2][1] - reaction_embeddings[3][1]))