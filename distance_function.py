from numpy.linalg import norm


def combine_substrate_and_transformation(substrate_emb1, transformation_emb2, substrate_weight):
    return substrate_weight * substrate_emb1 + (1-substrate_weight) * transformation_emb2


def rxn_distance(rxn_embedding1, rxn_embedding2, substrate_weight, **kwargs):
    rxn1 = combine_substrate_and_transformation(*rxn_embedding1, substrate_weight)
    rxn2 = combine_substrate_and_transformation(*rxn_embedding2, substrate_weight)

    return norm(rxn1-rxn2, **kwargs)
