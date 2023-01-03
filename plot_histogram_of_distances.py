import os
import pickle as pkl
from itertools import combinations
import random
from distance_function import rxn_distance

substrate_weighting = 0.25
model_name = "gcn_1024"

# Load data
current_dir = os.path.dirname(__file__)
# TODO: if file doesn't exist..
rxn_embeddings_file = os.path.join(current_dir, f"reaction_embeddings/embeddings/{model_name}_rxn_embeddings.pk")  
with open(rxn_embeddings_file, 'rb') as f:
    rxn_embeddings = pkl.load(f)

# Get the distribution of a sample of the reactions
distances = []
random_keys = random.sample(rxn_embeddings.keys(), 50)
for rxn_id1, rxn_id2 in combinations(random_keys, 2):
    r1 = rxn_embeddings[rxn_id1]
    r2 = rxn_embeddings[rxn_id2]
    r_dist = rxn_distance(r1, r2, substrate_weighting)
    distances.append(r_dist)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from sklearn.neighbors import KernelDensity

n_bins = 20
fig, axs = plt.subplots(1, 2, tight_layout=True)
N, bins, patches = axs.hist(distances, bins=n_bins, density=True)

# Now we format the y-axis to display percentage
axs.yaxis.set_major_formatter(PercentFormatter(xmax=1))
plt.show()

# Gaussian KDE
X = np.reshape(np.array(distances), (len(distances), 1))
X_plot = np.linspace(0, 10000, 1000).reshape(1000, 1)
kde = KernelDensity(kernel="gaussian", bandwidth=250).fit(X)
log_dens = kde.score_samples(X_plot)

fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="#AAAAFF")
ax.set_title("Gaussian Kernel Density")

ax.plot(X[:], np.full(X.shape[0], -0.00002), "+k", alpha=0.1)
plt.show()