import subprocess
from pathlib import Path
from tqdm import tqdm
import re
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import json
from itertools import chain
from operator import itemgetter
from multiprocessing import Array

"""
This script takes a list of pdbs and compares them to two reference pdbs.

This can be easily adapted to run pairwise tm-scores.

Contact Kyle Hippe khippe@anl.gov with questions
"""

pattern = re.compile("TM-score= ([+-]?[0-9]*[.]?[0-9]+)")

<<<<<<< HEAD

def run_tmalign(pdb1, pdb2) -> float:
    cmd = f"TMalign {str(pdb2)} {str(pdb1)}"
=======
def run_tmalign(reference, pdb_1, pattern) -> float:
    cmd = f"TMalign {str(pdb_1)} {str(reference)}"
>>>>>>> f8e2d9b6fa9a6941204cc509892cf783374e284f
    res = subprocess.run(cmd.split(), capture_output=True)
    score = [float(each) for each in pattern.findall(res.stdout.decode("utf-8"))]
    return score
    # return 1, 1


# def process_against_reference(isoform_path, all_pdbs, out_path):
#     pattern = re.compile("TM-score= ([+-]?[0-9]*[.]?[0-9]+)")

#     isoform_tmalign = partial(run_tmalign, isoform_path, pattern=pattern)
#     isoform_scores = {}
#     with ProcessPoolExecutor() as pool:
#         for file, score in tqdm(
#             pool.map(isoform_tmalign, all_pdbs), total=len(all_pdbs)
#         ):
#             isoform_scores[str(file)] = score

#     with open(out_path, "w") as f:
#         json.dump(isoform_scores, f)


def one_v_all(all_pdbs, j):
    pdb_j = all_pdbs[j]
    total_num = len(all_pdbs)
    ji_sim = [None] * (j + 1)
    ij_sim = [None] * (j + 1)
    for i in range(0, j + 1):
        if i == j:
            ij_sim[j] = (j, j, 1.0)
            ji_sim[j] = (j, j, 1.0)
        score_ji, score_ij = run_tmalign(pdb_j, all_pdbs[i])
        ji_sim[i] = (j, i, score_ji)
        ij_sim[i] = (i, j, score_ij)
    return j, ij_sim + ji_sim[:j] + ji_sim[j + 1 :]


def all_v_all(all_pdbs, out_path):
    total_num = len(all_pdbs)
    helper_fn = partial(one_v_all, all_pdbs)
    sim_batched = [None] * total_num
    with ProcessPoolExecutor() as pool:
        for idx, sim_list in tqdm(
            pool.map(helper_fn, range(total_num)), total=total_num
        ):
            sim_batched[idx] = sim_list

    sim_all = list(chain(*sim_batched))
    sim_all.sort(key=itemgetter(0, 1))

    with open(out_path, "w") as f:
        json.dump(sim_all, f)


if __name__ == "__main__":
<<<<<<< HEAD

    # tmscore_path = Path("TMscore")
    mdh_pdbs_path = Path(
        "/lus/eagle/projects/CVD-Mol-AI/hippekp/visualization_structures/mdh/mdh_structures_transfer"
    )
    # isoform_1_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/2pwz.pdb")
    # isoform_2_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/1nxu.pdb")
=======
    mdh_pdbs_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/mdh_pdbs")
    isoform_1_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/2pwz.pdb")
    isoform_2_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/1nxu.pdb")
>>>>>>> f8e2d9b6fa9a6941204cc509892cf783374e284f

    all_pdbs = list(mdh_pdbs_path.glob("*.pdb"))
    all_pdbs.sort()

    # all_pdbs = list(range(4))

    all_v_all(all_pdbs, "test.json")

    # process_against_reference(isoform_1_path, all_pdbs, "2wpz_tmalign_scores.json")
    # process_against_reference(isoform_2_path, all_pdbs, "1nxu_tmalign_scores.json")
