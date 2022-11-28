"""
This script takes a list of pdbs and compares them to two reference pdbs.

This can be easily adapted to run pairwise tm-scores.

Contact Kyle Hippe khippe@anl.gov with questions
"""
import subprocess
import itertools
from pathlib import Path
from tqdm import tqdm
import re
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import json
from itertools import chain
from operator import itemgetter


pattern = re.compile("TM-score= ([+-]?[0-9]*[.]?[0-9]+)")


def run_tmalign(pdbs) -> float:
    pdb1, pdb2 = pdbs
    cmd = f"TMalign {str(pdb1)} {str(pdb2)}"
    res = subprocess.run(cmd.split(), capture_output=True)
    # TODO: Parse this with an index lookup instead of regex
    score = [float(each) for each in pattern.findall(res.stdout.decode("utf-8"))]
    return score
    # return 1, 1


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


def all_v_all(all_pdbs, out_path, num_workers: int = 1):
    total_num = len(all_pdbs)
    total_num = total_num * (total_num - 1) // 2

    results = []
    chunksize = max(1, total_num // num_workers)
    with ProcessPoolExecutor(max_workers=num_workers) as pool:
        for scores in tqdm(
            pool.map(
                run_tmalign, itertools.combinations(all_pdbs, 2), chunksize=chunksize
            ),
            total=total_num,
        ):
            results.append(scores)

    with open(out_path, "w") as f:
        json.dump(results, f)


if __name__ == "__main__":

    mdh_pdbs_path = Path(
        "/lus/eagle/projects/CVD-Mol-AI/hippekp/visualization_structures/mdh/mdh_structures_transfer"
    )
    all_pdbs = list(mdh_pdbs_path.glob("*.pdb"))
    all_pdbs = all_pdbs[:10]  # TODO: testing
    num_workers = 64

    all_v_all(all_pdbs, "test.json", num_workers=num_workers)
