"""
This script takes a list of pdbs and compares them to two reference pdbs.

This can be easily adapted to run pairwise tm-scores.

Contact Kyle Hippe khippe@anl.gov with questions
"""
import subprocess
import itertools
import shutil
from pathlib import Path
from tqdm import tqdm
from typing import List, Tuple
import functools
import re
from concurrent.futures import ProcessPoolExecutor
import json
from itertools import chain
from operator import itemgetter


def run_tmalign(pdbs: Tuple[str, str], tmalign_path: str) -> Tuple[float, float]:
    pdb1, pdb2 = pdbs
    cmd = f"{tmalign_path} {pdb1} {pdb2} -outfmt 2"
    proc = subprocess.run(cmd.split(), capture_output=True)
    out = proc.stdout.decode("utf-8")
    print(out)
    # tm_score1, tm_score2
    return float(out[13][1:]), float(out[14][1:])


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


def all_v_all(
    all_pdbs: List[str],
    out_path: str,
    tmalign_path: str = "TMalign",
    num_workers: int = 1,
):
    total_num = len(all_pdbs)
    total_num = total_num * (total_num - 1) // 2

    pairs = list(itertools.combinations(all_pdbs, 2))
    fn = functools.partial(run_tmalign, tmalign_path=tmalign_path)
    results = []
    chunksize = max(1, total_num // num_workers)
    with ProcessPoolExecutor(max_workers=num_workers) as pool:
        for scores in tqdm(
            pool.map(fn, pairs, chunksize=chunksize),
            total=total_num,
        ):
            results.append(scores)

    with open(out_path, "w") as f:
        json.dump(results, f)


if __name__ == "__main__":

    input_pdb_dir = Path(
        "/lus/eagle/projects/CVD-Mol-AI/hippekp/visualization_structures/mdh/mdh_structures_transfer"
    )
    tmalign_path = Path(
        "/lus/eagle/projects/CVD-Mol-AI/braceal/conda/envs/tm-align/bin/TMalign"
    )

    node_local_path = Path("/tmp")
    node_local_pdb_dir = node_local_path / input_pdb_dir.name
    if not node_local_pdb_dir.exists():
        shutil.copytree(input_pdb_dir, node_local_pdb_dir)
    tmalign_path = str(shutil.copy(tmalign_path, node_local_path / tmalign_path.name))

    input_pdb_dir = node_local_pdb_dir
    all_pdbs = list(map(str, input_pdb_dir.glob("*.pdb")))
    all_pdbs = all_pdbs[:100]  # TODO: testing
    num_workers = 64

    all_v_all(all_pdbs, "test.json", num_workers=num_workers, tmalign_path=tmalign_path)
