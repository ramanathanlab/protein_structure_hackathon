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
from typing import List, Tuple, Dict, Any
import functools
from concurrent.futures import ProcessPoolExecutor
import json


def run_tmalign(pdbs: Tuple[str, str], tmalign_path: str) -> Tuple[float, float]:
    pdb1, pdb2 = pdbs
    cmd = f"{tmalign_path} {pdb1} {pdb2} -outfmt 2"
    proc = subprocess.run(cmd.split(), capture_output=True)
    out = str(proc.stdout).split("\\")
    tm_score1 = float(out[17].split("= ")[1].split(" (")[0])
    tm_score2 = float(out[18].split("= ")[1].split(" (")[0])
    # (if normalized by length of Chain_1, Chain_2)
    return tm_score1, tm_score2


def run_multiple_tmalign(kwargs: Dict[str, Any]) -> List[Tuple[str, str, float, float]]:

    # Parse arguments
    pdbs: List[Tuple[str, str]] = kwargs.get("pdbs")
    tmalign_path: str = kwargs.get("tmalign_path")
    verbose: bool = kwargs.get("verbose")

    results = []

    for pdb1, pdb2 in tqdm(pdbs, disable=not verbose):
        cmd = f"{tmalign_path} {pdb1} {pdb2} -outfmt 2"
        proc = subprocess.run(cmd.split(), capture_output=True)
        out = str(proc.stdout).split("\\")
        # (if normalized by length of Chain_1, Chain_2)
        tm_score1 = float(out[17].split("= ")[1].split(" (")[0])
        tm_score2 = float(out[18].split("= ")[1].split(" (")[0])
        results.append((pdb1, pdb2, tm_score1, tm_score2))
    return results


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
        for i, scores in tqdm(
            enumerate(pool.map(fn, pairs, chunksize=chunksize)),
            total=total_num,
        ):
            results.append((pairs[i], scores))

    with open(out_path, "w") as f:
        json.dump(results, f)


def all_v_all_v2(
    all_pdbs: List[str],
    out_path: str,
    tmalign_path: str = "TMalign",
    num_workers: int = 1,
):
    pairs = list(itertools.combinations(all_pdbs, 2))
    chunksize = max(1, len(pairs) // num_workers)
    chunks = [pairs[i * chunksize : (i + 1) * chunksize] for i in range(chunksize)]
    verbose = [True] + [False] * (chunksize - 1)
    kwargs = [
        {"pdbs": chunk, "verbose": v, "tmalign_path": tmalign_path}
        for chunk, v in zip(chunks, verbose)
    ]
    results = []
    with ProcessPoolExecutor(max_workers=num_workers) as pool:
        for result in pool.map(run_multiple_tmalign, kwargs):
            results.extend(result)

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
    # all_pdbs = all_pdbs[:100]  # TODO: testing
    num_workers = 64

    all_v_all_v2(
        all_pdbs, "test.json", num_workers=num_workers, tmalign_path=tmalign_path
    )
