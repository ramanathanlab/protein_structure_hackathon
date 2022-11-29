import os
import pickle
import re
import subprocess
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from pathlib import Path
from typing import List, Tuple, Union

from tqdm import tqdm

"""
This script takes a list of pdbs and compares them to two reference pdbs.

This can be easily adapted to run pairwise tm-scores.

Contact Kyle Hippe khippe@anl.gov with questions
"""

PathLike = Union[Path, str]

node_rank = int(os.environ.get("NODE_RANK", 0))  # zero indexed
num_nodes = int(os.environ.get("NRANKS", 1))


def run_tmalign(pdb1: Path, pdb2: Path, pattern: re.Pattern) -> float:
    cmd = f"tmalign {str(pdb1)} {str(pdb2)}"
    res = subprocess.run(cmd.split(), capture_output=True)
    # Score 1 is normalized by pdb 2 score 2 is normalized by pdb 1
    scores = pattern.findall(res.stdout.decode("utf-8"))

    return {(pdb1, pdb2): float(scores[0]), (pdb2, pdb1): float(scores[-1])}


def pairwise_processing(process_pdbs: List[Tuple[PathLike]], out_file: Path):
    pattern = re.compile("TM-score= ([+-]?[0-9]*[.]?[0-9]+)")

    pairwise_tmalign = partial(run_tmalign, pattern=pattern)
    futures = []
    with ProcessPoolExecutor() as pool:
        for (pdb1, pdb2) in tqdm(process_pdbs):
            futures.append(pool.submit(pairwise_tmalign, pdb1=pdb1, pdb2=pdb2))

    scores = []
    print_freq = 1000
    for i, fut in enumerate(as_completed(futures)):
        if i % print_freq == 0:
            print(f"Completed {i} iterations on node {node_rank}")
        scores.append(fut.result())

    pickle.dump(scores, out_file.open("wb"))


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--input_dir", type=Path, required=True, help="Folder with PDBs"
    )
    parser.add_argument(
        "-o", "--out_dir", type=Path, required=True, help="Path to save chunks to"
    )

    parser.add_argument(
        "-g",
        "--glob_pattern",
        default="*.pdb",
        help="Glob pattern to find input pdbs in folder",
    )

    args = parser.parse_args()

    pdbs = list(args.input_dir.glob(args.glob_pattern))

    combinations = []
    for i in range(len(pdbs)):
        for j in range(i, len(pdbs)):
            combinations.append((pdbs[i], pdbs[j]))

    # print(len(combinations))
    # exit()
    if num_nodes > 1:
        chunk_size = max(len(combinations) // num_nodes, 1)
        start_idx = node_rank * chunk_size
        end_idx = start_idx + chunk_size
        if node_rank + 1 == num_nodes:
            end_idx = len(combinations)

        print(
            f"Node {node_rank} / {num_nodes} starting at {start_idx}, ending at {end_idx} ({len(combinations)=})"
        )
        node_data = combinations[start_idx:end_idx]
        file_idx = f"_{start_idx}-{end_idx}"
    else:
        print(f"Node {node_rank} processing all {len(combinations)} combinations")
        node_data = combinations[:20000]
        file_idx = ""

    args.out_dir.mkdir(exist_ok=True, parents=True)
    out_file = args.out_dir / f"tmscore_combinations{file_idx}.pkl"

    pairwise_processing(node_data, out_file)
