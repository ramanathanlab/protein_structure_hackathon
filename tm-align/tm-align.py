import os
import pickle
import re
import subprocess
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import List, Tuple, Union


"""
This script takes a directory of PDBS and performs pairwise comparisons of TM-scores with them. 

TODO: add all gather for individual tm_score files

Contact Kyle Hippe khippe@anl.gov with questions
"""

PathLike = Union[Path, str]

node_rank = int(os.environ.get("NODE_RANK", 0))  # zero indexed
num_nodes = int(os.environ.get("NRANKS", 1))


def run_tmalign(pdb_pair: Tuple[Path], pattern: re.Pattern) -> float:
    pdb1, pdb2 = pdb_pair
    cmd = f"TMalign {str(pdb1)} {str(pdb2)}"
    res = subprocess.run(cmd.split(), capture_output=True)
    # Score 1 is normalized by pdb 2 score 2 is normalized by pdb 1
    scores = pattern.findall(res.stdout.decode("utf-8"))

    return {
        (str(pdb1), str(pdb2)): float(scores[0]),
        (str(pdb2), str(pdb1)): float(scores[-1]),
    }


def pairwise_processing(process_pdbs: List[Tuple[PathLike]], out_file: Path):
    pattern = re.compile("TM-score= ([+-]?[0-9]*[.]?[0-9]+)")

    pairwise_tmalign = partial(run_tmalign, pattern=pattern)
    scores = []

    print_freq = 1000
    chunk_size = 10000
    end_idx = len(process_pdbs) - (len(process_pdbs) % -chunk_size)
    assert end_idx > len(
        process_pdbs
    ), f"End idx not set properly {end_idx=}, {len(process_pdbs)=}"

    for chunk_start in range(0, end_idx, chunk_size):
        chunk_end = chunk_start + chunk_size
        with ProcessPoolExecutor() as pool:
            for i, res in enumerate(
                pool.map(pairwise_tmalign, process_pdbs[chunk_start:chunk_end])
            ):
                if i % print_freq == 0:
                    print(f"Completed {i+chunk_start} iterations on node {node_rank}")
                scores.append(res)

    print(f"Completed {i+chunk_start} iterations on node {node_rank}")

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

    parser.add_argument(
        "--nogather",
        action="store_false",
        help="If set, will NOT gather all files into a json file",
    )

    args = parser.parse_args()

    pdbs = list(args.input_dir.glob(args.glob_pattern))

    combinations = []
    for i in range(len(pdbs)):
        for j in range(i, len(pdbs)):
            combinations.append((pdbs[i], pdbs[j]))

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
        node_data = combinations
        file_idx = ""

    args.out_dir.mkdir(exist_ok=True, parents=True)
    out_file = args.out_dir / f"tmscore_combinations{file_idx}.pkl"

    pairwise_processing(node_data, out_file)

    if node_rank == 0 and not args.nogather:
        import json
        from collections import defaultdict

        scores = []

        for score_file in args.out_dir.glob("*.pkl"):
            scores.extend(pickle.load(score_file.open("rb")))

        scores_mapping = defaultdict(dict)

        for score_dict in scores:
            for pair, score in score_dict.items():
                pdb1, pdb2 = pair

                pdb1_name = Path(pdb1).stem
                pdb2_name = Path(pdb2).stem

                scores_mapping[pdb1_name][pdb2_name] = score

        gather_file = args.out_dir / f"{args.out_dir.stem}.json"
        json.dump(scores_mapping, gather_file.open("w"))
