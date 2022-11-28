import subprocess
from pathlib import Path
from tqdm import tqdm
import re
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import json

"""
This script takes a list of pdbs and compares them to two reference pdbs.

This can be easily adapted to run pairwise tm-scores.

Contact Kyle Hippe khippe@anl.gov with questions
"""


def run_tmalign(reference, pdb_1, pattern) -> float:
    cmd = f"tmalign {str(pdb_1)} {str(reference)}"
    res = subprocess.run(cmd.split(), capture_output=True)
    score = float(pattern.findall(res.stdout.decode("utf-8"))[-1])

    return pdb_1, score


def process_against_reference(isoform_path, all_pdbs, out_path):
    pattern = re.compile("TM-score= ([+-]?[0-9]*[.]?[0-9]+)")

    isoform_tmalign = partial(run_tmalign, isoform_path, pattern=pattern)
    isoform_scores = {}
    with ProcessPoolExecutor() as pool:
        for file, score in tqdm(
            pool.map(isoform_tmalign, all_pdbs), total=len(all_pdbs)
        ):
            isoform_scores[str(file)] = score

    with open(out_path, "w") as f:
        json.dump(isoform_scores, f)


if __name__ == "__main__":

    tmscore_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/tmscore/TMscore")
    mdh_pdbs_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/mdh_pdbs")
    isoform_1_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/2pwz.pdb")
    isoform_2_path = Path("/Users/kyle/Desktop/temp/mdh_visualization/1nxu.pdb")

    all_pdbs = list(mdh_pdbs_path.glob("*.pdb"))

    process_against_reference(isoform_1_path, all_pdbs, "2wpz_tmalign_scores.json")
    process_against_reference(isoform_2_path, all_pdbs, "1nxu_tmalign_scores.json")
