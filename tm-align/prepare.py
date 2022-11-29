from pathlib import Path
from itertools import combinations_with_replacement
from tqdm import tqdm
import subprocess
from multiprocessing import Pool
import os

mdh_pdbs_path = Path(
    "/lus/eagle/projects/CVD-Mol-AI/hippekp/visualization_structures/mdh/mdh_structures_transfer"
)
result_path = "chunk_scores"
os.makedirs(result_path, exist_ok=True)

blocksize = 1000

all_pdbs = list(map(str, mdh_pdbs_path.glob("*.pdb")))
all_pdbs.sort()
with open("./listOfFiles.list", "wb") as f:
    lines = [f"{each}\n".encode() for each in all_pdbs]
    f.writelines(lines)


all_pairs = list(combinations_with_replacement(all_pdbs, 2))


def run_against(chunk_file):
    # this runs in a separate thread
    try:
        subprocess.run(["/usr/bin/bash", "./script.sh", chunk_file], check=True)
    except e:
        print(chunk_file, e)


pool = Pool()

pbar = tqdm(total=len(all_pairs), desc="pair progress")

for i in tqdm(range(0, len(all_pairs), blocksize), desc="chunk dispatch"):

    start = i
    end = min(i + blocksize, len(all_pairs))
    chunk_file = f"./{result_path}/allPairs_{start}_{end-1}.list"

    with open(chunk_file, "wb") as f:
        block = all_pairs[start:end]
        lines = [f"{pair[0]} {pair[1]}\n".encode() for pair in block]
        f.writelines(lines)

    pool.apply_async(run_against, args=(chunk_file,), callback=pbar.update(len(block)))

pool.close()
pool.join()

pbar.close()
