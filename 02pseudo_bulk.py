"""
Collapse each per-tissue expression matrix into pseudo-bulk rows by
averaging cells that share a prefix of the "@name" identifier.

The cell identifier in the input is "tissue|age|batch|cell_id" (or
"tissue|age|cell_id" for FACS / droplet). --level controls which
prefix is used as the group key:

  --level tissue  -> "<tissue>"            (1 row per file)
  --level age     -> "<tissue>|<age>"      (1 row per (tissue, age))
  --level batch   -> "<tissue>|<age>|<batch>"  (1 row per (tissue, age, batch))

The first column of the output is the group key itself. The default
output directory is data02_bbknn_p_<level>/.
"""

import argparse
import glob
import os
from functools import partial
from multiprocessing import Pool

import numpy as np


LEVEL_DEPTH = {"tissue": 1, "age": 2, "batch": 3}


def pseudo_bulk(in_path, out_path, level):
    depth = LEVEL_DEPTH[level]

    counts: dict[str, int] = {}
    sums: dict[str, np.ndarray] = {}

    with open(in_path) as fp:
        header = next(fp)
        for line in fp:
            arr = line.rstrip("\n").split("\t")
            key = "|".join(arr[0].split("|")[:depth])
            v = np.fromiter((float(x) for x in arr[1:]), dtype=np.float64)
            if key in sums:
                sums[key] += v
                counts[key] += 1
            else:
                sums[key] = v.copy()
                counts[key] = 1

    with open(out_path, "w") as ofp:
        ofp.write(header)
        for key in sorted(sums.keys()):
            mean = sums[key] / counts[key]
            ofp.write(f"{key}\t" + "\t".join(map(str, mean)) + "\n")


def _worker(in_path, out_dir, level):
    out_path = os.path.join(out_dir, os.path.basename(in_path))
    print(f"{in_path} >> {out_path}")
    pseudo_bulk(in_path, out_path, level)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", default="data01_bbknn/*.txt",
                        help="Glob pattern for input per-tissue files")
    parser.add_argument("--level", choices=list(LEVEL_DEPTH), default="tissue",
                        help="Aggregation level: tissue / age / batch (default: tissue)")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory (default: data02_bbknn_p_<level>/)")
    parser.add_argument("--workers", type=int, default=16,
                        help="Process-pool size")
    args = parser.parse_args()

    out_dir = args.out_dir or f"data02_bbknn_p_{args.level}/"
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(glob.glob(args.input_glob))
    with Pool(args.workers) as pool:
        pool.map(partial(_worker, out_dir=out_dir, level=args.level), files)


if __name__ == "__main__":
    main()
