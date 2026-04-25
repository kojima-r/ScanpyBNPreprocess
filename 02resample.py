"""
Bootstrap-style resampling of per-tissue expression matrices.

For each input file, draw N resamples of size m (either the original
number of cells, "same", or its square root, "root"), average the
sampled rows, and write one row per resample to the output.

When --batched is set, the resampling is stratified by the
"age|batch" key encoded in the cell identifier — one averaged row is
emitted per (resample index, age|batch) combination.
"""

import argparse
import glob
import os
from functools import partial
from multiprocessing import Pool

import numpy as np


def read_matrix(filename):
    """Return (header_line, list_of_row_ids, ndarray)."""
    ids = []
    rows = []
    with open(filename) as fp:
        header = next(fp)
        for line in fp:
            arr = line.split("\t")
            ids.append(arr[0])
            rows.append(np.array(list(map(float, arr[1:]))))
    return header, ids, np.array(rows)


def resample_size(n, mode):
    return int(np.sqrt(n)) if mode == "root" else n


def resample_plain(in_path, out_path, n_resamples, size_mode):
    header, ids, data = read_matrix(in_path)
    m = resample_size(len(ids), size_mode)
    with open(out_path, "w") as ofp:
        ofp.write(header)
        for i in range(n_resamples):
            idx = np.random.randint(0, len(ids), size=m)
            v = np.mean(data[idx, :], axis=0)
            ofp.write(f"s{i}\t" + "\t".join(map(str, v)) + "\n")


def resample_batched(in_path, out_path, n_resamples, size_mode):
    """Stratify by 'age|batch' (positions 1 and 2 of the '|'-split id)."""
    ids = []
    rows = []
    with open(in_path) as fp:
        header = next(fp)
        for line in fp:
            arr = line.split("\t")
            parts = arr[0].split("|")
            ids.append(parts[1] + "|" + parts[2])
            rows.append(np.array(list(map(float, arr[1:]))))
    ids = np.array(ids)
    data = np.array(rows)
    unique_ids = np.unique(ids)
    m = resample_size(len(ids), size_mode)
    with open(out_path, "w") as ofp:
        ofp.write(header)
        for i in range(n_resamples):
            for target_id in unique_ids:
                pool = np.where(ids == target_id)[0]
                idx = np.random.choice(pool, size=m, replace=True)
                v = np.mean(data[idx, :], axis=0)
                ofp.write(f"s{i}|{target_id}\t" + "\t".join(map(str, v)) + "\n")


def _worker(in_path, out_dir, n_resamples, size_mode, batched):
    out_path = os.path.join(out_dir, os.path.basename(in_path))
    print(f"{in_path} >> {out_path}")
    fn = resample_batched if batched else resample_plain
    fn(in_path, out_path, n_resamples, size_mode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", default="01data_bbknn/*.txt",
                        help="Glob pattern for input per-tissue files")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory (default: 02data_bbknn[_b]/)")
    parser.add_argument("-n", "--n-resamples", type=int, default=10,
                        help="Number of resamples per input file")
    parser.add_argument("--size", choices=["same", "root"], default="same",
                        help="Resample size: same as input or its square root")
    parser.add_argument("--batched", action="store_true",
                        help="Stratify resampling by age|batch in the cell id")
    parser.add_argument("--workers", type=int, default=16,
                        help="Process-pool size")
    args = parser.parse_args()

    out_dir = args.out_dir or ("02data_bbknn_b/" if args.batched else "02data_bbknn/")
    os.makedirs(out_dir, exist_ok=True)

    files = sorted(glob.glob(args.input_glob))
    job = partial(_worker,
                  out_dir=out_dir,
                  n_resamples=args.n_resamples,
                  size_mode=args.size,
                  batched=args.batched)
    with Pool(args.workers) as pool:
        pool.map(job, files)


if __name__ == "__main__":
    main()
