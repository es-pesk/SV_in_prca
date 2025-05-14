#!/usr/bin/env python3
import argparse
import sys
import gzip
import io
import multiprocessing as mp

import numpy as np

def parse_depth_file(depth_path):
    """
    Yield (chrom, pos_array, depth_array) for each chromosome block
    in a 3-col bedGraph-like file: chr<tab>pos<tab>depth.
    """
    openf = gzip.open if depth_path.endswith((".gz", ".bgz")) else open
    with openf(depth_path, "rt", buffering=io.DEFAULT_BUFFER_SIZE*4) as fh:
        cur_chr, pos_list, dep_list = None, [], []
        for line in fh:
            chrom, pos, dep = line.rstrip("\n").split("\t")
            if chrom != cur_chr and cur_chr is not None:
                yield cur_chr, np.asarray(pos_list, dtype=np.int32), np.asarray(dep_list, dtype=np.float64)
                pos_list, dep_list = [], []
            cur_chr = chrom
            pos_list.append(int(pos))
            dep_list.append(float(dep))
        if cur_chr is not None:
            yield cur_chr, np.asarray(pos_list, dtype=np.int32), np.asarray(dep_list, dtype=np.float64)

def sliding_cov(depths, win, step):
    """
    Compute sliding CoV = sd/mean over windows of size `win` with shift `step`.
    Returns (cov_vals, center_offsets).
    """
    n = len(depths)
    if n < win:
        return np.array([], dtype=float), np.array([], dtype=int)

    # cumulative sums and sums of squares
    c  = np.concatenate(([0.], np.cumsum(depths,      dtype=np.float64)))
    c2 = np.concatenate(([0.], np.cumsum(depths**2,   dtype=np.float64)))

    # window start indices (0-based)
    idx = np.arange(0, n - win + 1, step, dtype=int)

    # sums and sums of squares over each window
    sum_  = c [idx + win] - c [idx]
    sum2  = c2[idx + win] - c2[idx]

    mean = sum_ / win
    var  = (sum2 / win) - mean**2

    # CoV = sd/mean; R returns NA if mean==0 or all NA
    cov = np.sqrt(np.maximum(var, 0.0)) / mean
    cov[mean == 0] = np.nan

    # compute centers exactly like R:
    # centers = starts + floor((win-1)/2)
    shift = (win - 1) // 2
    centers = idx + shift

    return cov, centers

def process_chrom(args):
    chrom, pos, dep, win, step = args
    cov, centers_idx = sliding_cov(dep, win, step)
    if cov.size == 0:
        return ""
    centers_pos = pos[centers_idx]
    lines = [
        f"{chrom}\t{p}\t{c:.6f}"
        for p, c in zip(centers_pos, cov)
    ]
    return "\n".join(lines)

def main():
    p = argparse.ArgumentParser(
        description="Sliding CoV (sd/mean) over depth file"
    )
    p.add_argument("depth",    help="depth file: chr\\tpos\\tdepth  (can be .gz)")
    p.add_argument("-w", "--win",   type=int, default=100, help="window size W")
    p.add_argument("-s", "--step",  type=int, default=50,  help="shift S")
    p.add_argument("-o", "--out",   default="-", help="output TSV (default stdout)")
    p.add_argument("-t", "--threads", type=int, default=6, help="max CPU cores (≤6)")
    args = p.parse_args()

    pool   = mp.Pool(min(args.threads, 6))
    tasks  = (
        (chrom, pos, dep, args.win, args.step)
        for chrom, pos, dep in parse_depth_file(args.depth)
    )

    writer = sys.stdout if args.out == "-" else open(args.out, "w")
    for chunk in pool.imap_unordered(process_chrom, tasks, chunksize=1):
        if chunk:
            writer.write(chunk + "\n")
    pool.close()
    pool.join()
    if writer is not sys.stdout:
        writer.close()

if __name__ == "__main__":
    main()
