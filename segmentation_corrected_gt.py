#!/usr/bin/env python3
import argparse
import sys
import gzip
import io
import os
import multiprocessing as mp
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")


def tsv_by_chrom(path):
    """
    Yield (chrom, pos_array, cov_array) from TSV `chr  pos  cov` (can be .gz).
    Lines with cov='NA' -> np.nan.
    """
    openf = gzip.open if path.endswith((".gz", ".bgz")) else open
    with openf(path, "rt", buffering=io.DEFAULT_BUFFER_SIZE * 4) as fh:
        cur_chr, pos_list, cov_list = None, [], []
        for line in fh:
            chrom, p, v = line.rstrip("\n").split("\t")
            if cur_chr is not None and chrom != cur_chr:
                yield cur_chr, np.asarray(pos_list, np.int32), np.asarray(cov_list, np.float32)
                pos_list, cov_list = [], []
            cur_chr = chrom
            pos_list.append(int(p))
            cov_list.append(np.nan if v.upper() == "NA" else float(v))
        if cur_chr is not None:
            yield cur_chr, np.asarray(pos_list, np.int32), np.asarray(cov_list, np.float32)


def process_chrom(args):
    (
        chrom, pos, cov,
        zthr, max_gap,
        geno_mode, ratio_thr,
        flank_size, flank_zmax, min_flank_bins
    ) = args

    # --- robust MAD / Z-score computation ---
    med = np.nanmedian(cov)
    mad = np.nanmedian(np.abs(cov - med)) * 1.4826

    # fallback to (scaled) std if MAD is zero/non-finite
    if not (np.isfinite(mad) and mad > 0):
        mad = np.nanstd(cov) * 1.4826

    if not (np.isfinite(mad) and mad > 0):
        # no variability (all-equal or all-NaN): produce zeros (no signals),
        # but preserve NaNs from original coverage.
        z = np.zeros_like(cov, dtype=np.float32)
        z[np.isnan(cov)] = np.nan
    else:
        z = (cov - med) / mad
        # keep NaNs where coverage is NaN
        z[np.isnan(cov)] = np.nan

    # bedGraph lines (per-base zscore; NA where z is nan)
    bg = [
        f"{chrom}\t{int(p)}\t{int(p)+1}\t{'NA' if np.isnan(z0) else f'{z0:.6f}'}"
        for p, z0 in zip(pos, z)
    ]

    # find signal positions: z not NaN and abs(z) > threshold
    mask = (~np.isnan(z)) & (np.abs(z) > zthr)
    sig_pos = pos[mask]

    # vectorized segmentation by max_gap
    seg_coords = []
    if sig_pos.size:
        if sig_pos.size == 1:
            seg_coords = [(int(sig_pos[0]), int(sig_pos[0]) + 1)]
        else:
            diffs = np.diff(sig_pos)
            # indices where gap > max_gap
            breaks = np.where(diffs > max_gap)[0]
            starts = np.concatenate(([sig_pos[0]], sig_pos[breaks + 1]))
            ends = np.concatenate((sig_pos[breaks], [sig_pos[-1]])) + 1
            seg_coords = [(int(s), int(e)) for s, e in zip(starts, ends)]

    # genotyping per-segment
#    seg = []
#    for s, e in seg_coords:
#        mask_seg = (pos >= s) & (pos < e)
#        seg_z_mean = np.nanmean(z[mask_seg])
#
#        if geno_mode == "coverage-minmax":
#            seg_cov = cov[mask_seg]
#            mn, mx = np.nanmin(seg_cov), np.nanmax(seg_cov)
#            metric = (mn + mx) / 2.0
#            genotype = "0/1" if np.nanmean(seg_cov) > metric else "1/1"
#            #genotype = "0/0" if np.nanmean(seg_cov) > metric else "0/1"

#        else:  # geno_mode == "zscore-ratio"
#            mask_flank = (
#                ((pos >= s - flank_size) & (pos < s)) |
#                ((pos >= e) & (pos < e + flank_size))
#            )
#            flank_z = z[mask_flank]
#            bad_flank = (
#                np.count_nonzero(~np.isnan(flank_z)) < min_flank_bins or
#                (np.nanmean(np.abs(flank_z)) > flank_zmax if np.count_nonzero(~np.isnan(flank_z)) > 0 else True)
#            )
#            if not bad_flank:
#                flank_mean = np.nanmean(np.abs(flank_z))
#                metric = np.abs(seg_z_mean) / flank_mean if flank_mean else np.nan
#                genotype = "1/1" if metric > ratio_thr else "0/1"
#            else:
#                # fallback to coverage minmax if flanks are bad
#                seg_cov = cov[mask_seg]
#                mn, mx = np.nanmin(seg_cov), np.nanmax(seg_cov)
#                metric = (mn + mx) / 2.0
#                genotype = "1/1" if np.nanmean(seg_cov) > metric else "0/1"
#                print(f"WARN\t{chrom}:{s}-{e}\tfallback=minmax", file=sys.stderr)
#
#        seg.append(
#            f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t{genotype}\t{metric:.6f}"
#        )
    # genotyping per-segment: сначала проверяем наличие сигнала в сегменте и валидность фланков
    seg = []
    EPS = 0.0 #1e-12 # псевдоноль
    MIN_SEG_BINS = 1  # минимум валидных точек внутри сегмента (можно поднять при желании)
    NOISE_RATIO = 2.0 # фланки считаем "рваными", если std(|z|) > NOISE_RATIO * mean(|z|)

    def flanks_ok_and_level(s, e):
        """
        Возвращает (ok, reason, flank_mean_cov, flank_mean_absz).
        ok=False и  {"few","noisy","nocov"} — причина отсутствия сигналла.
        """
        mask = (
            ((pos >= s - flank_size) & (pos < s)) |
            ((pos >= e) & (pos < e + flank_size))
        )
        zf = z[mask]
        cf = cov[mask]

        # 1) достаточно валидных бинов?
        n = int(np.count_nonzero(~np.isnan(zf)))
        if n < min_flank_bins:
            return False, "few", np.nan, np.nan

        # 2) нет ли шума
        absz = np.abs(zf)
        mean_absz = float(np.nanmean(absz))
        std_absz  = float(np.nanstd(absz))
        # если mean_absz == 0, std тоже должен быть 0
        if mean_absz == 0.0:
            is_noisy = std_absz > 0.0
        else:
            is_noisy = std_absz > NOISE_RATIO * mean_absz
        if is_noisy:
            return False, "noisy", np.nan, mean_absz

        # 3) во фланках должно быть ненулевое покрытие (не тех.дыра)
        mean_cov = float(np.nanmean(cf))
        if np.isnan(mean_cov) or mean_cov <= EPS:
            return False, "nocov", mean_cov, mean_absz

        return True, "", mean_cov, mean_absz

    for s, e in seg_coords:
        mask_seg = (pos >= s) & (pos < e)
        seg_cov = cov[mask_seg]
        seg_z = z[mask_seg]
        seg_z_mean = float(np.nanmean(seg_z))

        # --- Сегмент: есть ли вообще валидные точки? ---
        n_seg = int(np.count_nonzero(~np.isnan(seg_cov)))
        if n_seg < MIN_SEG_BINS:
            # внутри сегмента нет данных -> ./.
            seg.append(f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t./.\t{np.nan:.6f}")
            continue

        # --- Фланки: единая валидация ---
        ok, reason, flank_mean_cov, flank_mean_absz = flanks_ok_and_level(s, e)
        if not ok:
            print(f"WARN\t{chrom}:{s}-{e}\tno-call:bad-flanks:{reason}", file=sys.stderr)
            seg.append(f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t./.\t{np.nan:.6f}")
            continue

        # --- Генотипирование при валидных фланках ---
        if geno_mode == "coverage-minmax":
            mn = float(np.nanmin(seg_cov))
            mx = float(np.nanmax(seg_cov))
            if not (np.isfinite(mn) and np.isfinite(mx)):
                seg.append(f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t./.\t{np.nan:.6f}")
                continue

            metric   = (mn + mx) / 2.0
            mean_cov = float(np.nanmean(seg_cov))
            genotype = "0/1" if mean_cov > metric else "1/1"
            seg.append(f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t{genotype}\t{metric:.6f}")

        #else:  # geno_mode == "zscore-ratio"
        #    # отношение силы сигнала сегмента к «обычному» |z| во фланках
        #    if not (np.isfinite(flank_mean_absz) and flank_mean_absz > 0):
        #        seg.append(f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t./.\t{np.nan:.6f}")
        #        continue
        #    metric = abs(seg_z_mean) / flank_mean_absz
        #    if not np.isfinite(metric):
        #        seg.append(f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t./.\t{np.nan:.6f}")
        #        continue
        #    genotype = "1/1" if metric > ratio_thr else "0/1"
        #    seg.append(f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t{genotype}\t{metric:.6f}")

    return bg, seg


def main():
    p = argparse.ArgumentParser(
        description="Z segmentation of coverage and genotyping"
    )
    p.add_argument("cov_tsv", help="TSV: chr pos cov  (can be .gz)")
    p.add_argument("--bg-out", default="zscore.bedGraph")
    p.add_argument("--bed-out", default="zscore_segments.bed")
    p.add_argument("-z", "--zThreshold", type=float, default=10)
    p.add_argument("-g", "--maxGap", type=int, default=100)
    p.add_argument("-t", "--threads", type=int, default=6)
    p.add_argument("--geno-mode", choices=["coverage-minmax", "zscore-ratio"], default="coverage-minmax",
                   help="genotyping mode: 'coverage-minmax' or 'zscore-ratio'")
    p.add_argument("--ratio-threshold", type=float, default=2.0)
    p.add_argument("--flank-size", type=int, default=2000)
    p.add_argument("--flank-zmax", type=float, default=1.5)
    p.add_argument("--min-flank-bins", type=int, default=10)
    args = p.parse_args()

    # use threads as requested (no artificial cap)
    pool = mp.Pool(int(args.threads))
    tasks = (
        (
            chrom, pos, cov,
            args.zThreshold, args.maxGap,
            args.geno_mode, args.ratio_threshold,
            args.flank_size, args.flank_zmax, args.min_flank_bins
        )
        for chrom, pos, cov in tsv_by_chrom(args.cov_tsv)
    )

    with open(args.bg_out, "w") as bgf, open(args.bed_out, "w") as bedf:
        for bg_lines, seg_lines in pool.imap_unordered(process_chrom, tasks, chunksize=1):
            bgf.write("\n".join(bg_lines) + "\n")
            if seg_lines:
                bedf.write("\n".join(seg_lines) + "\n")
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()


