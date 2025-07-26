#!/usr/bin/env python3
import argparse, sys, gzip, io, os, multiprocessing as mp
import numpy as np


os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

def tsv_by_chrom(path):
    """
    Yield (chrom, pos_array, cov_array) из TSV `chr  pos  cov` (может быть .gz).
    Строки с cov="NA" → np.nan.
    """
    openf = gzip.open if path.endswith((".gz", ".bgz")) else open
    with openf(path, "rt", buffering=io.DEFAULT_BUFFER_SIZE*4) as fh:
        cur_chr, pos_list, cov_list = None, [], []
        for line in fh:
            chrom, p, v = line.rstrip("\n").split("\t")
            if cur_chr is not None and chrom != cur_chr:
                yield cur_chr, np.asarray(pos_list, np.int32), np.asarray(cov_list, np.float32)
                pos_list, cov_list = [], []
            cur_chr = chrom
            pos_list.append(int(p))
            # читаем “NA” как nan
            cov_list.append(np.nan if v.upper()=="NA" else float(v))
        if cur_chr is not None:
            yield cur_chr, np.asarray(pos_list, np.int32), np.asarray(cov_list, np.float32)


def process_chrom(args):
    chrom, pos, cov, zthr, max_gap, geno_mode, ratio_thr,
    flank_size, flank_zmax, min_flank_bins = args

    
    # 1)Z-score: медиана/MAD×1.4826 или std×1.4826
    med = np.nanmedian(cov)
    mad = np.nanmedian(np.abs(cov - med)) * 1.4826
    if mad == 0:
        mad = np.nanstd(cov) * 1.4826
    # z получим nan там, где cov был nan
    z = (cov - med) / mad

    # 2)bedGraph Z (nan → "NA")
    bg = []
    for p0, z0 in zip(pos, z):
        if np.isnan(z0):
            bg.append(f"{chrom}\t{p0}\t{p0+1}\tNA")
        else:
            bg.append(f"{chrom}\t{p0}\t{p0+1}\t{z0:.6f}")

    # 3)Сегменты по Z: Только |z|>zthr и z≠nan
    mask    = (~np.isnan(z)) & (np.abs(z) > zthr)
    sig_pos = pos[mask]

    seg = []
    if sig_pos.size:
        start = prev = sig_pos[0]
        for p1 in sig_pos[1:]:
            if p1 - prev > max_gap:
                seg_z = np.nanmean(z[(pos >= start) & (pos <= prev)])  # средний Z по сегменту
                seg.append(f\"{chrom}\\t{start}\\t{prev+1}\\t{seg_z:.6f}\") # средний Z по сегменту 4м столбцом в bed
                start = p1
            prev = p1
                
        # финальный сегмент
        seg_z = np.nanmean(z[(pos >= start) & (pos <= prev)])
        seg.append(f\"{chrom}\\t{start}\\t{prev+1}\\t{seg_z:.6f}\")       
        #seg.append(f"{chrom}\t{start}\t{prev+1}")
    ## блок генотипирования
    seg_lines = []
    for s, e in seg_coords:
        mask_seg = (pos >= s) & (pos < e)
        seg_z_mean = np.nanmean(z[mask_seg])
        
        # режим minmax
        if geno_mode == "minmax":                       
            seg_cov = cov[mask_seg]                     
            mn, mx = np.nanmin(seg_cov), np.nanmax(seg_cov)  
            metric = (mn + mx) / 2                   
            genotype = "1/1" if np.nanmean(seg_cov) > metric else "0/1"
            
        # режим Ratio сегмент/фланк (больше 1кб) + fallback 
        else:
            mask_flank = (
                ((pos >= s - flank_size) & (pos < s)) |
                ((pos >= e) & (pos < e + flank_size))
            )
            flank_z = z[mask_flank]
            
            bad_flank = (
                np.count_nonzero(~np.isnan(flank_z)) < min_flank_bins or
                np.nanmean(np.abs(flank_z)) > flank_zmax
            )

            if not bad_flank:
                flank_mean = np.nanmean(np.abs(flank_z))
                metric = np.abs(seg_z_mean) / flank_mean if flank_mean else np.nan
                genotype = "1/1" if metric > ratio_thr else "0/1"
            else:
                seg_cov = cov[mask_seg]
                mn, mx = np.nanmin(seg_cov), np.nanmax(seg_cov)
                metric = (mn + mx) / 2.0
                genotype = "1/1" if np.nanmean(seg_cov) > metric else "0/1"
                print(f"WARN\t{chrom}:{s}-{e}\tfallback=minmax", file=sys.stdout)

        seg.append(
            f"{chrom}\t{s}\t{e}\t{seg_z_mean:.6f}\t{genotype}\t{metric:.6f}"

    
    return bg, seg


def main():
    p = argparse.ArgumentParser(
        description="Сегменты геномного трека CoV - только Z-сегментация "
    )
    p.add_argument("cov_tsv", help="TSV: chr  pos  cov")
    p.add_argument("--bg-out", default="zscore.bedGraph",
                   help="output bedgraph")
    p.add_argument("--bed-out", default="zscore_segments.bed",
                   help="выходной BED с Z-сегментами [default: zscore_segments.bed]")
    p.add_argument("-z", "--zThreshold", type=float, default=10,
                   help="порог abs(Z) [default: 10]")
    p.add_argument("-g", "--maxGap", type=int, default=100,
                   help="макс. gap для объединения [default: 100]")
    p.add_argument("-t", "--threads", type=int, default=6,
                   help="число потоков (≤6) [default: 6]")

    # выбор метода генотипирования
    p.add_argument("--geno-mode", choices=["minmax", "ratio"], default="minmax",
                   help="Как определять 0/1 vs 1/1: minmax | ratio (с fallback)")
    p.add_argument("--ratio-threshold", type=float, default=2.0,
                   help="Порог для ratio‑метода [2.0]")

    # параметры фланков (для ratio)
    p.add_argument("--flank-size", type=int, default=2000, help="Длина фланков для ratio, bp [2000]")
    p.add_argument("--flank-zmax", type=float, default=1.5, help="Макс. |Z̄| во фланках, иначе fallback [1.5]")
    p.add_argument("--min-flank-bins", type=int, default=10, help="Мин. CoV‑окон во фланке [10]")


    
    args = p.parse_args()

    pool = mp.Pool(min(args.threads, 6))
    tasks = (
        (chrom, pos, cov,
            args.zThreshold, args.maxGap,
            args.geno_mode, args.ratio_threshold,
            args.flank_size, args.flank_zmax, args.min_flank_bins
        )
        for chrom, pos, cov in tsv_by_chrom(args.cov_tsv))

    with open(args.bg_out,  "w") as bgf, \
         open(args.bed_out, "w") as bedf:
        for bg_lines, seg_lines in pool.imap_unordered(process_chrom, tasks, chunksize=1):
            bgf.write("\n".join(bg_lines) + "\n")
            if seg_lines:
                bedf.write("\n".join(seg_lines) + "\n")

    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
