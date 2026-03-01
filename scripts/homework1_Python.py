#Write a Python script that runs in the bch709_vibe_coding conda environment.  
#Installed packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm, data.table, ggplot2, pheatmap, viridisLite, scales.

#Input:
#- ~/bch709_vibe_coding/data/mrna.fa.gz (yeast mRNA)

#Task:
#- extract sequence information from the UCSC yeast (Saccharomyces cerevisiae) mRNA FASTA file (mrna.fa.gz)
#- analyze GC content distribution, and produce a summary table and distribution graph
#- look for distinct GC-content subpopulations  

#Output:
#- create a summary table: accession, length, gc_content (4 decimals, sorted by gc_content desc),   GC-content subpopulations.

#- save the summary table to :  ~/bch709_vibe_coding/results/mrna_metrics.tsv
#- create a distribution graph: 	Histogram (1600×900 px, dpi 200)
#- save the distribution graph to :  ~/bch709_vibe_coding/results/gc_content_distribution.png
#- Plot specifications:

#X-axis: GC content (0–1)
#Histogram with density curve overlay
#Mean and median as vertical dashed lines
#Caption showing n, mean, median, sd


#Don't run the script. But, provide the full script, so I can run it in terminal and generate my results


# Interpretation: 
# The GC content of yeast mRNA is lower than eukaryotic mRNA generally. The potential reasons are as follows:
# - 1. High GC content can elevate mutation rates and recombination rates in yeast. This will cause the genome to be more unstable during evolution. Therefore, yeast prefers to have lower GC content in the mRNA.
# - 2. Lower GC content can reduce secondary structures in mRNA, facilitating more efficient tanslation because yeast has a fast growth rate and needs efficient translation.
# - reference: https://doi.org/10.1073/pnas.1807334115

#!/usr/bin/env python3

import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO


INPUT_FASTA = Path.home() / "bch709_vibe_coding" / "data" / "mrna.fa.gz"
OUTPUT_TABLE = Path.home() / "bch709_vibe_coding" / "results" / "mrna_metrics.tsv"
OUTPUT_PLOT = Path.home() / "bch709_vibe_coding" / "results" / "gc_content_distribution.png"


def gc_fraction(seq: str) -> float:
    s = seq.upper()
    valid = sum(base in {"A", "C", "G", "T"} for base in s)
    if valid == 0:
        return np.nan
    gc = s.count("G") + s.count("C")
    return gc / valid


def detect_modes(gc_values: np.ndarray, bins: int = 80) -> np.ndarray:
    counts, edges = np.histogram(gc_values, bins=bins, range=(0, 1))
    centers = (edges[:-1] + edges[1:]) / 2

    # Smooth histogram and find local maxima as candidate subpopulation modes
    kernel = np.array([1, 2, 3, 2, 1], dtype=float)
    kernel /= kernel.sum()
    smooth = np.convolve(counts, kernel, mode="same")

    peaks = []
    for i in range(1, len(smooth) - 1):
        if smooth[i] > smooth[i - 1] and smooth[i] >= smooth[i + 1]:
            peaks.append(i)

    if not peaks:
        return np.array([])

    min_height = 0.20 * np.max(smooth)
    peaks = [i for i in peaks if smooth[i] >= min_height]
    if not peaks:
        return np.array([])

    return np.sort(centers[peaks])


def kde_curve(gc_values: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    n = len(gc_values)
    if n <= 1:
        return np.zeros_like(x_grid)

    sd = float(np.std(gc_values, ddof=1))
    bw = max(1.06 * sd * (n ** (-1 / 5)), 0.01) if sd > 0 else 0.02

    z = (x_grid[:, None] - gc_values[None, :]) / bw
    density = np.mean(np.exp(-0.5 * z**2) / (bw * np.sqrt(2 * np.pi)), axis=1)
    return density


def main() -> None:
    OUTPUT_TABLE.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with gzip.open(INPUT_FASTA, "rt", encoding="utf-8") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            seq = str(rec.seq)
            rows.append(
                {
                    "accession": rec.id,
                    "length": len(seq),
                    "gc_raw": gc_fraction(seq),
                }
            )

    df = pd.DataFrame(rows).dropna(subset=["gc_raw"]).copy()
    gc = df["gc_raw"].to_numpy()

    n = len(gc)
    mean_gc = float(np.mean(gc))
    median_gc = float(np.median(gc))
    sd_gc = float(np.std(gc, ddof=1)) if n > 1 else 0.0

    modes = detect_modes(gc, bins=80)

    # Assign each transcript to a GC-content subpopulation based on mode midpoints
    if len(modes) >= 2:
        boundaries = (modes[:-1] + modes[1:]) / 2
        labels = [f"subpop_{i+1}" for i in range(len(modes))]
        df["gc_subpopulation"] = pd.cut(
            df["gc_raw"],
            bins=np.concatenate(([-np.inf], boundaries, [np.inf])),
            labels=labels,
            include_lowest=True,
        ).astype(str)
    else:
        df["gc_subpopulation"] = "subpop_1"

    # Save summary table
    out = df.sort_values("gc_raw", ascending=False).reset_index(drop=True).copy()
    out["gc_content"] = out["gc_raw"].map(lambda x: f"{x:.4f}")
    out = out[["accession", "length", "gc_content", "gc_subpopulation"]]
    out.to_csv(OUTPUT_TABLE, sep="\t", index=False)

    # Plot histogram + density + mean/median
    x_grid = np.linspace(0, 1, 500)
    density = kde_curve(gc, x_grid)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)  # 1600x900 px
    ax.hist(
        gc,
        bins=80,
        range=(0, 1),
        density=True,
        color="#74a9cf",
        edgecolor="white",
        linewidth=0.4,
        alpha=0.9,
        label="Histogram",
    )
    ax.plot(x_grid, density, color="#045a8d", linewidth=2.0, label="Density")
    ax.axvline(mean_gc, color="#d7301f", linestyle="--", linewidth=1.4, label=f"Mean: {mean_gc:.4f}")
    ax.axvline(median_gc, color="#238b45", linestyle="--", linewidth=1.4, label=f"Median: {median_gc:.4f}")

    ax.set_xlim(0, 1)
    ax.set_xlabel("GC content")
    ax.set_ylabel("Density")
    ax.set_title("Saccharomyces cerevisiae mRNA GC Content Distribution")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    caption = f"n={n} | mean={mean_gc:.4f} | median={median_gc:.4f} | sd={sd_gc:.4f}"
    fig.text(0.5, 0.01, caption, ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUTPUT_PLOT, dpi=200)
    plt.close(fig)

    print(f"Saved table: {OUTPUT_TABLE}")
    print(f"Saved plot:  {OUTPUT_PLOT}")
    print(caption)
    if len(modes) > 0:
        print("Detected GC mode(s):", ", ".join(f"{m:.4f}" for m in modes))
    else:
        print("Detected GC mode(s): none")


if __name__ == "__main__":
    main()
