#!/usr/bin/env python3

import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

INPUT_FASTA = Path("data/mrna.fa.gz")
OUTPUT_TABLE = Path("results/mrna_metrics.tsv")
OUTPUT_PLOT = Path("results/gc_content_distribution.png")


def gc_content_fraction(seq: str) -> float:
    seq = seq.upper()
    valid = sum(base in {"A", "C", "G", "T"} for base in seq)
    if valid == 0:
        return np.nan
    gc = seq.count("G") + seq.count("C")
    return gc / valid


def load_mrna_metrics(fasta_path: Path) -> pd.DataFrame:
    records = []
    with gzip.open(fasta_path, "rt", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)
            records.append(
                {
                    "accession": record.id,
                    "length": len(seq),
                    "gc_content": gc_content_fraction(seq),
                }
            )

    df = pd.DataFrame(records)
    df = df.dropna(subset=["gc_content"]).copy()
    df = df.sort_values("gc_content", ascending=False).reset_index(drop=True)
    return df


def detect_gc_modes(gc_values: np.ndarray, bins: int = 80):
    counts, edges = np.histogram(gc_values, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2

    # Smooth counts to reduce noise, then detect local maxima as candidate modes.
    kernel = np.array([1, 2, 3, 2, 1], dtype=float)
    kernel = kernel / kernel.sum()
    smooth = np.convolve(counts, kernel, mode="same")

    peak_indices = []
    for i in range(1, len(smooth) - 1):
        if smooth[i] > smooth[i - 1] and smooth[i] >= smooth[i + 1]:
            peak_indices.append(i)

    if not peak_indices:
        return []

    max_height = smooth.max()
    strong_peaks = [idx for idx in peak_indices if smooth[idx] >= 0.2 * max_height]

    # Merge peaks that are too close in GC space.
    merged = []
    min_sep = (centers.max() - centers.min()) * 0.04
    for idx in strong_peaks:
        gc = centers[idx]
        if not merged:
            merged.append(idx)
        elif abs(gc - centers[merged[-1]]) < min_sep:
            if smooth[idx] > smooth[merged[-1]]:
                merged[-1] = idx
        else:
            merged.append(idx)

    return [centers[idx] for idx in merged]


def save_histogram(df: pd.DataFrame, out_path: Path):
    gc = df["gc_content"].to_numpy()
    modes = detect_gc_modes(gc)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)
    ax.hist(gc, bins=80, color="#2b8cbe", edgecolor="white", linewidth=0.4)

    mean_gc = float(np.mean(gc))
    median_gc = float(np.median(gc))
    ax.axvline(mean_gc, color="#d7301f", linestyle="--", linewidth=1.2, label=f"Mean: {mean_gc:.4f}")
    ax.axvline(median_gc, color="#31a354", linestyle="-.", linewidth=1.2, label=f"Median: {median_gc:.4f}")

    for i, m in enumerate(modes, start=1):
        ax.axvline(m, color="#756bb1", linestyle=":", linewidth=1.0)
        ax.text(m, ax.get_ylim()[1] * 0.92, f"mode{i}\n{m:.3f}", rotation=90, va="top", ha="center", fontsize=7)

    ax.set_title("Yeast mRNA GC Content Distribution")
    ax.set_xlabel("GC content")
    ax.set_ylabel("Transcript count")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.2)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)

    print(f"n_sequences\t{len(df)}")
    print(f"mean_gc\t{mean_gc:.4f}")
    print(f"median_gc\t{median_gc:.4f}")
    if modes:
        print(f"gc_modes\t{','.join(f'{m:.4f}' for m in modes)}")
    else:
        print("gc_modes\tnone_detected")


def main():
    df = load_mrna_metrics(INPUT_FASTA)

    OUTPUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    table_df = df.copy()
    table_df["gc_content"] = table_df["gc_content"].map(lambda x: f"{x:.4f}")
    table_df.to_csv(OUTPUT_TABLE, sep="\t", index=False)

    save_histogram(df, OUTPUT_PLOT)


if __name__ == "__main__":
    main()
