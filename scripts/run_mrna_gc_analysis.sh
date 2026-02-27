#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${HOME}/bch709_vibe_coding"
INPUT_FASTA="${PROJECT_DIR}/data/mrna.fa.gz"
OUTPUT_TABLE="${PROJECT_DIR}/results/mrna_metrics.tsv"
OUTPUT_PLOT="${PROJECT_DIR}/results/gc_content_distribution.png"
CONDA_ENV="bch709_vibe_coding"

if [[ ! -f "${INPUT_FASTA}" ]]; then
  echo "ERROR: Input FASTA not found: ${INPUT_FASTA}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUTPUT_TABLE}")"

# Use /tmp cache to avoid mamba cache permission issues in restricted environments.
XDG_CACHE_HOME=/tmp conda run -n "${CONDA_ENV}" python - <<'PY'
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

project_dir = Path.home() / "bch709_vibe_coding"
input_fasta = project_dir / "data" / "mrna.fa.gz"
output_table = project_dir / "results" / "mrna_metrics.tsv"
output_plot = project_dir / "results" / "gc_content_distribution.png"


def gc_content(seq: str) -> float:
    s = seq.upper()
    valid = sum(base in {"A", "C", "G", "T"} for base in s)
    if valid == 0:
        return np.nan
    gc = s.count("G") + s.count("C")
    return gc / valid


rows = []
with gzip.open(input_fasta, "rt", encoding="utf-8") as handle:
    for rec in SeqIO.parse(handle, "fasta"):
        seq = str(rec.seq)
        rows.append(
            {
                "accession": rec.id,
                "length": len(seq),
                "gc_content_raw": gc_content(seq),
            }
        )

df = pd.DataFrame(rows).dropna(subset=["gc_content_raw"]).copy()
gc = df["gc_content_raw"].to_numpy()

mean_gc = float(np.mean(gc))
median_gc = float(np.median(gc))
sd_gc = float(np.std(gc, ddof=1)) if len(gc) > 1 else 0.0

# Detect subpopulations using histogram-peak modes.
counts, edges = np.histogram(gc, bins=80, range=(0, 1))
centers = (edges[:-1] + edges[1:]) / 2
kernel = np.array([1, 2, 3, 2, 1], dtype=float)
kernel /= kernel.sum()
smooth = np.convolve(counts, kernel, mode="same")
peaks = []
for i in range(1, len(smooth) - 1):
    if smooth[i] > smooth[i - 1] and smooth[i] >= smooth[i + 1]:
        peaks.append(i)

if peaks:
    threshold = 0.20 * np.max(smooth)
    peaks = [i for i in peaks if smooth[i] >= threshold]
peak_centers = centers[peaks] if peaks else np.array([])

if len(peak_centers) >= 2:
    peak_centers = np.sort(peak_centers)
    boundaries = (peak_centers[:-1] + peak_centers[1:]) / 2
    labels = [f"subpop_{i+1}" for i in range(len(peak_centers))]
    df["gc_subpopulation"] = pd.cut(
        df["gc_content_raw"],
        bins=np.concatenate(([-np.inf], boundaries, [np.inf])),
        labels=labels,
        include_lowest=True,
    ).astype(str)
else:
    df["gc_subpopulation"] = "subpop_1"

df = df.sort_values("gc_content_raw", ascending=False).reset_index(drop=True)
df["gc_content"] = df["gc_content_raw"].map(lambda x: f"{x:.4f}")
df = df[["accession", "length", "gc_content", "gc_subpopulation"]]
df.to_csv(output_table, sep="\t", index=False)

# Build a KDE-like density curve (Gaussian kernel) without SciPy dependency.
x = np.linspace(0, 1, 500)
if len(gc) > 1 and sd_gc > 0:
    bw = max(1.06 * sd_gc * (len(gc) ** (-1 / 5)), 0.01)
else:
    bw = 0.02
density = np.mean(
    np.exp(-0.5 * ((x[:, None] - gc[None, :]) / bw) ** 2) / (bw * np.sqrt(2 * np.pi)),
    axis=1,
)

fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)  # 1600x900 px
hist = ax.hist(
    gc,
    bins=80,
    range=(0, 1),
    density=True,
    color="#74a9cf",
    edgecolor="white",
    linewidth=0.4,
    alpha=0.85,
    label="Histogram",
)
ax.plot(x, density, color="#045a8d", linewidth=2.0, label="Density curve")
ax.axvline(mean_gc, color="#d7301f", linestyle="--", linewidth=1.4, label=f"Mean: {mean_gc:.4f}")
ax.axvline(median_gc, color="#238b45", linestyle="--", linewidth=1.4, label=f"Median: {median_gc:.4f}")

ax.set_xlim(0, 1)
ax.set_xlabel("GC content")
ax.set_ylabel("Density")
ax.set_title("Saccharomyces cerevisiae mRNA GC Content Distribution")
caption = f"n={len(gc)} | mean={mean_gc:.4f} | median={median_gc:.4f} | sd={sd_gc:.4f}"
fig.text(0.5, 0.01, caption, ha="center", va="bottom", fontsize=9)
ax.grid(axis="y", alpha=0.2)
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(output_plot, dpi=200)
plt.close(fig)

print(f"Saved table: {output_table}")
print(f"Saved plot:  {output_plot}")
print(caption)
if len(peak_centers) > 0:
    print("Detected subpopulation modes:", ", ".join(f"{m:.4f}" for m in np.sort(peak_centers)))
else:
    print("Detected subpopulation modes: none")
PY

echo "Done."
