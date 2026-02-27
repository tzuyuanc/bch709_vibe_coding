#!/usr/bin/env python3

import csv
import gzip
from collections import defaultdict
from pathlib import Path


CHROM_SIZES_PATH = Path("data/chrom.sizes")
GFF_PATH = Path("data/saccharomyces_cerevisiae.gff.gz")
OUTPUT_COUNTS_PATH = Path("results/chr_feature_counts.tsv")
DROPPED_SEQIDS_PATH = Path("results/dropped_seqids.txt")


def load_chrom_sizes(path: Path):
    chrom_order = []
    chrom_lengths = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            chrom, length_bp = line.split("\t")
            chrom_order.append(chrom)
            chrom_lengths[chrom] = int(length_bp)
    return chrom_order, chrom_lengths


def count_features(gff_path: Path, allowed_chroms):
    allowed_set = set(allowed_chroms)
    gene_counts = defaultdict(int)
    trna_counts = defaultdict(int)
    snorna_counts = defaultdict(int)
    exon_intervals = defaultdict(set)
    dropped_seqids = set()
    excluded_feature_lines = 0

    with gzip.open(gff_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid = fields[0]
            feature_type = fields[2]
            start = fields[3]
            end = fields[4]
            strand = fields[6]

            if seqid not in allowed_set:
                dropped_seqids.add(seqid)
                excluded_feature_lines += 1
                continue

            if feature_type == "gene":
                gene_counts[seqid] += 1
            elif feature_type == "exon":
                exon_intervals[seqid].add((start, end, strand))
            elif feature_type == "tRNA":
                trna_counts[seqid] += 1
            elif feature_type == "snoRNA":
                snorna_counts[seqid] += 1
            else:
                excluded_feature_lines += 1

    return (
        gene_counts,
        exon_intervals,
        trna_counts,
        snorna_counts,
        dropped_seqids,
        excluded_feature_lines,
    )


def write_outputs(
    chrom_order,
    chrom_lengths,
    gene_counts,
    exon_intervals,
    trna_counts,
    snorna_counts,
    dropped_seqids,
):
    OUTPUT_COUNTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for chrom in chrom_order:
        chrom_length_bp = chrom_lengths[chrom]
        per_mb_scale = chrom_length_bp / 1_000_000
        n_gene = gene_counts[chrom]
        n_exon_unique = len(exon_intervals[chrom])
        n_trna = trna_counts[chrom]
        n_snorna = snorna_counts[chrom]

        rows.append(
            {
                "chrom": chrom,
                "chrom_length_bp": chrom_length_bp,
                "n_gene": n_gene,
                "n_exon_unique": n_exon_unique,
                "n_tRNA": n_trna,
                "n_snoRNA": n_snorna,
                "gene_per_Mb": round(n_gene / per_mb_scale, 4),
                "exon_unique_per_Mb": round(n_exon_unique / per_mb_scale, 4),
                "tRNA_per_Mb": round(n_trna / per_mb_scale, 4),
                "snoRNA_per_Mb": round(n_snorna / per_mb_scale, 4),
            }
        )

    rows.sort(key=lambda row: row["gene_per_Mb"], reverse=True)

    with OUTPUT_COUNTS_PATH.open("w", newline="", encoding="utf-8") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "chrom_length_bp",
                "n_gene",
                "n_exon_unique",
                "n_tRNA",
                "n_snoRNA",
                "gene_per_Mb",
                "exon_unique_per_Mb",
                "tRNA_per_Mb",
                "snoRNA_per_Mb",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["chrom"],
                    row["chrom_length_bp"],
                    row["n_gene"],
                    row["n_exon_unique"],
                    row["n_tRNA"],
                    row["n_snoRNA"],
                    f'{row["gene_per_Mb"]:.4f}',
                    f'{row["exon_unique_per_Mb"]:.4f}',
                    f'{row["tRNA_per_Mb"]:.4f}',
                    f'{row["snoRNA_per_Mb"]:.4f}',
                ]
            )

    with DROPPED_SEQIDS_PATH.open("w", encoding="utf-8") as dropped_handle:
        for seqid in sorted(dropped_seqids):
            dropped_handle.write(f"{seqid}\n")


def print_top_rows(path: Path, n_rows: int = 5):
    with path.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle):
            if idx > n_rows:
                break
            print(line.rstrip("\n"))


def main():
    chrom_order, chrom_lengths = load_chrom_sizes(CHROM_SIZES_PATH)
    (
        gene_counts,
        exon_intervals,
        trna_counts,
        snorna_counts,
        dropped_seqids,
        excluded_feature_lines,
    ) = count_features(GFF_PATH, chrom_order)
    write_outputs(
        chrom_order,
        chrom_lengths,
        gene_counts,
        exon_intervals,
        trna_counts,
        snorna_counts,
        dropped_seqids,
    )
    print(f"dropped_seqids_count\t{len(dropped_seqids)}")
    print(f"excluded_feature_lines\t{excluded_feature_lines}")
    print_top_rows(OUTPUT_COUNTS_PATH, n_rows=5)


if __name__ == "__main__":
    main()
