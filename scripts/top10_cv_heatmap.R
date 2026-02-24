#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
condition_limit <- NA_integer_
if (length(args) >= 1 && nzchar(args[1])) {
  condition_limit <- suppressWarnings(as.integer(args[1]))
  if (is.na(condition_limit) || condition_limit <= 0) {
    stop("First argument must be a positive integer (number of conditions to use).")
  }
}

input_file <- "data/gasch2000.txt"
output_png <- "results/gene_top10_heatmap.png"
output_cv <- "results/top10_cv_table.tsv"

# Read safely: ignore comment lines and parse blanks as NA.
expr_raw <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  quote = "\"",
  comment.char = "#",
  na.strings = c("", "NA"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(expr_raw) < 4) {
  stop("Input table must have at least 4 columns (gene ID + metadata + conditions).")
}

# Column 1 is gene ID for labeling.
gene_ids <- as.character(expr_raw[[1]])

# Gasch matrix includes metadata columns before expression columns.
# Use columns 4..end as condition expression values.
all_condition_df <- expr_raw[, 4:ncol(expr_raw), drop = FALSE]
if (!is.na(condition_limit)) {
  if (condition_limit > ncol(all_condition_df)) {
    stop(sprintf(
      "Requested %d conditions but only %d are available.",
      condition_limit, ncol(all_condition_df)
    ))
  }
  condition_df <- all_condition_df[, seq_len(condition_limit), drop = FALSE]
  output_png <- sprintf("results/gene_top10_heatmap_first%d_conditions.png", condition_limit)
  output_cv <- sprintf("results/top10_cv_table_first%d_conditions.tsv", condition_limit)
} else {
  condition_df <- all_condition_df
}
condition_names <- colnames(condition_df)

# Coerce all condition columns to numeric (robust to occasional non-numeric tokens).
condition_num <- as.data.frame(lapply(condition_df, function(x) suppressWarnings(as.numeric(x))))
colnames(condition_num) <- condition_names

# Structure checks.
comment_lines <- sum(grepl("^\\s*#", readLines(input_file, warn = FALSE)))
gene_non_missing <- sum(!is.na(gene_ids) & nzchar(gene_ids))
gene_unique_non_missing <- length(unique(gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]))
gene_duplicate_non_missing <- gene_non_missing - gene_unique_non_missing
all_numeric_conditions <- all(vapply(condition_num, is.numeric, logical(1)))
na_exists <- any(is.na(condition_num))

cat("Structure check:\n")
cat(sprintf("- Rows: %d\n", nrow(expr_raw)))
cat(sprintf("- Columns: %d\n", ncol(expr_raw)))
cat(sprintf("- First column name: %s\n", colnames(expr_raw)[1]))
cat(sprintf("- Gene IDs non-missing: %d\n", gene_non_missing))
cat(sprintf("- Gene IDs unique (non-missing): %d\n", gene_unique_non_missing))
cat(sprintf("- Gene ID duplicates (non-missing): %d\n", gene_duplicate_non_missing))
cat(sprintf("- Condition columns available: %d\n", ncol(all_condition_df)))
cat(sprintf("- Condition columns used: %d\n", ncol(condition_num)))
if (!is.na(condition_limit)) {
  cat(sprintf("- Using first %d condition columns.\n", condition_limit))
}
cat(sprintf("- Condition columns numeric: %s\n", all_numeric_conditions))
cat(sprintf("- Any NA in condition data: %s\n", na_exists))
cat(sprintf("- Comment lines ignored via comment.char='#': %d\n", comment_lines))

# Decide whether log transform is needed.
mat <- as.matrix(condition_num)
neg_count <- sum(mat < 0, na.rm = TRUE)
value_range <- range(mat, na.rm = TRUE)
already_log <- neg_count > 0

if (already_log) {
  mat_log <- mat
  log_action <- "Detected negative values -> treating as already log-scale (no transform)."
} else {
  mat_log <- log2(mat + 1)
  log_action <- "No negative values detected -> applied log2(x + 1)."
}

cat(sprintf("- Value range before log decision: [%.4f, %.4f]\n", value_range[1], value_range[2]))
cat(sprintf("- Negative values count: %d\n", neg_count))
cat(sprintf("- Log handling: %s\n", log_action))

# CV per gene across all conditions, guarding mean==0.
row_mean <- rowMeans(mat_log, na.rm = TRUE)
row_sd <- apply(mat_log, 1, sd, na.rm = TRUE)
cv <- ifelse(row_mean == 0, NA_real_, row_sd / row_mean)

cv_df <- data.frame(
  gene = gene_ids,
  mean = row_mean,
  sd = row_sd,
  cv = cv,
  stringsAsFactors = FALSE
)

# Remove missing/empty gene IDs before ranking.
cv_df <- cv_df[!is.na(cv_df$gene) & nzchar(cv_df$gene), , drop = FALSE]

# Top 10 by CV descending.
ord <- order(cv_df$cv, decreasing = TRUE, na.last = NA)
top10 <- cv_df[ord, , drop = FALSE]
if (nrow(top10) > 10) top10 <- top10[1:10, , drop = FALSE]

write.table(top10, file = output_cv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("- Wrote top-10 CV table: %s\n", output_cv))

if (nrow(top10) == 0) {
  stop("No valid genes found for top-10 CV after filtering.")
}

# Build long format for heatmap with all conditions.
top_idx <- match(top10$gene, gene_ids)
top_mat <- mat_log[top_idx, , drop = FALSE]
rownames(top_mat) <- top10$gene

long_df <- as.data.frame(as.table(top_mat), stringsAsFactors = FALSE)
colnames(long_df) <- c("gene", "conditions", "expression")
long_df$gene <- factor(long_df$gene, levels = top10$gene)
long_df$conditions <- factor(long_df$conditions, levels = rev(colnames(top_mat)))

p <- ggplot(long_df, aes(x = gene, y = conditions, fill = expression)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_distiller(palette = "PuGn", direction = 1) +
  labs(
    title = "gene top 10",
    x = "gene",
    y = "conditions",
    fill = "expression"
  ) +
  theme_bw(base_family = "Times New Roman", base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black")
  )

ggsave(
  filename = output_png,
  plot = p,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat(sprintf("- Wrote heatmap: %s\n", output_png))
