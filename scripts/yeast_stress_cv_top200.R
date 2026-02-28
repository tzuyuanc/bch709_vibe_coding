#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(scales)
})

input_file <- "data/gasch2000.txt"
output_table <- "results/yeast_stress_cv_top200.tsv"
output_heatmap <- "results/yeast_stress_cv_top200_heatmap.png"

if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s", input_file))
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# Read TSV; preserve original column names and parse blanks as NA.
dt <- fread(
  input_file,
  sep = "\t",
  na.strings = c("", "NA"),
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(dt) < 5) {
  stop("Expected at least 5 columns: UID + metadata + expression conditions.")
}

# First column is gene_id/UID.
gene_id <- as.character(dt[[1]])

# Drop known metadata columns if present; always exclude first column.
drop_cols <- intersect(c(colnames(dt)[1], "NAME", "description", "GWEIGHT"), colnames(dt))
expr_df <- dt[, setdiff(colnames(dt), drop_cols), drop = FALSE]

if (ncol(expr_df) == 0) {
  stop("No condition columns available after removing metadata columns.")
}

# Coerce all condition columns to numeric.
expr_num <- as.data.frame(lapply(expr_df, function(x) suppressWarnings(as.numeric(x))))
colnames(expr_num) <- colnames(expr_df)
expr_mat <- as.matrix(expr_num)

# Filter genes: remove rows with all NA or zero variance.
all_na <- apply(expr_mat, 1, function(x) all(is.na(x)))
row_sd <- apply(expr_mat, 1, sd, na.rm = TRUE)
zero_var <- !is.na(row_sd) & row_sd == 0
keep <- !(all_na | zero_var)

expr_mat <- expr_mat[keep, , drop = FALSE]
gene_id <- gene_id[keep]

# Metrics across all conditions.
mean_expr <- rowMeans(expr_mat, na.rm = TRUE)
sd_expr <- apply(expr_mat, 1, sd, na.rm = TRUE)
cv <- sd_expr / abs(mean_expr)
cv[is.infinite(cv)] <- NA_real_

metrics <- data.frame(
  gene_id = gene_id,
  mean_expr = mean_expr,
  sd_expr = sd_expr,
  cv = cv,
  stringsAsFactors = FALSE
)

# Remove invalid rows after CV calculation and rank by CV descending.
metrics <- metrics[!is.na(metrics$gene_id) & nzchar(metrics$gene_id) & !is.na(metrics$cv), , drop = FALSE]
metrics <- metrics[order(metrics$cv, decreasing = TRUE), , drop = FALSE]
top_n <- min(200, nrow(metrics))

if (top_n == 0) {
  stop("No valid genes remain after filtering and CV calculation.")
}

top200 <- metrics[seq_len(top_n), , drop = FALSE]

# Round numeric output to 4 decimals and save.
top200_out <- top200
top200_out$mean_expr <- round(top200_out$mean_expr, 4)
top200_out$sd_expr <- round(top200_out$sd_expr, 4)
top200_out$cv <- round(top200_out$cv, 4)

write.table(
  top200_out[, c("gene_id", "mean_expr", "sd_expr", "cv"), drop = FALSE],
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Print top 10 summary table to console.
top10 <- head(top200_out[, c("gene_id", "mean_expr", "sd_expr", "cv"), drop = FALSE], 10)
print(top10, row.names = FALSE)

# Build heatmap matrix using top genes in CV-descending order.
idx <- match(top200$gene_id, gene_id)
heatmap_mat <- expr_mat[idx, , drop = FALSE]
rownames(heatmap_mat) <- top200$gene_id

# Preserve original condition order and fixed row order.
png(filename = output_heatmap, width = 1800, height = 1200, res = 200)
pheatmap(
  mat = heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 90,
  main = "Yeast stress response, CV top200 (Gasch et al. 2000)",
  color = viridisLite::viridis(100),
  border_color = NA,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)
dev.off()

cat(sprintf("\nSaved table: %s\n", output_table))
cat(sprintf("Saved heatmap: %s\n", output_heatmap))
