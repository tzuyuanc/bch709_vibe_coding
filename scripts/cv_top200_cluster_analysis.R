# Cluster 1: the expression of genes is inhibited by 37C to 25C shock and Heat Shock hs−2 the most, but is induced by Heat Shock hs−1 and heat shock 20 minutes.
# Cluster 2: the expression of genes is inhibited by 37C to 25C shoch and Heat Shock hs−1 strongly, but is induced by Heat Shock hs−2.
# Cluster 3: the expression of genes is induced by Heat Shock hs-2 and 37C to 25C shock strongly, but is inhibited by the other stresses (Heat Shock hs−1 is the most inhibitory).
# Cluster 4: the expression of genes is inhibited by Heat Shock hs−1, Heat Shock hs−2, and 37C to 25C shock, but is induced by the heat shock 20 minutes stress series.

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(viridisLite)
})

top200_file <- "results/yeast_stress_cv_top200.tsv"
gasch_file <- "data/gasch2000.txt"
output_heatmap <- "results/cv_top200_cluster_heatmap.pdf"
output_assign <- "results/cluster_assignment.tsv"
output_process <- "results/cluster_process_keywords.tsv"

if (!file.exists(top200_file)) {
  stop(sprintf("Missing input: %s", top200_file))
}
if (!file.exists(gasch_file)) {
  stop(sprintf("Missing input: %s", gasch_file))
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1) Load top-200 gene list.
top200_dt <- fread(top200_file, sep = "\t", data.table = FALSE, check.names = FALSE)
if (!("gene_id" %in% colnames(top200_dt))) {
  stop("results/yeast_stress_cv_top200.tsv must contain a 'gene_id' column.")
}
gene_list <- unique(as.character(top200_dt$gene_id))
gene_list <- gene_list[!is.na(gene_list) & nzchar(gene_list)]
if (length(gene_list) == 0) {
  stop("No valid gene_id values found in results/yeast_stress_cv_top200.tsv.")
}

# 2) Load Gasch matrix and keep numeric condition columns only.
gasch_dt <- fread(gasch_file, sep = "\t", data.table = FALSE, check.names = FALSE, na.strings = c("", "NA"))
if (ncol(gasch_dt) < 5) {
  stop("gasch2000.txt has too few columns.")
}

gasch_gene_id <- as.character(gasch_dt[[1]])
drop_cols <- intersect(c(colnames(gasch_dt)[1], "NAME", "description", "GWEIGHT"), colnames(gasch_dt))
cond_df <- gasch_dt[, setdiff(colnames(gasch_dt), drop_cols), drop = FALSE]
if (ncol(cond_df) < 30) {
  stop(sprintf("Need at least 30 condition columns; found %d.", ncol(cond_df)))
}

cond_num <- as.data.frame(lapply(cond_df, function(x) suppressWarnings(as.numeric(x))))
colnames(cond_num) <- colnames(cond_df)

# 3) Subset to first 30 condition columns.
cond30 <- cond_num[, seq_len(30), drop = FALSE]
mat30 <- as.matrix(cond30)

# 1 continued) Extract only top-200 genes by gene_id list.
row_idx <- match(gene_list, gasch_gene_id)
keep <- !is.na(row_idx)
if (!all(keep)) {
  warning(sprintf("Dropped %d gene_id values not found in gasch2000.txt.", sum(!keep)))
}
row_idx <- row_idx[keep]
gene_keep <- gene_list[keep]

if (length(row_idx) == 0) {
  stop("None of the top-200 genes were found in gasch2000.txt.")
}

mat_top <- mat30[row_idx, , drop = FALSE]
rownames(mat_top) <- gene_keep

# Drop rows with any NA or zero variance in first 30 conditions (required for dist/hclust).
row_all_na <- apply(mat_top, 1, function(x) all(is.na(x)))
row_sd <- apply(mat_top, 1, sd, na.rm = TRUE)
row_zero_sd <- !is.na(row_sd) & row_sd == 0
row_any_na <- apply(mat_top, 1, function(x) any(is.na(x)))
row_ok <- !(row_all_na | row_zero_sd | row_any_na)

if (!all(row_ok)) {
  warning(sprintf("Dropped %d rows (all-NA / zero-SD / contains-NA in first 30 conditions).", sum(!row_ok)))
}
mat_top <- mat_top[row_ok, , drop = FALSE]

if (nrow(mat_top) < 4) {
  stop(sprintf("Need at least 4 genes for k=4 clustering after filtering; found %d.", nrow(mat_top)))
}

# 4) Row-wise Z-score normalization.
row_mean <- rowMeans(mat_top)
row_sd2 <- apply(mat_top, 1, sd)
zmat <- sweep(mat_top, 1, row_mean, "-")
zmat <- sweep(zmat, 1, row_sd2, "/")

# 5) Hierarchical clustering: Euclidean + ward.D2.
d <- dist(zmat, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

# 6) cutree(k=4).
clusters <- cutree(hc, k = 4)
cluster_factor <- factor(clusters[rownames(zmat)], levels = 1:4)
annotation_row <- data.frame(cluster = cluster_factor)
rownames(annotation_row) <- rownames(zmat)

ann_colors <- list(
  cluster = c(
    "1" = "#1b9e77",
    "2" = "#d95f02",
    "3" = "#7570b3",
    "4" = "#e7298a"
  )
)

# Heatmap output.
pdf(output_heatmap, width = 8, height = 12)
pheatmap(
  mat = zmat,
  cluster_rows = hc,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 90,
  main = "Yeast stress response, CV top200 (Gasch et al. 2000)",
  color = viridisLite::viridis(100),
  border_color = NA,
  fontsize = 7,
  fontsize_row = 5,
  fontsize_col = 7
)
dev.off()

# Assignment table sorted by cluster ascending.
assign_df <- data.frame(
  gene_id = names(clusters),
  cluster = as.integer(clusters),
  stringsAsFactors = FALSE
)
assign_df <- assign_df[order(assign_df$cluster, assign_df$gene_id), , drop = FALSE]
write.table(assign_df, file = output_assign, sep = "\t", quote = FALSE, row.names = FALSE)

# Biological-process characterization using description keywords (if available).
desc_col <- intersect(c("description", "DESCRIPTION", "Description"), colnames(gasch_dt))
if (length(desc_col) > 0) {
  descriptions <- as.character(gasch_dt[[desc_col[1]]])
  desc_map <- data.frame(
    gene_id = gasch_gene_id,
    description = descriptions,
    stringsAsFactors = FALSE
  )
  ann <- merge(assign_df, desc_map, by = "gene_id", all.x = TRUE, sort = FALSE)

  # Simple text processing to highlight recurring biological-process terms.
  stop_words <- c(
    "protein", "putative", "required", "involved", "cell", "negative", "positive",
    "component", "subunit", "family", "activity", "process", "binding", "complex",
    "regulation", "response", "unknown", "function", "predicted", "probable",
    "yeast", "like", "domain", "containing", "associated", "factor"
  )

  extract_terms <- function(txt) {
    if (is.na(txt) || !nzchar(txt)) return(character(0))
    terms <- unlist(strsplit(tolower(gsub("[^a-z]+", " ", txt)), "\\s+"))
    terms <- terms[nchar(terms) >= 4]
    terms <- terms[!(terms %in% stop_words)]
    unique(terms)
  }

  process_rows <- list()
  for (k in sort(unique(ann$cluster))) {
    dsub <- ann$description[ann$cluster == k]
    words <- unlist(lapply(dsub, extract_terms))
    if (length(words) == 0) next
    freq <- sort(table(words), decreasing = TRUE)
    top_terms <- head(names(freq), 12)
    top_counts <- head(as.integer(freq), 12)
    process_rows[[length(process_rows) + 1]] <- data.frame(
      cluster = k,
      keyword = top_terms,
      count = top_counts,
      stringsAsFactors = FALSE
    )
  }

  if (length(process_rows) > 0) {
    process_df <- do.call(rbind, process_rows)
    process_df <- process_df[order(process_df$cluster, -process_df$count, process_df$keyword), , drop = FALSE]
    write.table(process_df, file = output_process, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("Saved process keywords: %s\n", output_process))
  } else {
    cat("No process keywords could be derived from description text.\n")
  }
} else {
  cat("Description column not found; skipped process keyword summary.\n")
}

cat(sprintf("Saved heatmap: %s\n", output_heatmap))
cat(sprintf("Saved assignments: %s\n", output_assign))
cat(sprintf("Genes clustered: %d\n", nrow(zmat)))
cat("Cluster sizes:\n")
print(table(assign_df$cluster))
