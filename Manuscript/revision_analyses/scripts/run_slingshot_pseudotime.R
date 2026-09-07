#!/usr/bin/env Rscript

# Clean Slingshot pseudotime export for the scSketch revision analyses.
# Trajectory inference is run on PCA coordinates, not UMAP coordinates, so the
# comparison is not circular with scSketch's 2D embedding interaction.

suppressPackageStartupMessages({
  required <- c("SingleCellExperiment", "Matrix", "scrapper", "slingshot")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  library(SingleCellExperiment)
  library(Matrix)
  library(scrapper)
  library(slingshot)
})

EXPRESSION_URL <- "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"
CELL_METADATA_URL <- "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"
GENE_ANNOTATION_URL <- "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"

OUTPUT_DIR <- file.path(
  "Manuscript",
  "revision_analyses",
  "outputs",
  "trajectory_task_dataset_v1_participant",
  "trajectory_inference_methods"
)
OUTPUT_FILE <- file.path(OUTPUT_DIR, "slingshot_pseudotime.csv")
LINEAGE_FILE <- file.path(OUTPUT_DIR, "slingshot_pseudotime_by_lineage.csv")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "slingshot_pseudotime_summary.csv")
RDS_FILE <- file.path(OUTPUT_DIR, "slingshot_result.rds")
QC_PLOT_FILE <- file.path(OUTPUT_DIR, "slingshot_pca_qc.pdf")
QC_PNG_FILE <- file.path(OUTPUT_DIR, "slingshot_pca_qc.png")

N_TOP_GENES <- 2000
N_PCS <- 50
N_NEIGHBORS <- 30
LEIDEN_RESOLUTION <- 1.0
set.seed(0)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading Packer embryo data...")
expression_matrix <- readRDS(url(EXPRESSION_URL))
cell_metadata <- readRDS(url(CELL_METADATA_URL))
gene_annotation <- readRDS(url(GENE_ANNOTATION_URL))

sce <- SingleCellExperiment(
  assays = list(counts = expression_matrix),
  colData = cell_metadata,
  rowData = gene_annotation
)
stopifnot(inherits(sce, "SingleCellExperiment"))

message("Normalizing counts and computing PCA...")
sce <- normalizeRnaCounts.se(
  sce,
  log = TRUE,
  output.name = "logcounts",
  factor.name = "sizeFactor"
)
hvg_result <- chooseRnaHvgs.se(sce, top = N_TOP_GENES)
hvgs <- rownames(hvg_result)[rowData(hvg_result)$hvg]
sce <- runPca.se(
  sce,
  features = hvgs,
  number = N_PCS,
  block = colData(sce)$batch,
  assay.type = "logcounts",
  num.threads = 1
)

message("Clustering cells for Slingshot lineage inference...")
clustered <- clusterGraph.se(
  sce,
  num.neighbors = N_NEIGHBORS,
  reddim.type = "PCA",
  method = "leiden",
  resolution = LEIDEN_RESOLUTION
)
cluster_labels <- as.character(colData(clustered)$clusters)
if (length(unique(cluster_labels)) < 2) {
  stop("Slingshot requires at least two clusters.", call. = FALSE)
}

embryo_time <- suppressWarnings(as.numeric(as.character(colData(sce)$embryo.time)))
names(embryo_time) <- colnames(sce)
cluster_summary <- aggregate(
  embryo_time,
  by = list(cluster = cluster_labels),
  FUN = function(x) median(x[is.finite(x)], na.rm = TRUE)
)
names(cluster_summary)[2] <- "median_embryo_time"
cluster_summary <- cluster_summary[is.finite(cluster_summary$median_embryo_time), ]
if (nrow(cluster_summary) == 0) {
  stop("Could not choose a Slingshot start cluster from embryo.time.", call. = FALSE)
}
start_cluster <- cluster_summary$cluster[which.min(cluster_summary$median_embryo_time)]

message("Running Slingshot on PCA coordinates with start.clus = ", start_cluster)
pca_coords <- reducedDim(sce, "PCA")[, seq_len(min(N_PCS, ncol(reducedDim(sce, "PCA"))))]
rownames(pca_coords) <- colnames(sce)
sling <- slingshot(
  pca_coords,
  clusterLabels = cluster_labels,
  start.clus = start_cluster
)

pseudotime_by_lineage <- slingPseudotime(sling)
shared_pseudotime <- slingAvgPseudotime(sling)
if (is.null(names(shared_pseudotime))) {
  names(shared_pseudotime) <- rownames(pseudotime_by_lineage)
}

plot_continuous_pca <- function(coords, value, title, curves) {
  finite <- is.finite(value)
  palette <- grDevices::hcl.colors(100, "viridis")
  color_index <- rep(1L, length(value))
  if (any(finite)) {
    scaled <- (value[finite] - min(value[finite])) /
      max(diff(range(value[finite])), .Machine$double.eps)
    color_index[finite] <- pmax(1L, pmin(100L, floor(scaled * 99) + 1L))
  }
  point_colors <- rep("grey85", length(value))
  point_colors[finite] <- palette[color_index[finite]]

  plot(
    coords[, 1],
    coords[, 2],
    col = point_colors,
    pch = 16,
    cex = 0.35,
    xlab = "PC1",
    ylab = "PC2",
    main = title
  )
  draw_slingshot_curves(curves)
}

draw_slingshot_curves <- function(curves) {
  for (curve in curves) {
    curve_points <- curve$s[curve$ord, , drop = FALSE]
    lines(curve_points[, 1], curve_points[, 2], col = "black", lwd = 2)
  }
}

curves <- slingCurves(sling)
cluster_palette <- grDevices::hcl.colors(length(unique(cluster_labels)), "Dark 3")
cluster_colors <- setNames(cluster_palette, sort(unique(cluster_labels)))

write_qc_plot <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 12, height = 4)
  } else {
    grDevices::png(file, width = 2400, height = 800, res = 200)
  }
  old_par <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  plot(
    pca_coords[, 1],
    pca_coords[, 2],
    col = cluster_colors[cluster_labels],
    pch = 16,
    cex = 0.35,
    xlab = "PC1",
    ylab = "PC2",
    main = "Slingshot clusters"
  )
  draw_slingshot_curves(curves)
  plot_continuous_pca(
    pca_coords,
    embryo_time[rownames(pseudotime_by_lineage)],
    "Embryo time",
    curves
  )
  plot_continuous_pca(
    pca_coords,
    shared_pseudotime[rownames(pseudotime_by_lineage)],
    "Slingshot average pseudotime",
    curves
  )
}

message("Writing PCA QC plots...")
write_qc_plot(QC_PLOT_FILE, "pdf")
write_qc_plot(QC_PNG_FILE, "png")

lineage_output <- data.frame(
  obs_name = rownames(pseudotime_by_lineage),
  as.data.frame(pseudotime_by_lineage, check.names = FALSE),
  check.names = FALSE
)
output <- data.frame(
  obs_name = rownames(pseudotime_by_lineage),
  method = "Slingshot_PCA_rooted",
  pseudotime = as.numeric(shared_pseudotime),
  embryo_time = embryo_time[rownames(pseudotime_by_lineage)],
  cluster = cluster_labels,
  start_cluster = start_cluster,
  is_finite = is.finite(shared_pseudotime),
  stringsAsFactors = FALSE
)

finite_mask <- is.finite(shared_pseudotime) &
  is.finite(embryo_time[rownames(pseudotime_by_lineage)])
summary <- data.frame(
  method = "Slingshot_PCA_rooted",
  n_cells = length(shared_pseudotime),
  n_finite_pseudotime = sum(is.finite(shared_pseudotime)),
  n_finite_pseudotime_and_embryo_time = sum(finite_mask),
  spearman_vs_embryo_time = unname(cor(
    shared_pseudotime[finite_mask],
    embryo_time[rownames(pseudotime_by_lineage)][finite_mask],
    method = "spearman"
  )),
  n_clusters = length(unique(cluster_labels)),
  start_cluster = start_cluster,
  n_lineages = ncol(pseudotime_by_lineage),
  stringsAsFactors = FALSE
)

write.csv(output, OUTPUT_FILE, row.names = FALSE)
write.csv(lineage_output, LINEAGE_FILE, row.names = FALSE)
write.csv(summary, SUMMARY_FILE, row.names = FALSE)
write.csv(cluster_summary, file.path(OUTPUT_DIR, "slingshot_cluster_summary.csv"), row.names = FALSE)
saveRDS(sling, RDS_FILE)

message("Wrote: ", OUTPUT_FILE)
message("Wrote: ", LINEAGE_FILE)
message("Wrote: ", SUMMARY_FILE)
message("Wrote: ", RDS_FILE)
message("Wrote: ", QC_PLOT_FILE)
message("Wrote: ", QC_PNG_FILE)
