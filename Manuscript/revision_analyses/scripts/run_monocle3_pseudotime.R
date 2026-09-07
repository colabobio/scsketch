#!/usr/bin/env Rscript

# Clean Monocle3 pseudotime export for the scSketch revision analyses.
# The script is intentionally non-interactive. Set USE_PARTITIONS to TRUE for
# Monocle3's default partition-aware graph, or FALSE to learn one global graph
# across partitions. Embryo time is used only to orient/root the trajectory;
# correlation with embryo time is therefore reported as a sanity check, not as
# a fully independent validation.

suppressPackageStartupMessages({
  required <- c("monocle3", "igraph", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "), "\n",
      "Install Monocle3 with:\n",
      "  remotes::install_github('cole-trapnell-lab/monocle3')",
      call. = FALSE
    )
  }
  library(monocle3)
})

EXPRESSION_URL <- "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"
CELL_METADATA_URL <- "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"
GENE_ANNOTATION_URL <- "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"

USE_PARTITIONS <- isTRUE(as.logical(Sys.getenv("MONOCLE3_USE_PARTITIONS", "FALSE")))
MONOCLE3_MODE <- if (USE_PARTITIONS) "partitioned" else "global"

OUTPUT_DIR <- file.path(
  "Manuscript",
  "revision_analyses",
  "outputs",
  "trajectory_task_dataset_v1_participant",
  "trajectory_inference_methods"
)
OUTPUT_PREFIX <- file.path(OUTPUT_DIR, paste0("monocle3_", MONOCLE3_MODE))
OUTPUT_FILE <- paste0(OUTPUT_PREFIX, "_pseudotime.csv")
SUMMARY_FILE <- paste0(OUTPUT_PREFIX, "_pseudotime_summary.csv")
ROOTS_FILE <- paste0(OUTPUT_PREFIX, "_roots.csv")
PARTITION_SUMMARY_FILE <- paste0(OUTPUT_PREFIX, "_partition_summary.csv")
RDS_FILE <- paste0(OUTPUT_PREFIX, "_cds_ordered.rds")
QC_PREFIX <- paste0(OUTPUT_PREFIX, "_qc")

NUM_DIM <- 50
ROOT_QUANTILE <- 0.05
ROOT_MIN_CELLS <- 20
set.seed(0)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading Packer embryo data...")
expression_matrix <- readRDS(url(EXPRESSION_URL))
cell_metadata <- readRDS(url(CELL_METADATA_URL))
gene_annotation <- readRDS(url(GENE_ANNOTATION_URL))

cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

message("Running Monocle3 preprocessing...")
cds <- preprocess_cds(cds, num_dim = NUM_DIM)
cds <- align_cds(
  cds,
  alignment_group = "batch",
  residual_model_formula_str = paste(
    "~ bg.300.loading + bg.400.loading + bg.500.1.loading +",
    "bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading"
  )
)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = USE_PARTITIONS)

select_root_node_from_cells <- function(cds, root_candidate_cells) {
  closest_vertices <- as.matrix(
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  )
  closest_vertices <- closest_vertices[colnames(cds), , drop = FALSE]
  graph_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  nearest_vertex <- closest_vertices[root_candidate_cells, 1]
  root_vertex_index <- as.integer(names(which.max(table(nearest_vertex))))
  graph_nodes[root_vertex_index]
}

select_partition_roots <- function(
  cds,
  time_column = "embryo.time",
  root_quantile = 0.05,
  root_min_cells = 20
) {
  if (!time_column %in% colnames(colData(cds))) {
    stop("Missing required cell metadata column: ", time_column, call. = FALSE)
  }

  embryo_time <- suppressWarnings(as.numeric(as.character(colData(cds)[[time_column]])))
  names(embryo_time) <- colnames(cds)
  partition_ids <- partitions(cds)

  roots <- character()
  root_rows <- list()

  for (partition_id in sort(unique(as.character(partition_ids)))) {
    partition_cells <- names(partition_ids)[as.character(partition_ids) == partition_id]
    valid_cells <- partition_cells[
      is.finite(embryo_time[partition_cells])
    ]
    if (length(valid_cells) == 0) {
      next
    }

    # Use the early tail of each partition rather than a single absolute
    # earliest time point. This makes root selection less sensitive to one
    # sparsely sampled or outlying cell while preserving biological direction.
    cutoff <- as.numeric(stats::quantile(
      embryo_time[valid_cells],
      probs = root_quantile,
      na.rm = TRUE,
      names = FALSE
    ))
    root_candidate_cells <- valid_cells[embryo_time[valid_cells] <= cutoff]

    if (length(root_candidate_cells) < root_min_cells && length(valid_cells) > length(root_candidate_cells)) {
      ordered_cells <- valid_cells[order(embryo_time[valid_cells])]
      root_candidate_cells <- ordered_cells[
        seq_len(min(root_min_cells, length(ordered_cells)))
      ]
    }

    root_node <- select_root_node_from_cells(cds, root_candidate_cells)

    roots <- c(roots, root_node)
    root_rows[[length(root_rows) + 1]] <- data.frame(
      partition = as.character(partition_id),
      root_pr_node = root_node,
      root_quantile = root_quantile,
      root_time_cutoff = cutoff,
      n_valid_partition_cells = length(valid_cells),
      n_root_candidate_cells = length(root_candidate_cells),
      root_candidate_min_embryo_time = min(embryo_time[root_candidate_cells], na.rm = TRUE),
      root_candidate_max_embryo_time = max(embryo_time[root_candidate_cells], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  if (length(roots) == 0) {
    stop("Could not infer any Monocle3 root nodes.", call. = FALSE)
  }

  list(
    roots = unique(roots),
    root_table = do.call(rbind, root_rows)
  )
}

select_global_root <- function(
  cds,
  time_column = "embryo.time",
  root_quantile = 0.05,
  root_min_cells = 20
) {
  if (!time_column %in% colnames(colData(cds))) {
    stop("Missing required cell metadata column: ", time_column, call. = FALSE)
  }

  embryo_time <- suppressWarnings(as.numeric(as.character(colData(cds)[[time_column]])))
  names(embryo_time) <- colnames(cds)
  valid_cells <- names(embryo_time)[is.finite(embryo_time)]
  if (length(valid_cells) == 0) {
    stop("Could not infer a global Monocle3 root: embryo time is missing.", call. = FALSE)
  }

  # In global-graph mode, use a single early root so Monocle3 reports one
  # pseudotime scale across the connected graph.
  cutoff <- as.numeric(stats::quantile(
    embryo_time[valid_cells],
    probs = root_quantile,
    na.rm = TRUE,
    names = FALSE
  ))
  root_candidate_cells <- valid_cells[embryo_time[valid_cells] <= cutoff]

  if (length(root_candidate_cells) < root_min_cells && length(valid_cells) > length(root_candidate_cells)) {
    ordered_cells <- valid_cells[order(embryo_time[valid_cells])]
    root_candidate_cells <- ordered_cells[
      seq_len(min(root_min_cells, length(ordered_cells)))
    ]
  }

  root_node <- select_root_node_from_cells(cds, root_candidate_cells)
  root_table <- data.frame(
    graph_mode = "global",
    root_pr_node = root_node,
    root_quantile = root_quantile,
    root_time_cutoff = cutoff,
    n_valid_cells = length(valid_cells),
    n_root_candidate_cells = length(root_candidate_cells),
    root_candidate_min_embryo_time = min(embryo_time[root_candidate_cells], na.rm = TRUE),
    root_candidate_max_embryo_time = max(embryo_time[root_candidate_cells], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    roots = root_node,
    root_table = root_table
  )
}

root_info <- if (USE_PARTITIONS) {
  select_partition_roots(
    cds,
    root_quantile = ROOT_QUANTILE,
    root_min_cells = ROOT_MIN_CELLS
  )
} else {
  select_global_root(
    cds,
    root_quantile = ROOT_QUANTILE,
    root_min_cells = ROOT_MIN_CELLS
  )
}

message(
  "Ordering cells in ", MONOCLE3_MODE,
  " mode with root principal graph node(s): ",
  paste(root_info$roots, collapse = ", ")
)
cds <- order_cells(cds, root_pr_nodes = root_info$roots)

pseudotime_values <- pseudotime(cds)
embryo_time <- suppressWarnings(as.numeric(as.character(colData(cds)$embryo.time)))
names(embryo_time) <- colnames(cds)
embryo_time_ordered <- embryo_time[names(pseudotime_values)]
finite_mask <- is.finite(pseudotime_values) & is.finite(embryo_time_ordered)
partition_ids <- as.character(partitions(cds))
names(partition_ids) <- colnames(cds)

partition_summary <- do.call(rbind, lapply(sort(unique(partition_ids)), function(partition_id) {
  partition_cells <- names(partition_ids)[partition_ids == partition_id]
  partition_mask <- partition_cells[
    is.finite(pseudotime_values[partition_cells]) &
      is.finite(embryo_time[partition_cells])
  ]

  data.frame(
    graph_mode = MONOCLE3_MODE,
    partition = partition_id,
    n_cells = length(partition_cells),
    n_finite_pseudotime_and_embryo_time = length(partition_mask),
    spearman_vs_embryo_time = if (length(partition_mask) >= 3) {
      unname(cor(
        pseudotime_values[partition_mask],
        embryo_time[partition_mask],
        method = "spearman"
      ))
    } else {
      NA_real_
    },
    min_embryo_time = min(embryo_time[partition_cells], na.rm = TRUE),
    max_embryo_time = max(embryo_time[partition_cells], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

root_selection_note <- if (USE_PARTITIONS) {
  paste(
    "Roots were selected from early embryo-time cells per partition;",
    "global embryo-time correlation is diagnostic only."
  )
} else {
  paste(
    "A single root was selected from early embryo-time cells;",
    "use_partition = FALSE forces one graph across partitions."
  )
}

# Add explicit columns used by QC plots. Monocle3 can color by built-in
# pseudotime/partitions/clusters, but named columns make the saved plots and CSV
# provenance easier to audit.
colData(cds)$monocle_partition <- as.character(partitions(cds))
colData(cds)$monocle_cluster <- as.character(clusters(cds))
colData(cds)$monocle_pseudotime <- pseudotime_values[colnames(cds)]

save_monocle_qc_plot <- function(cds, color_by, label, filename_stub) {
  plot <- plot_cells(
    cds,
    color_cells_by = color_by,
    label_cell_groups = FALSE,
    label_branch_points = TRUE,
    label_roots = TRUE,
    label_leaves = FALSE,
    label_principal_points = TRUE,
    graph_label_size = 2,
    cell_size = 0.35
  ) + ggplot2::ggtitle(label)

  ggplot2::ggsave(
    paste0(QC_PREFIX, "_", filename_stub, ".png"),
    plot = plot,
    width = 6,
    height = 5,
    dpi = 300
  )
  ggplot2::ggsave(
    paste0(QC_PREFIX, "_", filename_stub, ".pdf"),
    plot = plot,
    width = 6,
    height = 5
  )
}

message("Writing Monocle3 QC plots...")
save_monocle_qc_plot(cds, "embryo.time", "Monocle3 UMAP colored by embryo time", "embryo_time")
save_monocle_qc_plot(cds, "monocle_partition", "Monocle3 UMAP colored by partition", "partition")
save_monocle_qc_plot(cds, "monocle_cluster", "Monocle3 UMAP colored by cluster", "cluster")
save_monocle_qc_plot(cds, "monocle_pseudotime", "Monocle3 UMAP colored by pseudotime", "pseudotime")

output <- data.frame(
  obs_name = names(pseudotime_values),
  method = "Monocle3",
  graph_mode = MONOCLE3_MODE,
  pseudotime = as.numeric(pseudotime_values),
  embryo_time = embryo_time_ordered,
  is_finite = is.finite(pseudotime_values),
  stringsAsFactors = FALSE
)

summary <- data.frame(
  method = "Monocle3",
  graph_mode = MONOCLE3_MODE,
  use_partitions = USE_PARTITIONS,
  n_cells = length(pseudotime_values),
  n_finite_pseudotime = sum(is.finite(pseudotime_values)),
  n_finite_pseudotime_and_embryo_time = sum(finite_mask),
  spearman_vs_embryo_time = unname(cor(
    pseudotime_values[finite_mask],
    embryo_time_ordered[finite_mask],
    method = "spearman"
  )),
  root_quantile = ROOT_QUANTILE,
  root_min_cells = ROOT_MIN_CELLS,
  root_pr_nodes = paste(root_info$roots, collapse = ";"),
  root_selection_note = root_selection_note,
  stringsAsFactors = FALSE
)

write.csv(output, OUTPUT_FILE, row.names = FALSE)
write.csv(summary, SUMMARY_FILE, row.names = FALSE)
write.csv(root_info$root_table, ROOTS_FILE, row.names = FALSE)
write.csv(partition_summary, PARTITION_SUMMARY_FILE, row.names = FALSE)
saveRDS(cds, RDS_FILE)

message("Wrote: ", OUTPUT_FILE)
message("Wrote: ", SUMMARY_FILE)
message("Wrote: ", ROOTS_FILE)
message("Wrote: ", PARTITION_SUMMARY_FILE)
message("Wrote: ", RDS_FILE)
message("Wrote Monocle3 QC plots with prefix: ", QC_PREFIX)
