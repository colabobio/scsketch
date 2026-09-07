#!/usr/bin/env Rscript

if (!requireNamespace("monocle3", quietly = TRUE)) {
  stop(
    paste(
      "The R package 'monocle3' is required to export the tutorial UMAP.",
      "Install it in R before running this script:",
      "remotes::install_github('cole-trapnell-lab/monocle3')",
      sep = "\n"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(monocle3))

base_dir <- normalizePath(
  file.path(getwd(), "Manuscript", "revision_analyses"),
  mustWork = FALSE
)
data_dir <- file.path(base_dir, "data", "monocle_rds")
output_dir <- file.path(base_dir, "outputs", "monocle_umap")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_base_url <- "https://depts.washington.edu/trapnell-lab/software/monocle3/celegans/data"
files <- c(
  expression = "packer_embryo_expression.rds",
  coldata = "packer_embryo_colData.rds",
  rowdata = "packer_embryo_rowData.rds"
)

for (filename in files) {
  destination <- file.path(data_dir, filename)
  if (!file.exists(destination)) {
    download.file(
      url = paste0(data_base_url, "/", filename),
      destfile = destination,
      mode = "wb",
      quiet = FALSE
    )
  }
}

expression_matrix <- readRDS(file.path(data_dir, files[["expression"]]))
cell_metadata <- readRDS(file.path(data_dir, files[["coldata"]]))
gene_annotation <- readRDS(file.path(data_dir, files[["rowdata"]]))

cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(
  cds,
  alignment_group = "batch",
  residual_model_formula_str = paste(
    "~ bg.300.loading + bg.400.loading + bg.500.1.loading +",
    "bg.500.2.loading + bg.r17.loading + bg.b01.loading +",
    "bg.b02.loading"
  )
)

set.seed(0)
cds <- reduce_dimension(cds)

umap <- reducedDims(cds)[["UMAP"]]
coords <- data.frame(
  cell = rownames(umap),
  UMAP1 = umap[, 1],
  UMAP2 = umap[, 2],
  check.names = FALSE
)

output_path <- file.path(output_dir, "monocle_tutorial_umap.csv")
write.csv(coords, output_path, row.names = FALSE)
message("Wrote Monocle UMAP coordinates: ", output_path)
