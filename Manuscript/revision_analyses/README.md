# scSketch Bioinformatics Revision Analyses

This folder contains reproducibility materials for the Bioinformatics major
revision of the scSketch Application Note.

## Dataset

`data/trajectory_task_dataset_v1_participant.h5ad` contains the C. elegans
embryo dataset used for the trajectory-inference comparison. The filename is
retained for path stability from earlier exploratory analyses.

The dataset is a 6,188-cell subset of the Packer and Zhu et al. C. elegans
embryo single-cell RNA-seq dataset:

Packer JS, Zhu Q, Huynh C, et al. A lineage-resolved molecular atlas of
C. elegans embryogenesis at single-cell resolution. Science. 2019;365(6459):
eaax1971. doi:10.1126/science.aax1971.

The source data are distributed in the Monocle3 C. elegans embryo tutorial as:

- `packer_embryo_expression.rds`
- `packer_embryo_colData.rds`
- `packer_embryo_rowData.rds`

For the scSketch revision analyses, these data were converted to AnnData format
and used to compare saved scSketch directional selections against pseudotime
estimates from PAGA/DPT, Slingshot, and Monocle3.

The `.h5ad` file is included in the archived release so Python/scverse users can
download the exact AnnData input used by the notebooks without rerunning the R
conversion step. The original RDS files remain the upstream source of record;
the AnnData file is a convenience conversion for reproducing the Python
scSketch, PAGA/DPT, and linked-view analyses.

## Main Reproducibility Workflow

1. `paga_scsketch_comparison_baseline.ipynb`
   Computes the Scanpy PAGA/DPT baseline and exports the scSketch session and
   PAGA/DPT outputs.

2. `scripts/run_slingshot_pseudotime.R`
   Computes Slingshot pseudotime from PCA coordinates and exports pseudotime
   and QC files.

3. `scripts/run_monocle3_pseudotime.R`
   Computes Monocle3 pseudotime using the global graph setting used in the
   primary manuscript comparison. The partition-aware run is retained as QC.

4. `trajectory_method_scsketch_comparison.ipynb`
   Joins saved scSketch selections to PAGA/DPT, Slingshot, and Monocle3
   pseudotime outputs by cell identifier and generates Supplementary Figure S3
   and Supplementary Table S2.

5. `umap_pca_linked_view_suppfig.ipynb`
   Generates the linked UMAP/PCA diagnostic figure used to contextualize
   embedding-dependent interpretation.

## Primary Outputs

The trajectory-method comparison outputs are written under:

`outputs/trajectory_task_dataset_v1_participant/trajectory_method_scsketch_comparison/`

The main reviewer-facing files are:

- `trajectory_method_revision_summary_primary.csv`
- `trajectory_method_revision_summary_figure.pdf`
- `trajectory_method_revision_summary_figure.png`
- `trajectory_method_selection_agreement.csv`
- `trajectory_method_gene_overlap_summary.csv`
- `trajectory_method_pseudotime_qc_summary.csv`
