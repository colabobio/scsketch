# Performance retest

This folder contains a small benchmark harness to retest scSketch performance on a fixed set of datasets and selection sizes.

## What it measures

Compute-only wall time (no UI / no network):

- **Directional Search**
  - `directional_corr_p_s`: correlation + p-values (`test_direction`)
  - `directional_lord_s`: LORD update (`lord_test`) timed separately
- **Differential Expression**
  - `de_global_stats_s`: one-time global sums/sumsq over `adata.X` (mirrors scSketch caching)
  - `de_compute_s`: per-selection Welch t-test from summary stats

## Datasets and selection sizes

Defaults are hard-coded in `benchmarks/perf_retest.py`:

- `pbmc3k.h5ad`: `[500, 1000, 2000]`
- `nygc_pbmc_161k_lite.h5ad`: `[20661, 25531, 50000]`
- `1M_20260121.h5ad`: `[25531, 143672, 376364]`

These files can be downloaded from:

* pbmc3k.h5ad (PBMC 3K): https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad
* nygc_pbmc_161k_lite.h5ad (161k): https://cellxgene.cziscience.com/collections/b0cf0afa-ec40-4d65-b570-ed4ceacc6813
* 1M_20260121.h5ad (1.2M): https://datasets.cellxgene.cziscience.com/c3d1a5e6-780b-4fe9-a39b-1864f927e87b.h5ad

Rename the 2nd and 3rd to nygc_pbmc_161k_lite.h5ad and 1M_20260121.h5ad respectively and move them under the root folder of this repo.

Selections are generated deterministically from the embedding by picking cells close to the principal axis and ordering them along that axis (to mimic a directional brush).

## Run

From the repo root:

```bash
uv run python benchmarks/perf_retest.py --out benchmarks/perf_results.csv
```

If you only want one dataset:

```bash
uv run python benchmarks/perf_retest.py --datasets pbmc161k --out benchmarks/perf_results.csv
```

## Optional: benchmark a real saved selection

You can export the indices from a scSketch notebook session and benchmark that exact selection.

In a notebook cell (after you have a saved selection):

```python
import numpy as np

sel = np.asarray(sketch.active_selection.points, dtype=int)
np.savez(
    "benchmarks/selection_pbmc161k_sel1.npz",
    dataset="pbmc161k",
    name="sel1",
    indices=sel,
)
```

Then benchmark it:

```bash
uv run python benchmarks/perf_retest.py --datasets pbmc161k --selection-npz benchmarks/selection_pbmc161k_sel1.npz --out benchmarks/perf_results.csv
```

To run *only* the exported selection (skip synthetic sizes), add:

```bash
uv run python benchmarks/perf_retest.py --datasets pbmc161k --selection-npz benchmarks/selection_pbmc161k_sel1.npz --no-synthetic --out benchmarks/perf_results.csv
```
