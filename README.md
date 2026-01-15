# scSketch project

scSketch is an interactive exploration tool of single-cell embeddings (UMAP, tSNE, etc.) for Python notebooks. It is based on the [jupyter-scatter widget](https://jupyter-scatter.dev/) by [Fritz Lekschas](https://lekschas.de/) and it reimplements the earlier [SCIViwewer visualizer](https://github.com/colabobio/sciviewer).

It's currently provided as single notebook (scSketch.ipynb) with all the widget components and calculations in it. Plans involve to make it available as a package for easier installation.

## Usage

### Quick Start

The easiest way to try scSketch is with the built-in demo (no installation required):

```bash
uvx scsketch demo
```

This single command will automatically install scSketch and all dependencies, then launch the demo notebook.

Alternatively, if you've cloned the repository, you can run the demo notebook directly with juv:

```bash
git clone https://github.com/colabobio/scsketch.git
cd scsketch
uvx juv run demo.ipynb
```

### Installing as a Package

```bash
# Install in development mode
pip install -e .

# Or install from PyPI (when published)
pip install scsketch
```

Then use in any notebook:

```python
from scsketch import ScSketch

sketch = ScSketch(data=df, categorical_columns=cols)
sketch.show()
```

### Running the original notebook with juv

To run the original inline notebook, first install [juv](https://github.com/manzt/juv) and then call:

```bash
juv run scSketch.ipynb
```

### Creating a conda environment

Assuming that conda is already installed on the system, create an environment with all required dependencies:

```conda create --name scsketch --file requirements.txt --channel conda-forge python=3.12```

Once the environemnt is sucesfully created, activate it:

```conda activate scsketch```

and once in there, install the jupyter-scatter-scsketch package using pip:

```pip install jupyter-scatter-scsketch```

Once all of this is completed, launch JupyterLab in the current folder:

```jupyter lab```

and open the notebook.

## Development

### Publish a New Version

To bump the version use one of the following commands:

1. `uvx bump-my-version bump minor` (e.g., v0.1.0 → v0.2.0)
2. `uvx bump-my-version bump patch` (e.g., v0.1.0 → v0.1.1)
3. `uvx bump-my-version bump major` (e.g., v0.1.0 → v1.0.0)

Afterward do `git push --follow-tags`. Github actions will handle the rest.