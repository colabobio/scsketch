# scSketch project

scSketch is an interactive exploration tool of single-cell embeddings (UMAP, tSNE, etc.) for Python notebooks. It is based on the [jupyter-scatter widget](https://jupyter-scatter.dev/) by [Fritz Lekschas](https://lekschas.de/) and it reimplements the earlier [SCIViwewer visualizer](https://github.com/colabobio/sciviewer).

It's currently provided as single notebook (scSketch.ipynb) with all the widget components and calculations in it. Plans involve to make it available as a package for easier installation.

## Usage

The scSketch.ipynb requires a number of Python packages to be installed in order to run. There are two ways to install the dependencies, either with juv or creating a conda environment, as explained below

### Running the notebook with juv

To run the notebook, first install [juv](https://github.com/manzt/juv) and then call: 

```juv run scSketch.ipynb```

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
