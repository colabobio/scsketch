# scSketch

Interactive tooling for exploring single-cell embeddings in Jupyter (built on top of `jupyter-scatter` + `anywidget`).

## Install (Pip)

If you already have JupyterLab installed, this is all you need:

```sh
pip install scsketch
```

If you do not have JupyterLab yet, install it (or use the convenience extra):

```sh
pip install jupyterlab
# or: pip install "scsketch[lab]"
```

## Run

- Launch JupyterLab: `jupyter lab`
- Open `demo.ipynb` for an end-to-end example (it downloads example data from the internet).

## Install (Conda / Fallback)

If `pip install` fails on your system (common with scientific packages), use conda:

```sh
git clone https://github.com/colabobio/scsketch.git
cd scsketch

conda env create -f environment.yml
conda activate scsketch

jupyter lab demo.ipynb
```

## Troubleshooting

- Widgets don’t show up in JupyterLab:
  - Make sure `ipywidgets` and `jupyterlab_widgets` are installed in the same environment as `jupyter lab`.
  - Restart JupyterLab after installing dependencies.
  - If you’re on JupyterLab 3, you may need to run `jupyter lab build` once.

## Development

With [uv](https://github.com/astral-sh/uv):

```sh
uv run jupyter lab demo.ipynb
```

Or with editable installs:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
jupyter lab demo.ipynb
```
