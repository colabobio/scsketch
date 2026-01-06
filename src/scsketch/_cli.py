from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import importlib.metadata


IS_WINDOWS = sys.platform.startswith("win")


def _get_version() -> str:
    try:
        return importlib.metadata.version("scsketch")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _default_cache_dir() -> Path:
    if IS_WINDOWS:
        base = os.environ.get("LOCALAPPDATA") or os.environ.get("APPDATA")
        if base:
            return Path(base) / "scsketch"
        return Path.home() / "AppData" / "Local" / "scsketch"
    base = os.environ.get("XDG_CACHE_HOME")
    if base:
        return Path(base) / "scsketch"
    return Path.home() / ".cache" / "scsketch"


def _download(url: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url) as response, dst.open("wb") as f:
            f.write(response.read())
    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to download {url}: {e}") from e


def _download_demo_notebook(cache_dir: Path, ref_candidates: list[str]) -> Path:
    dst = cache_dir / "demo.ipynb"
    errors: list[str] = []

    for ref in ref_candidates:
        url = f"https://github.com/colabobio/scsketch/raw/{ref}/demo.ipynb"
        try:
            _download(url, dst)
            return dst
        except RuntimeError as e:
            errors.append(str(e))

    detail = "\n".join(f"- {e}" for e in errors)
    raise RuntimeError(
        "Could not download `demo.ipynb` from GitHub. Tried:\n" + detail
    )


def _run_jupyterlab_with_uv(notebook_path: Path, version: str) -> None:
    if shutil.which("uv") is None:
        raise RuntimeError(
            "Missing dependency: `uv` is not installed (https://github.com/astral-sh/uv)."
        )

    package_spec = "." if version == "unknown" else f"scsketch=={version}"

    command = [
        "uv",
        "tool",
        "run",
        "--python",
        "3.12",
        "--from",
        "jupyter-core",
        "--with",
        "jupyterlab",
        "--with",
        package_spec,
        "jupyter",
        "lab",
        str(notebook_path),
    ]

    if IS_WINDOWS:
        completed = subprocess.run(command)
        raise SystemExit(completed.returncode)
    os.execvp(command[0], command)


def _run_jupyterlab_in_current_env(notebook_path: Path) -> None:
    command = ["jupyter", "lab", str(notebook_path)]
    if IS_WINDOWS:
        completed = subprocess.run(command)
        raise SystemExit(completed.returncode)
    os.execvp(command[0], command)


def _demo(args: argparse.Namespace) -> None:
    version = _get_version()

    if args.local:
        repo_root = Path(__file__).resolve().parents[2]
        notebook_path = repo_root / "demo.ipynb"
        if not notebook_path.exists():
            raise RuntimeError(f"Local notebook not found: {notebook_path}")
    else:
        cache_dir = Path(args.cache_dir) if args.cache_dir else _default_cache_dir()
        ref_candidates = []
        if args.ref:
            ref_candidates.append(args.ref)
        if version != "unknown":
            ref_candidates.append(f"refs/tags/v{version}")
        ref_candidates.append("refs/heads/main")
        notebook_path = _download_demo_notebook(cache_dir=cache_dir, ref_candidates=ref_candidates)

    if args.no_uv:
        _run_jupyterlab_in_current_env(notebook_path)
        return

    try:
        _run_jupyterlab_with_uv(notebook_path, version=version)
    except RuntimeError as e:
        if args.fallback_to_env:
            _run_jupyterlab_in_current_env(notebook_path)
            return
        raise


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="scsketch")
    subparsers = parser.add_subparsers(dest="command", required=True)

    demo = subparsers.add_parser("demo", help="Open the demo notebook in JupyterLab")
    demo.add_argument(
        "--cache-dir",
        help="Cache directory for downloaded demo notebook (default: OS cache dir).",
    )
    demo.add_argument(
        "--ref",
        help="GitHub ref to download demo notebook from (e.g. refs/heads/main).",
    )
    demo.add_argument(
        "--local",
        action="store_true",
        help="Use local `demo.ipynb` from a source checkout (no download).",
    )
    demo.add_argument(
        "--no-uv",
        action="store_true",
        help="Run `jupyter lab` in the current environment instead of using `uv tool run`.",
    )
    demo.add_argument(
        "--fallback-to-env",
        action="store_true",
        help="If `uv` is missing, fall back to running `jupyter lab` in the current environment.",
    )
    demo.set_defaults(func=_demo)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
