"""External API calls — Reactome pathway data and SVG diagram retrieval.

All HTTP requests to third-party services live here so they can be mocked
in tests and replaced independently of UI code.
"""

from __future__ import annotations

import logging
import base64
from typing import Optional

import requests

logger = logging.getLogger(__name__)

# ── Reactome base URL ────────────────────────────────────────────────────────
_REACTOME_BASE = "https://reactome.org/ContentService"


def fetch_pathways(gene: str) -> list[dict]:
    """Fetch Reactome pathways for a human gene symbol via UniProt mapping.

    Parameters
    ----------
    gene:
        HGNC gene symbol (e.g. ``"TP53"``).

    Returns
    -------
    List of ``{"Pathway": str, "stId": str}`` dicts, empty on error.
    """
    url = f"{_REACTOME_BASE}/data/mapping/UniProt/{gene}/pathways?species=9606"
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        pathways = response.json()
        return [{"Pathway": entry["displayName"], "stId": entry["stId"]} for entry in pathways]
    except requests.exceptions.RequestException as exc:
        logger.warning("Error fetching Reactome pathways for %s: %s", gene, exc)
        return []


def fetch_pathway_svg(pathway_id: str) -> Optional[str]:
    """Fetch a Reactome pathway diagram as a base64-encoded SVG string.

    Parameters
    ----------
    pathway_id:
        Reactome stable identifier (e.g. ``"R-HSA-109581"``).

    Returns
    -------
    Base64-encoded UTF-8 string of the SVG content, or ``None`` on error / empty response.
    """
    svg_url = f"{_REACTOME_BASE}/exporter/diagram/{pathway_id}.svg"
    try:
        response = requests.get(svg_url, timeout=15)
        response.raise_for_status()
        svg_text = response.text.strip()
        if len(svg_text) < 50:
            logger.warning("Empty SVG returned from Reactome for %s", pathway_id)
            return None
        return base64.b64encode(svg_text.encode("utf-8")).decode("utf-8")
    except requests.exceptions.RequestException as exc:
        logger.warning("Error fetching SVG diagram for %s: %s", pathway_id, exc)
        return None
