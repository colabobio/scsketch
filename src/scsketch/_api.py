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
    List of ``{"name": str, "stId": str}`` dicts, empty on error.
    """
    url = f"{_REACTOME_BASE}/data/mapping/UniProt/{gene}/pathways?species=9606"
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        pathways = response.json()
        return [{"name": entry["displayName"], "stId": entry["stId"]} for entry in pathways]
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


def fetch_pathway_participants(pathway_id: str) -> list[str]:
    """Return UniProt identifiers of all proteins participating in a Reactome pathway.

    Parameters
    ----------
    pathway_id:
        Reactome stable identifier (e.g. ``"R-HSA-109581"``).

    Returns
    -------
    List of UniProt identifier strings, empty on error.
    """
    url = f"{_REACTOME_BASE}/data/participants/{pathway_id}"
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        data = response.json()
        return [
            ref["identifier"]
            for entry in data
            if "refEntities" in entry
            for ref in entry["refEntities"]
            if "identifier" in ref
        ]
    except requests.exceptions.RequestException as exc:
        logger.warning("Error fetching participants for pathway %s: %s", pathway_id, exc)
        return []


def gene_symbols_to_uniprot(gene_symbols: list[str]) -> list[str]:
    """Convert HGNC gene symbols to primary Swiss-Prot UniProt IDs via MyGene.info.

    Parameters
    ----------
    gene_symbols:
        List of gene symbols (e.g. ``["TP53", "BRCA1"]``).

    Returns
    -------
    List of UniProt accession strings (primary Swiss-Prot only), empty on error.
    """
    mapping: dict[str, str] = {}
    try:
        for gene in gene_symbols:
            url = f"https://mygene.info/v3/query?q={gene}&fields=uniprot.Swiss-Prot&species=human"
            response = requests.get(url, timeout=15)
            response.raise_for_status()
            hits = response.json().get("hits", [])
            for hit in hits:
                if "uniprot" in hit and isinstance(hit["uniprot"], dict):
                    primary = hit["uniprot"].get("Swiss-Prot")
                    if primary is not None:
                        mapping[gene] = primary[0] if isinstance(primary, list) else primary
                        break
    except requests.exceptions.RequestException as exc:
        logger.warning("Error fetching UniProt IDs: %s", exc)
    return list(mapping.values())
