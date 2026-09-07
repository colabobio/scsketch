"""External API calls — Reactome pathway data and SVG diagram retrieval.

All HTTP requests to third-party services live here so they can be mocked
in tests and replaced independently of UI code.
"""

from __future__ import annotations

import base64
import logging
from functools import lru_cache
from typing import Optional

import requests

logger = logging.getLogger(__name__)

# ── Reactome base URL ────────────────────────────────────────────────────────
_REACTOME_BASE = "https://reactome.org/ContentService"
_MYGENE_BASE = "https://mygene.info/v3"
_GENE_DESCRIPTION_FIELDS = "symbol,name,summary,entrezgene,ensembl.gene,taxid"
_REACTOME_SPECIES_ALIASES = {
    "human": "9606",
    "homo sapiens": "9606",
}


def _normalize_reactome_species(species: str | int | None = "human") -> str:
    """Return the Reactome species query value, preserving explicit overrides."""
    if species is None:
        return _REACTOME_SPECIES_ALIASES["human"]
    text = str(species).strip()
    if not text:
        return _REACTOME_SPECIES_ALIASES["human"]
    return _REACTOME_SPECIES_ALIASES.get(text.lower(), text)


@lru_cache(maxsize=512)
def fetch_gene_description(
    gene: str,
    species: str | int = "human",
) -> Optional[dict]:
    """Fetch a short gene annotation from MyGene.info.

    Parameters
    ----------
    gene:
        Gene symbol, Ensembl gene ID, or WormBase gene ID (e.g. ``"TP53"``,
        ``"ENSG00000141510"``, or ``"WBGene00010957"``).
    species:
        Species filter passed to MyGene.info query lookups. Defaults to human to
        match the existing Reactome integration.

    Returns
    -------
    Dict with selected MyGene fields, or ``None`` on error / no match.
    """
    gene = (gene or "").strip()
    if not gene:
        return None

    params = {"fields": _GENE_DESCRIPTION_FIELDS}
    try:
        gene_upper = gene.upper()
        if gene_upper.startswith("ENS"):
            response = requests.get(
                f"{_MYGENE_BASE}/gene/{gene}",
                params=params,
                timeout=15,
            )
            if response.status_code == 404:
                return None
            response.raise_for_status()
            payload = response.json()
        else:
            query = (
                f"wormbase:{gene}"
                if gene_upper.startswith("WBGENE")
                else f"symbol:{gene}"
            )
            response = requests.get(
                f"{_MYGENE_BASE}/query",
                params={
                    "q": query,
                    "fields": _GENE_DESCRIPTION_FIELDS,
                    "species": species,
                    "size": 1,
                },
                timeout=15,
            )
            response.raise_for_status()
            hits = response.json().get("hits", [])
            if not hits:
                return None
            payload = hits[0]

        if not isinstance(payload, dict):
            return None
        if not (payload.get("summary") or payload.get("name") or payload.get("symbol")):
            return None
        return payload
    except requests.exceptions.RequestException as exc:
        logger.warning("Error fetching MyGene description for %s: %s", gene, exc)
        return None


def fetch_pathways(gene: str, species: str | int = "human") -> list[dict]:
    """Fetch Reactome pathways for a gene symbol via UniProt mapping.

    Parameters
    ----------
    gene:
        Gene symbol (e.g. ``"TP53"``).
    species:
        Reactome species filter. Defaults to human. Numeric NCBI taxon IDs
        are passed through, and ``"human"`` / ``"Homo sapiens"`` are mapped
        to ``9606`` for backward-compatible behavior.

    Returns
    -------
    List of ``{"name": str, "stId": str}`` dicts, empty on error.
    """
    url = f"{_REACTOME_BASE}/data/mapping/UniProt/{gene}/pathways"
    try:
        response = requests.get(
            url,
            params={"species": _normalize_reactome_species(species)},
            timeout=15,
        )
        response.raise_for_status()
        pathways = response.json()
        return [
            {"name": entry["displayName"], "stId": entry["stId"]}
            for entry in pathways
        ]
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
    Base64-encoded UTF-8 string of the SVG content, or ``None`` on error /
    empty response.
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
        logger.warning(
            "Error fetching participants for pathway %s: %s",
            pathway_id,
            exc,
        )
        return []


def gene_symbols_to_uniprot(
    gene_symbols: list[str],
    species: str | int = "human",
) -> list[str]:
    """Convert gene symbols to primary Swiss-Prot UniProt IDs via MyGene.info.

    Parameters
    ----------
    gene_symbols:
        List of gene symbols (e.g. ``["TP53", "BRCA1"]``).
    species:
        Species filter passed to MyGene.info. Defaults to human.

    Returns
    -------
    List of UniProt accession strings (primary Swiss-Prot only), empty on error.
    """
    mapping: dict[str, str] = {}
    try:
        for gene in gene_symbols:
            response = requests.get(
                f"{_MYGENE_BASE}/query",
                params={
                    "q": gene,
                    "fields": "uniprot.Swiss-Prot",
                    "species": species,
                },
                timeout=15,
            )
            response.raise_for_status()
            hits = response.json().get("hits", [])
            for hit in hits:
                if "uniprot" in hit and isinstance(hit["uniprot"], dict):
                    primary = hit["uniprot"].get("Swiss-Prot")
                    if primary is not None:
                        mapping[gene] = (
                            primary[0] if isinstance(primary, list) else primary
                        )
                        break
    except requests.exceptions.RequestException as exc:
        logger.warning("Error fetching UniProt IDs: %s", exc)
    return list(mapping.values())
