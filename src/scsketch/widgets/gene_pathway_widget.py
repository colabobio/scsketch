"""Gene/pathway selector widget.

Renders a gene dropdown; on gene selection fetches Reactome pathways via the
Python kernel (using :mod:`scsketch._api`) and populates a pathway dropdown.
All HTTP calls are handled in Python — the JS side is purely presentational.
"""

from pathlib import Path

import traitlets
from anywidget import AnyWidget
from traitlets import List, Unicode

from .._api import fetch_pathway_participants, gene_symbols_to_uniprot, fetch_pathways

_STATIC = Path(__file__).parent.parent / "static" / "widgets"


class GenePathwayWidget(AnyWidget):
    """A Jupyter widget to select genes, browse their Reactome pathways, and view pathway diagrams."""

    _esm = _STATIC / "gene_pathway_widget.js"
    _css = _STATIC / "gene_pathway_widget.css"

    # ── Synced state ─────────────────────────────────────────────────────────
    genes = List([]).tag(sync=True)
    selected_gene = Unicode("").tag(sync=True)
    pathways = List([]).tag(sync=True)
    selected_pathway = Unicode("").tag(sync=True)
    pathway_image_url = Unicode("").tag(sync=True)
    participant_proteins = List([]).tag(sync=True)
    matched_proteins = List([]).tag(sync=True)

    # ── Observers ────────────────────────────────────────────────────────────

    @traitlets.observe("selected_gene")
    def _on_gene_change(self, change):
        """Fetch Reactome pathways for the newly selected gene."""
        gene = change["new"]
        if not gene:
            self.pathways = []
            return
        raw = fetch_pathways(gene)
        self.pathways = raw

    @traitlets.observe("selected_pathway")
    def _on_pathway_change(self, change):
        """Update the pathway image URL and fetch participant proteins."""
        pathway_id = change["new"]
        if not pathway_id:
            self.pathway_image_url = ""
            self.participant_proteins = []
            self.matched_proteins = []
            return

        self.pathway_image_url = (
            f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.png"
        )

        participants = fetch_pathway_participants(pathway_id)
        self.participant_proteins = participants

        uniprot_ids = gene_symbols_to_uniprot(self.genes)
        self.matched_proteins = list(set(participants).intersection(uniprot_ids))

