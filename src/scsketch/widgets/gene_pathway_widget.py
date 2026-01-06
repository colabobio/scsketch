import requests
import traitlets
from anywidget import AnyWidget


class GenePathwayWidget(AnyWidget):
    """A Jupyter Anywidget to select genes, view their pathways, and display pathway images."""

    _esm = """
    function render({ model, el }) {
        const geneDropdown = document.createElement("select");
        model.get("genes").forEach(gene => {
            const option = document.createElement("option");  
            option.value = gene;
            option.textContent = gene;
            geneDropdown.appendChild(option);  
        });
        el.appendChild(geneDropdown);

        const pathwayDropdown = document.createElement("select");
        pathwayDropdown.style.display = "none"; 
        el.appendChild(pathwayDropdown);

        const pathwayImage = document.createElement("img");
        pathwayImage.style.display = "none";  
        pathwayImage.style.maxWidth = "100%"; 
        pathwayImage.alt = "Pathway Image";
        el.appendChild(pathwayImage);

        geneDropdown.addEventListener("change", () => {
            const selectedGene = geneDropdown.value;
            model.set("selected_gene", selectedGene);
            model.save_changes();
        });

        pathwayDropdown.addEventListener("change", () => {
            const selectedPathwayId = pathwayDropdown.value;  
            model.set("selected_pathway", selectedPathwayId);
            model.save_changes();  
        });

        model.on("change:pathways", () => {
            const pathways = model.get("pathways");
            pathwayDropdown.innerHTML = ""; 
            if (pathways.length > 0) {
                pathwayDropdown.style.display = "block"; 
                pathways.forEach(pathway => {
                    const option = document.createElement("option");
                    option.value = pathway.stId;  
                    option.textContent = pathway.name;
                    pathwayDropdown.appendChild(option);
                });
            } else {
                pathwayDropdown.style.display = "none"; 
            }
        });

        model.on("change:pathway_image_url", () => {
            const imageUrl = model.get("pathway_image_url");
            if (imageUrl) {
                pathwayImage.src = imageUrl;
                pathwayImage.style.display = "block"; 
            } else {
                pathwayImage.style.display = "none"; 
            }
        });
    }
    export default { render };
    """

    # List of genes
    genes = traitlets.List([]).tag(sync=True)
    selected_gene = traitlets.Unicode('').tag(sync=True)
    pathways = traitlets.List([]).tag(sync=True)
    selected_pathway = traitlets.Unicode('').tag(sync=True)
    pathway_image_url = traitlets.Unicode('').tag(sync=True)
    participant_proteins = traitlets.List([]).tag(sync=True)
    matched_proteins = traitlets.List([]).tag(sync=True)

    @traitlets.observe('selected_gene')
    def fetch_pathways(self, change):
        """Fetch pathways for the selected gene from Reactome API"""
        gene = change['new']
        if not gene:
            self.pathways = []
            return

        try:
            url = f'https://reactome.org/ContentService/data/mapping/UniProt/{gene}/pathways?species=9606'
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()

            pathways = [
                {'name': entry['displayName'], 'stId': entry['stId']}
                for entry in data
                if 'stId' in entry
            ]
            self.pathways = pathways
        except requests.exceptions.RequestException as e:
            print(f'Error fetching pathways for {gene}: {e}')
            self.pathways = []

    @traitlets.observe('selected_pathway')
    def fetch_pathway_image(self, change):
        """Fetch and display a pathway diagram image based on selected pathway ID"""
        pathway_id = change['new']
        if not pathway_id:
            self.pathway_image_url = ''
            return

        self.pathway_image_url = (
            f'https://reactome.org/ContentService/exporter/diagram/{pathway_id}.png'
        )

    @traitlets.observe('selected_pathway')
    def fetch_participants(self, change):
        """Fetch participants in the selected pathway and find overlaps with user's genes."""
        pathway_id = change['new']
        if not pathway_id:
            self.participant_proteins = []
            self.matched_proteins = []
            return

        try:
            participants_url = (
                f'https://reactome.org/ContentService/data/participants/{pathway_id}'
            )
            response = requests.get(participants_url)
            response.raise_for_status()
            participants_data = response.json()

            self.participant_proteins = [
                ref['identifier']
                for entry in participants_data
                if 'refEntities' in entry
                for ref in entry['refEntities']
                if 'identifier' in ref
            ]

            uniprot_ids = self.get_uniprot_ids(self.genes)

            print(
                f'Participant Proteins from Reactome API: {self.participant_proteins}'
            )
            print(f"UniProt IDs for User's Genes: {uniprot_ids}")

            matched = list(
                set(self.participant_proteins).intersection(set(uniprot_ids))
            )
            self.matched_proteins = matched

            if matched:
                print(f'Matched Proteins Found in Pathway {pathway_id}: {matched}')
            else:
                print(f'No matched proteins found in Pathway {pathway_id}')

        except requests.exceptions.RequestException as e:
            print(f'Error fetching participants for pathway {pathway_id}: {e}')
            self.participant_proteins = []

    def get_uniprot_ids(self, gene_symbols):
        """Convert gene symbols to UniProt IDs using MyGene.info API and ensure only primary IDs are used"""
        uniprot_mapping = {}

        try:
            for gene in gene_symbols:
                url = f'https://mygene.info/v3/query?q={gene}&fields=uniprot.Swiss-Prot&species=human'
                response = requests.get(url)
                response.raise_for_status()
                data = response.json().get('hits', [])

                if data:
                    for hit in data:
                        if 'uniprot' in hit and isinstance(hit['uniprot'], dict):
                            if 'Swiss-Prot' in hit['uniprot']:
                                # Store only primary UniProt ID
                                primary_id = hit['uniprot']['Swiss-Prot']
                                if isinstance(primary_id, list):
                                    primary_id = primary_id[
                                        0
                                    ]  # Use the first one if multiple exist
                                uniprot_mapping[gene] = primary_id

            print(f'Gene Symbol to UniProt Mapping: {uniprot_mapping}')
            return list(uniprot_mapping.values())

        except requests.exceptions.RequestException as e:
            print(f'Error fetching UniProt IDs: {e}')
            return []
