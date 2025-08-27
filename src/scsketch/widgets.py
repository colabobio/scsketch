import requests
import traitlets
from anywidget import AnyWidget
from traitlets import List, Dict, Unicode, Int, Bool

#Additional UI Widgets - To visualize the results of scSketch, we need a few additional widgets:
#Directional Search Interactive Table Widget: a widget to visualiaze results of directional analysis of embedding and see what pathways the most upregulated and downregulated genes are a part of in Reactome.
#Label: a widget to display a selection of points
#Divider: a widget to visually add some clarity between groups of histograms

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
            print(f'❌ Error fetching pathways for {gene}: {e}')
            self.pathways = []

    @traitlets.observe('selected_pathway')
    def fetch_pathway_image(self, change):
        """Fetch the pathway image and participant proteins from Reactome API based on selected pathway ID"""
        pathway_id = change['new']
        if not pathway_id:
            self.pathway_image_url = ''
            return

        image_url = (
            f'https://reactome.org/ContentService/exporter/diagram/{pathway_id}.png'
        )
        self.pathway_image_url = image_url

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
                f'📌 Participant Proteins from Reactome API: {self.participant_proteins}'
            )
            print(f"📌 UniProt IDs for User's Genes: {uniprot_ids}")

            matched = list(
                set(self.participant_proteins).intersection(set(uniprot_ids))
            )
            self.matched_proteins = matched

            if matched:
                print(f'⚠️ Matched Proteins Found in Pathway {pathway_id}: {matched}')
            else:
                print(f'✅ No matched proteins found in Pathway {pathway_id}')

        except requests.exceptions.RequestException as e:
            print(f'❌ Error fetching participants for pathway {pathway_id}: {e}')
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

            print(f'✅ Gene Symbol to UniProt Mapping: {uniprot_mapping}')
            return list(uniprot_mapping.values())

        except requests.exceptions.RequestException as e:
            print(f'❌ Error fetching UniProt IDs: {e}')
            return []

#In the following we're going to create these widgets, which are all implemented with Trevor Manz's fantastic AnyWidget library.

# Directional Search Interactive Table Widget

class CorrelationTable(AnyWidget):
    _esm = """
    function render({ model, el }) {
      const container = document.createElement("div");
      const searchInput = document.createElement("input");
      searchInput.type = "text";
      searchInput.placeholder = "Search genes...";
      searchInput.style.width = "100%";
      searchInput.style.padding = "8px";
      searchInput.style.marginBottom = "8px";
      searchInput.style.boxSizing = "border-box";

      const table = document.createElement("table");
      table.classList.add("correlation-table");

      container.appendChild(searchInput);
      container.appendChild(table);
      el.appendChild(container);

      const pathwayTable = document.createElement("table");
      pathwayTable.classList.add("pathway-table");
      pathwayTable.style.display = "none";
      el.appendChild(pathwayTable);

      const pathwayImage = document.createElement("img");
      pathwayImage.style.display = "none";  
      pathwayImage.style.maxWidth = "100%";
      pathwayImage.alt = "Pathway Image";
      el.appendChild(pathwayImage);

      let rowsCache = [];
      const MAX_ROWS = 200; //minimum visible rows at a time

      const initializeTable = () => {
        const headerRow = document.createElement("tr");
        ["Gene", "R", "p"].forEach(col => {
          const th = document.createElement("th");
          th.textContent = col;
          headerRow.appendChild(th);
        });
        table.appendChild(headerRow);

        rowsCache = model.get("data").map(row => {
          const tr = document.createElement("tr");
          tr.dataset.gene = row["Gene"].toLowerCase();
          tr.style.cursor = "pointer";
          tr.onclick = () => {
            model.set("selected_gene", row["Gene"]);
            model.save_changes();
          };

          ["Gene", "R", "p"].forEach(col => {
            const td = document.createElement("td");
            td.textContent = row[col];
            tr.appendChild(td);
          });

          table.appendChild(tr);
          return tr; //caching the row
        });
      };

      initializeTable();

      let previousLength = 0;
      
      const updateTable = () => {
        const filterText = searchInput.value.toLowerCase();
        let visibleCount = 0;
        
        requestAnimationFrame(() => {
          rowsCache.forEach(row => {
            if (visibleCount < MAX_ROWS && row.dataset.gene.includes(filterText)) {
              row.style.display = "table-row";
              visibleCount++;
            } else {
              row.style.display = "none";
            }
          });
        });
      };

      function debounce(func, wait) {
        let timeout; 
        return (...args) => {
          clearTimeout(timeout);
          timeout = setTimeout(() => func.apply(this,args), wait);
        };
      }

      searchInput.addEventListener("input", debounce(() => {
        const currentLength = searchInput.value.length;
        debounce(updateTable, currentLength < previousLength ? 300 : 200)();
        previousLength = currentLength; 
      }, 50));

      model.on("change:pathways", () => {
        const pathways = model.get("pathways");
        pathwayTable.innerHTML = "";
        if (pathways.length > 0) {
          pathwayTable.style.display = "table";

          const headerRow = document.createElement("tr");
          ["Pathway"].forEach(header => {
            const th = document.createElement("th");
            th.textContent = header;
            headerRow.appendChild(th);
          });
          pathwayTable.appendChild(headerRow);

          pathways.forEach(pathway => {
            const row = document.createElement("tr");
            row.style.cursor = "pointer";
            row.onclick = () => {
              model.set("selected_pathway", pathway.stId);
              model.save_changes();
            };

            const td = document.createElement("td");
            td.textContent = pathway.name;
            row.appendChild(td);
            pathwayTable.appendChild(row);
          });

        } else {
          pathwayTable.style.display = "none";
        }
      });

      model.on("change:pathway_image_url", () => {
        const imageUrl = model.get("pathway_image_url");
        pathwayImage.src = imageUrl;
        pathwayImage.style.display = imageUrl ? "block" : "none";
      });

    }
    export default { render };
    """

    _css = """
    .correlation-table {
      width: 100%;
      border-collapse: collapse;
      margin-top: 10px;
    }
    .correlation-table th, .correlation-table td {
      border: 1px solid #ddd;
      padding: 8px;
      text-align: left;
    }
    .correlation-table th {
      background-color: #333;
      color: white;
    }
    .correlation-table tr:hover {
      background-color: #eee;
    }
    .pathway-table {
      width: 100%;
      border-collapse: collapse;
      margin-top: 10px;
    }
    .pathway-table th, .pathway-table td {
      border: 1px solid #ddd;
      padding: 8px;
      text-align: left;
    }
    .pathway-table th {
      background-color: #555;
      color: white;
    }
    .pathway-table tr:hover {
      background-color: #f2f2f2;
    }
    """

    data = List(Dict()).tag(sync=True)
    selected_gene = traitlets.Unicode("").tag(sync=True)
    pathways = traitlets.List([]).tag(sync=True)
    selected_pathway = traitlets.Unicode("").tag(sync=True)
    pathway_image_url = traitlets.Unicode("").tag(sync=True)
    participant_proteins = traitlets.List([]).tag(sync=True)
    matched_proteins = traitlets.List([]).tag(sync=True)

    def get_uniprot_ids(self, gene_symbols):
        """Convert gene symbols to UniProt IDs using MyGene.info API and ensure only primary IDs are used"""
        uniprot_mapping = {}

        try:
            for gene in gene_symbols:
                url = f"https://mygene.info/v3/query?q={gene}&fields=uniprot.Swiss-Prot&species=human"
                response = requests.get(url)
                response.raise_for_status()
                data = response.json().get("hits", [])

                if data:
                    for hit in data:
                        if "uniprot" in hit and isinstance(hit["uniprot"], dict):
                            if "Swiss-Prot" in hit["uniprot"]:
                                # Store only primary UniProt ID
                                primary_id = hit["uniprot"]["Swiss-Prot"]
                                if isinstance(primary_id, list):
                                    primary_id = primary_id[
                                        0
                                    ]  # Use the first one if multiple exist
                                uniprot_mapping[gene] = primary_id

            print(f"Gene Symbol to UniProt Mapping: {uniprot_mapping}")
            return list(uniprot_mapping.values())

        except requests.exceptions.RequestException as e:
            print(f"Error fetching UniProt IDs: {e}")
            return []

class PathwayTable(AnyWidget):
    _esm = """
    function render({ model, el }) {
      const table = document.createElement("table");
      table.classList.add("pathway-table");

      const update = () => {
        const pathways = model.get("data") || [];

        table.innerHTML = "";

        if (pathways.length === 0) {
          table.innerHTML = "<tr><td>No pathways available</td></tr>";
          return;
        }

        const headerRow = document.createElement("tr");
        ["Pathway"].forEach(col => {
          const th = document.createElement("th");
          th.textContent = col;
          headerRow.appendChild(th);
        });
        table.appendChild(headerRow);

        pathways.forEach(pathway => {
          const row = document.createElement("tr");
          row.style.cursor = "pointer";
          row.onclick = () => {
            model.set("selected_pathway", pathway.stId);
            model.save_changes();
          };

          const td = document.createElement("td");
          td.textContent = pathway.Pathway;
          row.appendChild(td);
          table.appendChild(row);
        });

        el.innerHTML = "";
        el.appendChild(table);
      };

      model.on("change:data", update);
      update();
    }
    export default { render };
    """

    _css = """
    .pathway-table {
      width: 100%;
      border-collapse: collapse;
      margin-top: 10px;
    }
    .pathway-table th, .pathway-table td {
      border: 1px solid #ddd;
      padding: 8px;
      text-align: left;
    }
    .pathway-table th {
      background-color: #555;
      color: white;
    }
    .pathway-table tr:hover {
      background-color: #f2f2f2;
    }
    """

    data = traitlets.List([]).tag(sync=True)
    selected_pathway = traitlets.Unicode("").tag(sync=True)


class InteractiveSVG(AnyWidget):
    _esm = """
    export function render({ model, el }) {
        el.style.position = 'relative';
        el.style.overflow = 'hidden';
        el.style.border = '1px solid #ddd';
        el.style.width = '100%';
        el.style.height = '800px';

        const container = document.createElement('div');
        container.style.width = '100%';
        container.style.height = '100%';
        container.style.overflow = 'auto';
        container.style.cursor = 'grab';
        container.style.position = 'relative';

        const img = document.createElement('img');
        img.style.transformOrigin = 'top left';
        img.style.width = '100%';
        img.style.height = 'auto';
        img.style.userSelect = 'none';

        let scale = 1;
        const scaleStep = 0.1;
        const minScale = 0.1;
        const maxScale = 10;

        const zoomInBtn = document.createElement('button');
        zoomInBtn.innerHTML = '+';
        zoomInBtn.style.position = 'absolute';
        zoomInBtn.style.top = '10px';
        zoomInBtn.style.right = '50px';
        zoomInBtn.style.zIndex = '10';

        const zoomOutBtn = document.createElement('button');
        zoomOutBtn.innerHTML = '−';
        zoomOutBtn.style.position = 'absolute';
        zoomOutBtn.style.top = '10px';
        zoomOutBtn.style.right = '10px';
        zoomOutBtn.style.zIndex = '10';

        function applyTransform() {
            img.style.transform = `scale(${scale})`;
        }

        zoomInBtn.onclick = () => {
            scale = Math.min(scale + scaleStep, maxScale);
            applyTransform();
        };

        zoomOutBtn.onclick = () => {
            scale = Math.max(scale - scaleStep, minScale);
            applyTransform();
        };

        container.appendChild(img);
        el.appendChild(container);
        el.appendChild(zoomInBtn);
        el.appendChild(zoomOutBtn);

        // Drag-to-pan
        let isDragging = false;
        let startX, startY, scrollLeft, scrollTop;

        container.addEventListener('mousedown', (e) => {
            isDragging = true;
            startX = e.pageX - container.offsetLeft;
            startY = e.pageY - container.offsetTop;
            scrollLeft = container.scrollLeft;
            scrollTop = container.scrollTop;
            container.style.cursor = 'grabbing';
            e.preventDefault();
        });

        container.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            const x = e.pageX - container.offsetLeft;
            const y = e.pageY - container.offsetTop;
            const walkX = x - startX;
            const walkY = y - startY;
            container.scrollLeft = scrollLeft - walkX;
            container.scrollTop = scrollTop - walkY;
        });

        window.addEventListener('mouseup', () => {
            isDragging = false;
            container.style.cursor = 'grab';
        });

        // Mouse wheel for zooming
        container.addEventListener('wheel', (e) => {
            e.preventDefault();
            const oldScale = scale;
            if (e.deltaY < 0) {
                scale = Math.min(scale + scaleStep, maxScale);
            } else {
                scale = Math.max(scale - scaleStep, minScale);
            }

            // Calculate zoom towards the mouse cursor
            const rect = container.getBoundingClientRect();
            const offsetX = (e.clientX - rect.left) + container.scrollLeft;
            const offsetY = (e.clientY - rect.top) + container.scrollTop;
            const dx = (offsetX / oldScale) * (scale - oldScale);
            const dy = (offsetY / oldScale) * (scale - oldScale);

            applyTransform();
            container.scrollLeft += dx;
            container.scrollTop += dy;
        });

        const update = () => {
            const svgContent = model.get('svg_content');
            img.src = `data:image/svg+xml;base64,${svgContent}`;
            scale = 1;
            applyTransform();
        };

        model.on('change:svg_content', update);
        update();
    }
    """

    svg_content = Unicode("").tag(sync=True)

#label widget - The next widget we're going to create is for representing a selection of points as a label. Nothing fancy here. The key thing we're going to use this for is to (a) tell you which points you have selected, (b) allow you to delete a selection, and (c) enable you to zoom to the selected points upon clicking on the label.


class Label(AnyWidget):
    _esm = """
    function render({ model, el }) {
      const label = document.createElement("div");
      label.classList.add(
        'jupyter-widgets',
        'jupyter-scatter-label'
      );
      label.tabIndex = 0;
      
      const update = () => {
        label.textContent = model.get('name');

        for (const [key, value] of Object.entries(model.get('style'))) {
          label.style[key] = value;
        }
      }
      
      model.on('change:name', update);
      model.on('change:style', update);

      update();

      const createFocusChanger = (value) => () => {
        model.set('focus', value);
        model.save_changes();
      }

      const focusHandler = createFocusChanger(true);
      const blurHandler = createFocusChanger(false);

      label.addEventListener('focus', focusHandler);
      label.addEventListener('blur', blurHandler);

      el.appendChild(label);

      const updateFocus = () => {
        if (model.get('focus')) {
          label.focus();
        }
      }
      
      model.on('change:focus', updateFocus);

      window.requestAnimationFrame(() => {
        updateFocus();
      });

      return () => {
        label.removeEventListener('focus', focusHandler);
        label.removeEventListener('blur', blurHandler);
      }
    }
    export default { render };
    """

    _css = """
    .jupyter-scatter-label {
      display: flex;
      align-items: center;
      width: 100%;
      height: var(--jp-widgets-inline-height);
      padding: var(--jp-widgets-input-padding) calc(var(--jp-widgets-input-padding)* 2);
      border-top-left-radius: var(--jp-border-radius);
      border-rop-right-radius: 0;
      border-bottom-left-radius: var(--jp-border-radius);
      border-bottom-right-radius: 0;
    }
    .jupyter-scatter-label:focus {
      font-weight: bold;
      outline: 1px solid var(--jp-widgets-input-focus-border-color);
      outline-offset: 1px;
    }
    """

    name = Unicode("").tag(sync=True)
    style = Dict({}).tag(sync=True)
    focus = Bool(False).tag(sync=True)

#Divider Widget - And finally, the technically most challenging (ahhh boring) widget for rendering a dividing horizontal line. Please don't waste time looking at the code as there's nothing to see here.

class Div(AnyWidget):
    _esm = """
    function render({ model, el }) {
      const div = document.createElement("div");
      div.classList.add(
        'jupyter-widgets',
        'jupyter-scatter-div'
      );
      
      const update = () => {
        for (const [key, value] of Object.entries(model.get('style'))) {
          div.style[key] = value;
        }
      }
      
      model.on('change', update);

      update();

      el.appendChild(div);
    }
    export default { render };
    """

    style = Dict({}).tag(sync=True)
