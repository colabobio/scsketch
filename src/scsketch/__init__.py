import importlib.metadata
import pathlib

import anywidget
import traitlets

try:
    __version__ = importlib.metadata.version("scsketch")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


class Widget(anywidget.AnyWidget):
    _esm = pathlib.Path(__file__).parent / "static" / "widget.js"
    _css = pathlib.Path(__file__).parent / "static" / "widget.css"
    value = traitlets.Int(0).tag(sync=True)



#Additional UI Widgets - To visualize the results of scSketch, we need a few additional widgets:
#Directional Search Interactive Table Widget: a widget to visualiaze results of directional analysis of embedding and see what pathways the most upregulated and downregulated genes are a part of in Reactome.
#Label: a widget to display a selection of points
#Divider: a widget to visually add some clarity between groups of histograms


class GenePathwayWidget(anywidget.AnyWidget):
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
        pathwayImage.alt = "Pathway Image";ƒ
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
from anywidget import AnyWidget
from traitlets import List, Dict, Unicode, Int
# from dropdown_widget import GenePathwayWidget
import traitlets
import requests


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


from anywidget import AnyWidget
from traitlets import Unicode


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
from anywidget import AnyWidget
from traitlets import Bool, Dict, Unicode


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

from anywidget import AnyWidget
from traitlets import Bool, Dict, Unicode


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

#Widget Composition - Finally, we're going to instantiate the scatter plot and all the other widgets and link them using their traits. The output is the UI you've been waiting for :)
def view(adata):
    print("I am in view function")
    from IPython.display import display, HTML
    import numpy as np
    
    from collections import OrderedDict
    from dataclasses import dataclass, field
    
    # from dimbridge.predicate_engine import compute_predicate_sequence
    from IPython.display import HTML
    from ipywidgets import Checkbox, Dropdown, GridBox, HBox, Layout, IntText, Text, VBox
    from itertools import cycle
    from jscatter import Scatter, glasbey_light, link, okabe_ito, Line
    from jscatter.widgets import Button
    from numpy import histogram, isnan
    from matplotlib.colors import to_hex
    from scipy.spatial import ConvexHull
####################################################################################################################################
    import pandas as pd
    
    # UMAP coordinates
    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["x", "y"], index=adata.obs_names)
    
    # add metadata (I will pick columns I am interested in)
    metadata_cols = [
        "dpi",
        "strain",
        "percent.cmv.log10",
        "seurat_clusters",
        "virus.presence",
        "Infection_localminima",
        "UL123_define_infection",
        "Infection_state",
        "Infection_state_bkgd",
    ]
    metadata_df = adata.obs[metadata_cols]
    
    # Extract gene expression from adata.raw
    gene_exp_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        columns=adata.var_names,
        index=adata.obs_names,
    )
    
    # combine into a single dataframe
    df = pd.concat([umap_df, metadata_df, gene_exp_df], axis=1)
    df.head()
    
    # Ensure all column names are unique by standardizing capitalization
    # df.columns = df.columns.str.lower()
    df = df.loc[:, ~df.columns.duplicated()]
    
    
    categorical_cols = [
        "dpi",
        "strain",
        "seurat_clusters",
        "virus.presence",
        "Infection_localminima",
        "UL123_define_infection",
        "Infection_state",
        "Infection_state_bkgd",
    ]
    
    for col in categorical_cols:
        df[col] = df[col].astype(str)
    
    # Generate explicit color map for categorical variables
    from itertools import cycle
    from jscatter import glasbey_light
    
    categorical_color_maps = {}
    
    for col in categorical_cols:
        unique_categories = df[col].unique()
        categorical_color_maps[col] = dict(zip(unique_categories, cycle(glasbey_light)))
################################################################################################################################
    
    all_colors = okabe_ito.copy()
    available_colors = [color for color in all_colors]
    
    color_by = "seurat_clusters"
    
    df[color_by] = df[color_by].astype(str)
    
    non_protein_cols = ["x", "y"]
    protein_cols = [col for col in df.columns if col not in non_protein_cols]
    
    color_map = dict(zip(df[color_by].unique(), cycle(glasbey_light[1:])))
    color_map["Non-robust"] = (0.2, 0.2, 0.2, 1.0)
    
    continuous_color_maps = [
        ["#00dadb", "#da00db"],
        ["#00dadb", "#a994dc", "#da00db"],
        ["#00dadb", "#8faddc", "#bd77dc", "#da00db"],
        ["#00dadb", "#7eb9dc", "#a994dc", "#c567dc", "#da00db"],
        ["#00dadb", "#72c0db", "#9aa3dc", "#b583dc", "#ca5cdb", "#da00db"],
        ["#00dadb", "#69c4db", "#8faddc", "#a994dc", "#bd77dc", "#cd54db", "#da00db"],
        [
            "#00dadb",
            "#62c7db",
            "#86b4dc",
            "#9e9fdc",
            "#b288dc",
            "#c16edc",
            "#cf4ddb",
            "#da00db",
        ],
        [
            "#00dadb",
            "#5ccadb",
            "#7eb9dc",
            "#96a7dc",
            "#a994dc",
            "#b87fdc",
            "#c567dc",
            "#d048db",
            "#da00db",
        ],
        [
            "#00dadb",
            "#57ccdb",
            "#78bddc",
            "#8faddc",
            "#a19ddc",
            "#b08bdc",
            "#bd77dc",
            "#c861db",
            "#d144db",
            "#da00db",
        ],
    ]
    
    color_by_default = "seurat_clusters"  # initial default coloring
    
    scatter = Scatter(
        data=df,
        x="x",
        y="y",
        background_color="#111111",
        axes=False,
        height=720,
        color_by=color_by,
        color_map=color_map,
        tooltip=True,
        legend=True,
        tooltip_properties=[
            "dpi",
            "strain",
            "percent.cmv.log10",
            "seurat_clusters",
            "virus.presence",
            "Infection_localminima",
            "UL123_define_infection",
            "Infection_state",
            "Infection_state_bkgd",
        ],
        tooltip_histograms_size="large",
        legend_labels={key: key for key in df["seurat_clusters"].unique()},
    )
    
    # Remove non-robust cell populations from the histogram as it's uninteresting
    color_histogram = scatter.widget.color_histogram.copy()
    
    if "Non-robust" in scatter._color_categories:
        color_histogram[scatter._color_categories["Non-robust"]] = 0
    
    scatter.widget.color_histogram = color_histogram
    
    
    @dataclass
    class Selection:
        """Class for keeping track of a selection."""
    
        index: int
        name: str
        points: np.ndarray
        color: str
        lasso: Line
        hull: Line
    
    
    @dataclass
    class Selections:
        """Class for keeping track of selections."""
    
        selections: list[Selection] = field(default_factory=list)
    
        def all_points(self) -> np.ndarray:
            return np.unique(
                np.concatenate(
                    list(map(lambda selection: selection.points, self.selections))
                )
            )
    
        def all_hulls(self) -> list[Line]:
            return [s.hull for s in self.selections]
    
    
    @dataclass
    class Lasso:
        """Class for keeping track of the lasso polygon."""
    
        polygon: Line | None = None
    
    
    lasso = Lasso()
    selections = Selections()
    
    
    def update_annotations():
        lasso_polygon = [] if lasso.polygon is None else [lasso.polygon]
        scatter.annotations(selections.all_hulls() + lasso_polygon)
    
    
    def lasso_selection_polygon_change_handler(change):
        if change["new"] is None:
            lasso.polygon = None
        else:
            points = change["new"].tolist()
            points.append(points[0])
            lasso.polygon = Line(points, line_color=scatter.widget.color_selected)
        update_annotations()
    
    
    scatter.widget.observe(
        lasso_selection_polygon_change_handler, names=["lasso_selection_polygon"]
    )
    
    selection_name = Text(value="", placeholder="Select some points…", disabled=True)
    selection_name.layout.width = "100%"
    
    selection_add = Button(
        description="",
        tooltip="Save Selection",
        disabled=True,
        icon="plus",
        width=36,
        rounded=["top-right", "bottom-right"],
    )
    
    selection_subdivide = Checkbox(value=False, description="Subdivide", indent=False)
    
    selection_num_subdivisions = IntText(
        value=5,
        min=2,
        max=10,
        step=1,
        description="Parts",
    )
    
    selection_subdivide_wrapper = HBox([selection_subdivide, selection_num_subdivisions])
    
    selections_elements = VBox(layout=Layout(grid_gap="2px"))
    
    selections_predicates_css = """
    <style>
    .jupyter-scatter-dimbridge-selections-predicates {
        position: absolute !important;
    }
    
    .jupyter-scatter-dimbridge-selections-predicates-wrapper {
        position: relative;
    }
    </style>
    """
    
    display(HTML(selections_predicates_css))
    
    selections_predicates = VBox(
        layout=Layout(
            top="4px",
            left="0px",
            right="0px",
            bottom="4px",
            grid_gap="4px",
        )
    )
    selections_predicates.add_class("jupyter-scatter-dimbridge-selections-predicates")
    
    selections_predicates_wrapper = VBox(
        [selections_predicates],
        layout=Layout(
            height="100%",
        ),
    )
    selections_predicates_wrapper.add_class(
        "jupyter-scatter-dimbridge-selections-predicates-wrapper"
    )
    
    compute_predicates = Button(
        description="Compute Directional Search",
        style="primary",
        disabled=True,
        full_width=True,
    )
    
    compute_predicates_between_selections = Checkbox(
        value=False, description="Compare Between Selections", indent=False
    )
    
    compute_predicates_wrapper = VBox([compute_predicates])
    
    
    def add_selection_element(selection: Selection):
        hex_color = to_hex(selection.color)
    
        selection_name = Label(
            name=selection.name,
            style={"background": hex_color},
        )
    
        selection_remove = Button(
            description="",
            tooltip="Remove Selection",
            icon="trash",
            width=36,
            background=hex_color,
            rounded=["top-right", "bottom-right"],
        )
    
        element = GridBox(
            [
                selection_name,
                selection_remove,
            ],
            layout=Layout(grid_template_columns="1fr 40px"),
        )
    
        def focus_handler(change):
            if change["new"]:
                scatter.zoom(to=selection.points, animation=500, padding=2)
            else:
                scatter.zoom(to=None, animation=500, padding=0)
    
        selection_name.observe(focus_handler, names=["focus"])
    
        def remove_handler(change):
            selections_elements.children = [
                e for e in selections_elements.children if e != element
            ]
            selections.selections = [s for s in selections.selections if s != selection]
            update_annotations()
            compute_predicates.disabled = len(selections.selections) == 0
    
        selection_remove.on_click(remove_handler)
    
        selections_elements.children = selections_elements.children + (element,)
    
    
    def add_subdivided_selections():
        lasso_polygon = scatter.widget.lasso_selection_polygon
        lasso_points = lasso_polygon.shape[0]
    
        lasso_mid = int(lasso_polygon.shape[0] / 2)
        lasso_spine = (lasso_polygon[:lasso_mid, :] + lasso_polygon[lasso_mid:, :]) / 2
    
        lasso_part_one = lasso_polygon[:lasso_mid, :]
        lasso_part_two = lasso_polygon[lasso_mid:, :][::-1]
    
        n_split_points = selection_num_subdivisions.value + 1
    
        sub_lassos_part_one = split_line_equidistant(lasso_part_one, n_split_points)
        sub_lassos_part_two = split_line_equidistant(lasso_part_two, n_split_points)
    
        base_name = selection_name.value
        if len(base_name) == 0:
            base_name = f"Selection {len(selections.selections) + 1}"
    
        color_map = continuous_color_maps[selection_num_subdivisions.value]
    
        for i, part_one in enumerate(sub_lassos_part_one):
            polygon = np.vstack((part_one, sub_lassos_part_two[i][::-1]))
            idxs = np.where(points_in_polygon(df[["x", "y"]].values, polygon))[0]
            points = df.iloc[idxs][["x", "y"]].values
            hull = ConvexHull(points)
            hull_points = np.vstack((points[hull.vertices], points[hull.vertices[0]]))
            color = color_map[i]
            name = f"{base_name}.{i + 1}"
    
            lasso_polygon = polygon.tolist()
            lasso_polygon.append(lasso_polygon[0])
    
            selection = Selection(
                index=len(selections.selections) + 1,
                name=name,
                points=idxs,
                color=color,
                lasso=Line(lasso_polygon),
                hull=Line(hull_points, line_color=color, line_width=2),
            )
            selections.selections.append(selection)
            add_selection_element(selection)
    
    
    def add_selection():
        idxs = scatter.selection()
        points = df.iloc[idxs][["x", "y"]].values
        hull = ConvexHull(points)
        hull_points = np.vstack((points[hull.vertices], points[hull.vertices[0]]))
        color = available_colors.pop(0)
    
        name = selection_name.value
        if len(name) == 0:
            name = f"Selection {len(selections.selections) + 1}"
    
        lasso_polygon = scatter.widget.lasso_selection_polygon.tolist()
        lasso_polygon.append(lasso_polygon[0])
    
        selection = Selection(
            index=len(selections.selections) + 1,
            name=name,
            points=idxs,
            color=color,
            lasso=Line(lasso_polygon),
            hull=Line(hull_points, line_color=color, line_width=2),
        )
        selections.selections.append(selection)
        add_selection_element(selection)
    
    
    def selection_add_handler(event):
        lasso.polygon = None
    
        if scatter.widget.lasso_type == "brush" and selection_subdivide.value:
            add_subdivided_selections()
        else:
            add_selection()
    
        compute_predicates.disabled = False
    
        scatter.selection([])
        update_annotations()
    
        if len(selections.selections) > 1:
            compute_predicates_wrapper.children = (
                compute_predicates_between_selections,
                compute_predicates,
            )
        else:
            compute_predicates_wrapper.children = (compute_predicates,)
    
    
    selection_add.on_click(selection_add_handler)
    
    
    def selection_handler(change):
        if len(change["new"]) > 0:
            selection_add.disabled = False
            selection_name.disabled = False
            selection_name.placeholder = "Name selection…"
            new_index = 1
            if len(selections.selections) > 0:
                new_index = selections.selections[-1].index + 1
            selection_name.value = f"Selection {new_index}"
        else:
            selection_add.disabled = True
            selection_name.disabled = True
            selection_name.placeholder = "Select some points…"
            selection_name.value = ""
    
    
    scatter.widget.observe(selection_handler, names=["selection"])
    
    
    def clear_predicates(event):
        compute_predicates.style = "primary"
        compute_predicates.description = "Compute Predicates"
        compute_predicates.on_click(compute_predicates_handler)
    
        selections_predicates.children = ()
    
        if len(selections.selections) > 1:
            compute_predicates_wrapper.children = (
                compute_predicates_between_selections,
                compute_predicates,
            )
        else:
            compute_predicates_wrapper.children = (compute_predicates,)
    
    
    import ipywidgets as widgets
    from IPython.display import display
    
    
    def fetch_pathways(gene):
        """Fetch Reactome pathways for a given gene symbol."""
        url = f"https://reactome.org/ContentService/data/mapping/UniProt/{gene}/pathways?species=9606"
        try:
            response = requests.get(url)
            response.raise_for_status()
            pathways = response.json()
            return [
                {"Pathway": entry["displayName"], "stId": entry["stId"]}
                for entry in pathways
            ]
        except requests.exceptions.RequestException as e:
            print(f"Error fetching Reactome pathways for {gene}: {e}")
            return []
    
    
    def fetch_pathway_image(pathway_id):
        """Fetch Reactome pathway diagram image."""
        url = f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.png"
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.content
        except requests.exceptions.RequestException as e:
            print(f"Error fetching pathway image: {e}")
            return None
    
    
    search_gene = widgets.Text()
    
    
    def show_directional_results(directional_results):
        # Display the results of the directional analysis as a table.
        # Args:directional_results (list): List of computed correlations from directional analysis.
    
        compute_predicates.style = ""
        compute_predicates.description = "Clear Results"
        compute_predicates.on_click(clear_predicates)  # Attach clear button
    
        all_results = []
    
        for i, result in enumerate(directional_results):
            for entry in result:
                all_results.append(
                    {
                        "Gene": entry["attribute"],
                        "R": round(entry["interval"][0], 4),
                        "p": round(entry["interval"][1], 6),
                    }
                )
    
        # Convert to DataFrame and sort by absolute R-value
        results_df = pd.DataFrame(all_results)
        # Filter out rows where 'R' or 'p' are NaN
        results_df = results_df.dropna(subset=["R", "p"])
        # Sort after removing NaNs
        results_df = results_df.sort_values(by="R", ascending=False).reset_index(drop=True)
    
        # create interactive table with click support
        # existing gene correlation table widget (already displayed):
        gene_table_widget = CorrelationTable(data=results_df.to_dict(orient="records"))
    
        # create new widgets explicitly for Reactome pathway table and diagram
        pathway_table_widget = PathwayTable(data=[])  # initially empty
    
        pathway_table_container.layout.display = "none"
        reactome_diagram_container.layout.display = "none"
    
        # link the selected gene in table to GenePathwayWidget
        # handlers for interactive selections
        def on_gene_click(change):
            gene = change["new"]
            pathways = fetch_pathways(gene)
            pathway_table_widget.data = pathways
            pathway_table_container.layout.display = "block" if pathways else "none"
            reactome_diagram_container.layout.display = "none"
    
        import base64
        import requests
        from ipywidgets import HTML
    
        # instantiate the widget only once, outside the click handler
        interactive_svg_widget = InteractiveSVG()
    
        def on_pathway_click(change):
            pathway_id = change["new"]
            svg_url = (
                f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.svg"
            )
    
            try:
                response = requests.get(svg_url)
                response.raise_for_status()
                svg_content = response.text
                svg_base64 = base64.b64encode(svg_content.encode("utf-8")).decode("utf-8")
    
                interactive_svg_widget.svg_content = svg_base64
                reactome_diagram_container.layout.display = "block"
                reactome_diagram_container.children = [interactive_svg_widget]
    
            except requests.exceptions.RequestException as e:
                print(f"Error fetching SVG diagram: {e}")
    
        # connect handlers
        gene_table_widget.observe(on_gene_click, names=["selected_gene"])
        pathway_table_widget.observe(
            on_pathway_click, names=["selected_pathway"]
        )  # use selected_pathway traitlet
    
        # Show in the UI
        selections_predicates.children = [gene_table_widget]
    
        pathway_table_container.children = [
            widgets.HTML("<b>Reactome Pathways</b>"),
            pathway_table_widget,
        ]
    
        # reactome_diagram_container.children = [pathway_image_widget]
    
        print("Showing directional results...")
    
    
    #############################Part 1
    
    import numpy as np
    import scipy.stats as ss
    import pandas as pd
    
    
    def compute_directional_analysis(df, selections):
        # Computes the correlation of gene expression along a directional axis.
        # Args:
        # df (pd.DataFrame): Dataframe containing gene expression data and spatial coordinates.
        # selections (Selections): The selected points for directional analysis.
        # Returns:
        #     list: A list of dictionaries containing the computed correlations.
    
        if len(selections.selections) == 0:
            return []
    
        results = []
    
        for selection in selections.selections:
            selected_indices = selection.points
            selected_embeddings = df.iloc[selected_indices][["x", "y"]].values
    
            # Ensure we have at least two points for a valid direction vector
            if selected_embeddings.shape[0] < 2:
                continue
    
            # Compute direction vector
            v = selected_embeddings[-1] - selected_embeddings[0]
            v = v / np.linalg.norm(v)  # Normalize
    
            # Compute projections
            start_point = selected_embeddings[0]
            projections = np.array(
                [np.dot(pt - start_point, v) for pt in selected_embeddings]
            )
    
            # Get gene expression data
            columns_to_drop = [
                col
                for col in [
                    "x",
                    "y",
                    "dpi",
                    "strain",
                    "percent.cmv.log10",
                    "seurat_clusters",
                    "virus.presence",
                    "Infection_localminima",
                    "UL123_define_infection",
                    "Infection_state",
                    "Infection_state_bkgd",
                ]
                if col in df.columns
            ]
            selected_expression = df.iloc[selected_indices].drop(columns=columns_to_drop)
    
            # Compute correlations
            correlations = []
            for gene in selected_expression.columns:
                r, p = ss.pearsonr(projections, selected_expression[gene])
                correlations.append(
                    {
                        "attribute": gene,
                        "interval": (r, p),
                        "quality": abs(
                            r
                        ),  # Use absolute correlation as a measure of quality
                    }
                )
    
            results.append(correlations)
    
        return results
    
    
    ######################Part 2
    
    
    def compute_predicates_handler(event):
        if len(selections.selections) == 0:
            return
    
        compute_predicates.disabled = True
        compute_predicates.description = "Computing Directional Analysis…"
    
        # Compute directional correlations
        directional_results = compute_directional_analysis(df, selections)
    
        # Display in a table instead of histogram
        show_directional_results(directional_results)
    
        compute_predicates.disabled = False
        from IPython.display import display
        import ipywidgets as widgets
    
        debug_output = widgets.Output()
        display(debug_output)
        with debug_output:
            print("Running directional analysis...")
    
    
    compute_predicates.on_click(compute_predicates_handler)
    
    
    compute_predicates.on_click(compute_predicates_handler)
    
    add = GridBox(
        [
            selection_name,
            selection_add,
        ],
        layout=Layout(grid_template_columns="1fr 40px"),
    )
    
    complete_add = VBox([add], layout=Layout(grid_gap="4px"))
    
    
    def lasso_type_change_handler(change):
        if change["new"] == "brush":
            complete_add.children = (add, selection_subdivide_wrapper)
        else:
            complete_add.children = (add,)
    
    
    scatter.widget.observe(lasso_type_change_handler, names=["lasso_type"])
    
    metadata_cols_lower = [col.lower() for col in metadata_cols]
    
    color_by = Dropdown(
        options=[("Seurat Clusters", "seurat_clusters")]
        + [
            (col.capitalize(), col)
            for col in metadata_cols_lower
            if col not in ["x", "y", "seurat_clusters"]
        ]
        + [(gene, gene) for gene in adata.var_names],
        # [(col.replace('.', ' ').capitalize(), col) for col in categorical_cols if col != 'seurat_clusters'] +
        # [(p, p) for p in protein_cols],
        value=color_by,
        description="Color By:",
    )
    
    
    def color_by_change_handler(change):
        selected_col = change["new"]
        if selected_col in categorical_color_maps:
            scatter.color(by=selected_col, map=categorical_color_maps[selected_col])
        else:
            scatter.color(by=selected_col, map="plasma")
    
        # cmap = color_map if change['new'] == 'seurat_clusters' else 'plasma'
        # scatter.color(by=change['new'], map=cmap)
    
    
    color_by.observe(color_by_change_handler, names=["value"])
    
    # Main scatterplot and color selection
    plot_wrapper = VBox([scatter.show(), color_by])
    
    pathway_table_container = VBox(
        [],
        layout=Layout(
            overflow_y="auto",
            height="400px",
            border="1px solid #ddd",
            padding="10px",
            display="none",
        ),
    )
    
    reactome_diagram_container = VBox(
        [], layout=Layout(overflow_y="auto", height="400px", padding="10px", display="none")
    )
    
    # Sidebar with selection controls
    sidebar = GridBox(
        [
            complete_add,
            selections_elements,
            selections_predicates_wrapper,
            compute_predicates_wrapper,
        ],
        layout=Layout(
            # grid_template_rows='min-content max-content 1fr min-content',
            grid_template_rows="min-content max-content 1fr min-content",
            overflow_y="auto",
            height="800px",
            grid_gap="4px",
            # height='100%',
        ),
    )
    
    # Pathway table (right panel)
    pathway_table_container.layout = Layout(
        overflow_y="auto",
        height="800px",
        border="1px solid #ddd",
        padding="10px",
        display="none",  # initially hidden until gene selection
    )
    
    # Pathway diagram (bottom panel)
    reactome_diagram_container.layout = Layout(
        overflow_y="auto",
        height="800px",
        border="1px solid #ddd",
        padding="10px",
        display="none",  # initially hidden until pathway selection
    )
    
    # Combine top three panels
    top_layout = GridBox(
        [
            plot_wrapper,
            sidebar,
            pathway_table_container,
        ],
        layout=Layout(
            grid_template_columns="2fr 1fr 1fr",
            grid_gap="10px",
            height="auto",
        ),
    )
    
    from IPython.display import display, HTML
    
    display(
        HTML(
            """
    <style>
    .jp-OutputArea-output, .jp-Cell-outputArea, .jp-Notebook {
        overflow: auto !important;
        max-height: none !important;
    }
    </style>
    """
        )
    )
    
    # Final combined layout
    combined_gene_pathway_panel = GridBox(
        [
            VBox([sidebar], layout=Layout(overflow_y="auto", height="800px")),
            VBox(
                [pathway_table_container], layout=Layout(overflow_y="auto", height="800px")
            ),
        ],
        layout=Layout(
            grid_template_columns="3fr 2fr",  # Gene table 60% and pathway table 40%
            grid_gap="5px",
            # overflow='hidden',
            height="800px",
        ),
    )
    
    # Update the top-level GridBox to include only two columns now
    top_layout_updated = GridBox(
        [
            plot_wrapper,
            combined_gene_pathway_panel,  # combined gene/pathway panel
        ],
        layout=Layout(
            grid_template_columns="3fr 2fr",  # Scatterplot 60%, combined panel 40%
            grid_gap="10px",
            # overflow='hidden',
            height="800px",
        ),
    )
    
    # Final updated layout with pathway diagram at the bottom
    final_layout_updated = VBox(
        [
            top_layout_updated,
            reactome_diagram_container,
        ],
        layout=Layout(grid_gap="10px", width="100%", height="auto"),
    )
    
    # Display the final layout
    display(final_layout_updated)
    
    
