import requests
import traitlets
from anywidget import AnyWidget
from traitlets import List, Dict, Unicode, Int, Bool

#Additional UI Widgets - To visualize the results of scSketch, we need a few additional widgets:
#Directional Search Interactive Table Widget: a widget to visualiaze results of directional analysis of embedding and see what pathways the most upregulated and downregulated genes are a part of in Reactome.
#Label: a widget to display a selection of points
#Divider: a widget to visually add some clarity between groups of histograms

class GeneProjectionPlot(AnyWidget):
    """Interactive scatterplot of gene expression vs projection."""
    _esm = """
    function render({ model, el }) {
      el.style.display = "flex";
      el.style.flexDirection = "column";
      el.style.alignItems = "center";

        

      const canvas = document.createElement("canvas");
      el.appendChild(canvas);
    
      // Create the drawing context FIRST
      const ctx = canvas.getContext("2d");
       
      // --- DPI-aware responsive sizing (fixed for Retina + Firefox) ---
      let currentDPR = window.devicePixelRatio || 1;
    
      function sizeCanvas() {
        const r = el.getBoundingClientRect();
        const dpr = window.devicePixelRatio || 1;
        
        // Detect Retina / XDR environment
        const isHighDPI = dpr >= 1.7;   // MacBook Pro Retina / XDR screens
        const screenW = window.innerWidth;
        const isSmallScreen = screenW < 1600;
        
        // --- Responsive aspect tuning ---
        // XDR: make flatter; otherwise use standard
        const baseAspect = isHighDPI ? (isSmallScreen ? 0.32 : 0.38)
                                       : (isSmallScreen ? 0.40 : 0.50);
        // --- Determine available height from parent ---
        let parentMax = 0;
        try {
          // read computed CSS height if parent sets one
          const style = window.getComputedStyle(el);
          const maxH = parseFloat(style.maxHeight || "0");
          const explicitH = parseFloat(style.height || "0");
          parentMax = Math.max(maxH, explicitH);
        } catch {}
        
        // --- Calculate adaptive height (smaller content, respects parent height) ---
        let parentH = el.parentElement?.getBoundingClientRect()?.height || 180;
        let targetHeight = Math.min(parentH * 0.95, r.width * baseAspect);
        
        // Cap final height; keep modest on Retina
        targetHeight = Math.max(180, targetHeight);
        
        // Apply CSS sizing
        canvas.style.width = "100%";
        canvas.style.height = `${targetHeight}px`;
        el.style.height = `${targetHeight + 10}px`;
        el.style.overflow = "hidden";
        // --- Apply CSS size ---
        canvas.style.width = "100%";
        canvas.style.height = `${targetHeight}px`;
        
        // --- Clip outer container so nothing overflows ---
        el.style.maxHeight = `${targetHeight + 20}px`;
        el.style.overflow = "hidden";
        
        // --- Apply CSS size ---
        canvas.style.width = "100%";
        canvas.style.height = `${targetHeight}px`;
        
        // --- Constrain outer container (scroll if overflow) ---
        el.style.maxHeight = `${targetHeight + 50}px`;
        el.style.overflow = "auto";
        
        const cssW = canvas.clientWidth;
        const cssH = canvas.clientHeight;
        
        // --- Retina-aware internal buffer ---
        const pxWidth = Math.round(cssW * dpr);
        const pxHeight = Math.round(cssH * dpr);
        canvas.width = pxWidth;
        canvas.height = pxHeight;
        
        // --- Proper scaling transform (prevent right/bottom cutoff) ---
        ctx.setTransform(1, 0, 0, 1, 0, 0);   // reset transform first
        ctx.scale(dpr, dpr);                   // scale drawing coordinates for Retina
        ctx.clearRect(0, 0, cssW, cssH);       // clear only the visible region
        
        // --- Global font/marker scale for small Retina displays ---
        const scale = isHighDPI ? 0.78 : 1.0;
        window._fontScale = scale;
      }
        
      // --- Observe parent resize ---
      const ro = new ResizeObserver(() => {
        sizeCanvas();
        drawPlot();
      });
      ro.observe(el);
    
      // --- Watch for DPI changes ---
      function watchDPR() {
        let last = currentDPR;
        const loop = () => {
          const newDPR = window.devicePixelRatio || 1;
          if (Math.abs(newDPR - last) > 0.05) {
            last = newDPR;
            console.log(`[DPI] detected ${newDPR}`);
            setTimeout(() => {
              sizeCanvas();
              drawPlot();
            }, 180);
          }
          requestAnimationFrame(loop);
        };
        requestAnimationFrame(loop);
      }
      watchDPR();
    
      // --- Initial sizing ---
      sizeCanvas();
      
      function drawPlot() {

        const cssW = canvas.clientWidth;
        const cssH = canvas.clientHeight;

        ctx.clearRect(0, 0, cssW, cssH);
        
        // Inner device pixel buffer
        // ctx.strokeStyle = "rgba(180,0,180,0.6)";
        // ctx.lineWidth = 1;
        // ctx.strokeRect(0, 0, cssW, cssH);
        
        const data = model.get("data") || [];
        const gene = model.get("gene");
        console.log("drawPlot called", gene, data.length);
    
        const points = data.map(d => ({ x: d.projection, y: d.expression }));

        const fontScale = window._fontScale || 1.0;
        
        // === Title with gene name ===
        ctx.fillStyle = "black";
        ctx.font = `${(cssW < 700 ? 13 : 16) * fontScale}px sans-serif`;
        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        const titleText = gene ? `Gene Expression vs Projection for ${gene}` : "Gene Expression vs Projection";
        
        // If title too long, split it
        const maxWidth = cssW - 40;
        if (ctx.measureText(titleText).width > maxWidth) {
          const parts = titleText.split(" for ");
          ctx.fillText(parts[0], cssW / 2, 8);
          if (parts[1]) ctx.fillText("for " + parts[1], cssW / 2, 26);
        } else {
          ctx.fillText(titleText, cssW / 2, 10);
        }
    
        if (points.length === 0) {
          ctx.fillText("No data", 10, 40);
          return;
        }
    
        const xs = points.map(p => p.x);
        const ys = points.map(p => p.y);
        let minX = Math.min(...xs), maxX = Math.max(...xs);
        let minY = Math.min(...ys), maxY = Math.max(...ys);
    
        if (!Number.isFinite(minX)) minX = 0;
        if (!Number.isFinite(maxX)) maxX = 1;
        if (!Number.isFinite(minY)) minY = 0;
        if (!Number.isFinite(maxY)) maxY = 1;
        if (maxX === minX) { maxX += 1e-6; minX -= 1e-6; }
        if (maxY === minY) { maxY += 1e-6; minY -= 1e-6; }
    
        // === Layout constants (shrink inner chart to fit compact box) ===
        const shrinkFactor = 0.98;  // controls how much to shrink everything inside
        
        const marginLeft = 45;
        const marginRight = 25;
        const marginBottom = 40;
        const marginTop = 25;
        
        const innerWidth = (cssW - marginLeft - marginRight) * shrinkFactor;
        const innerHeight = (cssH - marginTop - marginBottom) * shrinkFactor;
                
        // --- Scaling functions (center shrunken plot within canvas) ---
        const xOffset = 0;
        const yOffset = 0;
        
        const xScale = (x) =>
          marginLeft + ((x - minX) / (maxX - minX)) * innerWidth;
        const yScale = (y) =>
          cssH - marginBottom - ((y - minY) / (maxY - minY)) * innerHeight;
        // --- Light gray grid & axes ---
        ctx.strokeStyle = "#cccccc";
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(marginLeft, cssH - marginBottom);
        ctx.lineTo(cssW - marginRight, cssH - marginBottom); // x-axis
        ctx.moveTo(marginLeft, marginTop);
        ctx.lineTo(marginLeft, cssH - marginBottom); // y-axis
        ctx.stroke();
    
        // --- Tick marks and numeric labels (0, 0.5, 1.0) ---
        const tickVals = [0.0, 0.5, 1.0];
        ctx.font = `${11 * fontScale}px sans-serif`;
        ctx.fillStyle = "#444";
        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        tickVals.forEach(t => {
          const x = xScale(t);
          ctx.beginPath();
          ctx.moveTo(x, cssH - marginBottom);
          ctx.lineTo(x, cssH - marginBottom + 4);
          ctx.stroke();
          ctx.fillText(t.toFixed(1), x, cssH - marginBottom + 6);
        });
    
        ctx.textAlign = "right";
        ctx.textBaseline = "middle";
        tickVals.forEach(t => {
          const y = yScale(t);
          ctx.beginPath();
          ctx.moveTo(marginLeft - 4, y);
          ctx.lineTo(marginLeft, y);
          ctx.stroke();
          ctx.fillText(t.toFixed(1), marginLeft - 8, y);
        });
    
        // --- Points ---
        ctx.fillStyle = "#777";
        const MAX_POINTS = 8000;
        let renderPts = points;
        if (points.length > MAX_POINTS) {
          const step = Math.ceil(points.length / MAX_POINTS);
          renderPts = points.filter((_, i) => i % step === 0);
        }
        renderPts.forEach(p => {
          const x = xScale(p.x);
          const y = yScale(p.y);
          ctx.beginPath();
          ctx.arc(x, y, 2 * fontScale, 0, 2 * Math.PI);
          ctx.fill();
        });
        
        // --- Axis labels (closer to axes) ---
        ctx.fillStyle = "black";
        ctx.font = `${12 * fontScale}px sans-serif`;
        ctx.textAlign = "center";
        ctx.fillText("Projection", marginLeft + innerWidth / 2, cssH - marginBottom + 35);
        
        ctx.save();
        ctx.translate(marginLeft - 40, marginTop + innerHeight / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.textAlign = "center";
        ctx.fillText("Expression", 0, 0);
        ctx.restore();
        // Outer CSS pixel bounds
        // ctx.strokeStyle = "rgba(255,0,0,0.4)";
        // ctx.lineWidth = 2;
        // ctx.strokeRect(0, 0, cssW, cssH);
      }
      model.on("change:data", drawPlot);
      model.on("change:gene", drawPlot);
      drawPlot();
    }
    export default { render };
    """

    data = List(Dict()).tag(sync=True)  # [{"projection": float, "expression": float}]
    gene = Unicode("").tag(sync=True)

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
        const data = model.get("data") || [];

        //Always show these columns, in this order:
        const columns = ["Gene", "R", "p", "Selection"];  //removed "alpha_i" and "reject"
        
        //Header
        const headerRow = document.createElement("tr");
        columns.forEach(col => {
          const th = document.createElement("th");
          th.textContent = col;
          headerRow.appendChild(th);
        });
        table.appendChild(headerRow);

        rowsCache = data.map(row => {
          const tr = document.createElement("tr");
          const geneVal = (row["Gene"] ?? "").toString();
          tr.dataset.gene = geneVal.toLowerCase();
          tr.style.cursor = "pointer";
          tr.onclick = () => {
            if (geneVal){      
              model.set("selected_gene", geneVal);
              model.save_changes();
            }
          };

          columns.forEach(col => {
            const td = document.createElement("td");
            const val = row[col];
            
            if (col === "R" || col === "alpha_i") {
            // format to 4 decimal places if numeric
              const num = Number(val);
              td.textContent = Number.isFinite(num) ? num.toFixed(4) : (val ?? "");
            } else if (col === "p"){
                const num = Number(val);
                td.textContent = Number.isFinite(num) ? num.toExponential(3) : (val ?? "");
            } else if (col === "reject") {
                td.textContent = typeof val === "boolean" ? (val ? "Pass" : "") : (val ?? "");
            } else {
                td.textContent = (val ?? "").toString();
            }
            
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

        // --- container styling ---
        el.style.position = 'relative';
        el.style.overflow = 'hidden';
        el.style.border = '1px solid #ddd';
        el.style.width = '100%';
        el.style.height = 'auto';
        el.style.minHeight = '0px';

        const container = document.createElement('div');
        container.style.width = '100%';
        container.style.maxHeight = '800px';
        container.style.height = 'auto';
        container.style.overflow = 'auto';
        container.style.cursor = 'grab';
        container.style.position = 'relative';

        const img = document.createElement('img');
        img.style.transformOrigin = 'top left';
        img.style.userSelect = 'none';

        // IMPORTANT: give the SVG real layout size
        img.style.width = '100%';
        img.style.height = 'auto';
        img.style.display = 'block';
        img.style.minHeight = '300px';

        container.appendChild(img);
        el.appendChild(container);

        // --- Zoom buttons ---
        const zoomIn = document.createElement('button');
        zoomIn.textContent = '+';
        zoomIn.style.position = 'absolute';
        zoomIn.style.top = '10px';
        zoomIn.style.right = '50px';
        zoomIn.style.zIndex = '10';

        const zoomOut = document.createElement('button');
        zoomOut.textContent = '−';
        zoomOut.style.position = 'absolute';
        zoomOut.style.top = '10px';
        zoomOut.style.right = '10px';
        zoomOut.style.zIndex = '10';

        el.appendChild(zoomIn);
        el.appendChild(zoomOut);

        // --- Zoom state ---
        let scale = 1;
        const step = 0.1;

        function applyZoom() {
            img.style.transform = `scale(${scale})`;
        }

        zoomIn.onclick = () => { scale = Math.min(scale + step, 10); applyZoom(); };
        zoomOut.onclick = () => { scale = Math.max(scale - step, 0.1); applyZoom(); };

        // --- Drag-to-pan ---
        let dragging = false, startX = 0, startY = 0, scrollLeft = 0, scrollTop = 0;

        container.addEventListener('mousedown', (e) => {
            dragging = true;
            startX = e.pageX - container.offsetLeft;
            startY = e.pageY - container.offsetTop;
            scrollLeft = container.scrollLeft;
            scrollTop = container.scrollTop;
            container.style.cursor = 'grabbing';
        });

        container.addEventListener('mousemove', (e) => {
            if (!dragging) return;
            container.scrollLeft = scrollLeft - (e.pageX - container.offsetLeft - startX);
            container.scrollTop = scrollTop - (e.pageY - container.offsetTop - startY);
        });

        window.addEventListener('mouseup', () => {
            dragging = false;
            container.style.cursor = 'grab';
        });

        // --- Wheel zoom ---
        container.addEventListener('wheel', (e) => {
            e.preventDefault();
            if (e.deltaY < 0) scale = Math.min(scale + step, 10);
            else scale = Math.max(scale - step, 0.1);
            applyZoom();
        });

        // --- Blob URL rendering (Firefox–safe) ---
        let currentURL = null;

        const update = () => {
            const svgContent = model.get("svg_content");
            if (!svgContent) return;

            img.src = `data:image/svg+xml;base64, ${svgContent}`;
            scale = 1;
            applyZoom();
        };

        model.on("change:svg_content", update);
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