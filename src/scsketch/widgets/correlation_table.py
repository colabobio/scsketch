import requests
import traitlets
from anywidget import AnyWidget
from traitlets import Dict, List, Unicode


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
