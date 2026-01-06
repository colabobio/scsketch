from anywidget import AnyWidget
import traitlets


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
