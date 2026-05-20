function render({ model, el }) {
  const container = document.createElement("div");
  const searchInput = document.createElement("input");
  searchInput.type = "text";
  searchInput.placeholder = "Search genes...";
  searchInput.classList.add("ct-search-input");

  const tableWrap = document.createElement("div");
  tableWrap.classList.add("ct-table-wrap");

  const table = document.createElement("table");
  table.classList.add("correlation-table");
  const thead = document.createElement("thead");
  const tbody = document.createElement("tbody");
  table.appendChild(thead);
  table.appendChild(tbody);

  container.appendChild(searchInput);
  tableWrap.appendChild(table);
  container.appendChild(tableWrap);
  el.appendChild(container);

  let rowsCache = [];
  const MAX_ROWS = 200; // minimum visible rows at a time

  const columnClass = (col) => {
    if (col === "Gene") {
      return "ct-col-gene";
    }
    if (col === "Selection") {
      return "ct-col-selection";
    }
    if (col === "p") {
      return "ct-col-pvalue";
    }
    return "ct-col-metric";
  };

  const initializeTable = () => {
    const data = model.get("data") || [];

    // Always show these columns, in this order:
    // Default is directional mode: ["Gene", "R", "p", "Selection"].
    const columns = model.get("columns") || ["Gene", "R", "p", "Selection"];

    // Header
    const headerRow = document.createElement("tr");
    columns.forEach(col => {
      const th = document.createElement("th");
      th.textContent = col;
      th.title = col;
      th.classList.add(columnClass(col));
      headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);

    rowsCache = data.map(row => {
      const tr = document.createElement("tr");
      const geneVal = (row["Gene"] ?? "").toString();
      tr.dataset.gene = geneVal.toLowerCase();
      tr.style.cursor = "pointer";
      tr.onclick = () => {
        if (geneVal) {
          model.set("selected_gene", geneVal);
          model.save_changes();
        }
      };

      columns.forEach(col => {
        const td = document.createElement("td");
        const val = row[col];
        td.classList.add(columnClass(col));

        if (col === "R" || col === "T" || col === "alpha_i") {
          // format to 4 decimal places if numeric
          const num = Number(val);
          td.textContent = Number.isFinite(num) ? num.toFixed(4) : (val ?? "");
        } else if (col === "p") {
          const num = Number(val);
          td.textContent = Number.isFinite(num) ? num.toExponential(3) : (val ?? "");
        } else if (col === "reject") {
          td.textContent = typeof val === "boolean" ? (val ? "Pass" : "") : (val ?? "");
        } else {
          td.textContent = (val ?? "").toString();
        }

        td.title = td.textContent;
        tr.appendChild(td);
      });

      tbody.appendChild(tr);
      return tr; // caching the row
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
      timeout = setTimeout(() => func.apply(this, args), wait);
    };
  }

  searchInput.addEventListener("input", debounce(() => {
    const currentLength = searchInput.value.length;
    debounce(updateTable, currentLength < previousLength ? 300 : 200)();
    previousLength = currentLength;
  }, 50));
}
export default { render };
