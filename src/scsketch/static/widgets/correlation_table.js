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
    if (col === "Discovery Score") {
      return "ct-col-score";
    }
    if (col === "R" || col === "T" || col === "alpha_i") {
      return "ct-col-stat";
    }
    return "ct-col-metric";
  };

  const initializeTable = () => {
    const data = model.get("data") || [];

    const columns = model.get("columns") || ["Gene", "R", "Discovery Score", "Selection"];

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
      const geneId = (row["_gene_id"] ?? row["Gene"] ?? "").toString();
      const geneVal = (row["Gene"] ?? geneId).toString();
      tr.dataset.gene = `${geneVal} ${geneId}`.toLowerCase();
      tr.style.cursor = "pointer";
      tr.onclick = () => {
        if (geneId) {
          model.set("selected_gene", geneId);
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
        } else if (col === "Discovery Score") {
          const num = Number(val);
          td.textContent = Number.isFinite(num) ? Math.round(num).toString() : (val ?? "");
        } else if (col === "p") {
          const num = Number(val);
          td.textContent = Number.isFinite(num) ? num.toExponential(3) : (val ?? "");
        } else if (col === "reject") {
          td.textContent = typeof val === "boolean" ? (val ? "Pass" : "") : (val ?? "");
        } else {
          td.textContent = (val ?? "").toString();
        }

        td.title =
          col === "Gene" && geneId && geneId !== td.textContent
            ? `${td.textContent} (${geneId})`
            : td.textContent;
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

    requestAnimationFrame(() => {
      rowsCache.forEach(row => {
        if (row.dataset.gene.includes(filterText)) {
          row.style.display = "table-row";
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
