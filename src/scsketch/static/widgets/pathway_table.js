function render({ model, el }) {
  const table = document.createElement("table");
  table.classList.add("pathway-table");
  el.appendChild(table);

  const render = () => {
    table.innerHTML = "";
    const data = model.get("data") || [];

    if (data.length === 0) return;

    const headerRow = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = "Pathway";
    headerRow.appendChild(th);
    table.appendChild(headerRow);

    data.forEach(pathway => {
      const row = document.createElement("tr");
      row.style.cursor = "pointer";
      row.onclick = () => {
        model.set("selected_pathway", pathway.stId);
        model.save_changes();
      };

      const td = document.createElement("td");
      td.textContent = pathway.name;
      row.appendChild(td);
      table.appendChild(row);
    });
  };

  model.on("change:data", render);
  render();
}
export default { render };
