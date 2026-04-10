function render({ model, el }) {
  el.classList.add("gene-pathway-widget");

  const geneDropdown = document.createElement("select");
  model.get("genes").forEach((gene) => {
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
    model.set("selected_gene", geneDropdown.value);
    model.save_changes();
  });

  pathwayDropdown.addEventListener("change", () => {
    model.set("selected_pathway", pathwayDropdown.value);
    model.save_changes();
  });

  model.on("change:pathways", () => {
    const pathways = model.get("pathways");
    pathwayDropdown.innerHTML = "";
    if (pathways.length > 0) {
      pathwayDropdown.style.display = "block";
      pathways.forEach((pathway) => {
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
