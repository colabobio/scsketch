export function render({ model, el }) {
  el.classList.add("isvg-host");

  const container = document.createElement('div');
  container.classList.add("isvg-container");

  const img = document.createElement('img');
  img.classList.add("isvg-image");

  container.appendChild(img);
  el.appendChild(container);

  // --- Zoom buttons ---
  const zoomIn = document.createElement('button');
  zoomIn.textContent = '+';
  zoomIn.classList.add("isvg-btn", "isvg-btn--zoom-in");

  const zoomOut = document.createElement('button');
  zoomOut.textContent = '−';
  zoomOut.classList.add("isvg-btn", "isvg-btn--zoom-out");

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
    container.classList.add("isvg-container--grabbing");
  });

  container.addEventListener('mousemove', (e) => {
    if (!dragging) return;
    container.scrollLeft = scrollLeft - (e.pageX - container.offsetLeft - startX);
    container.scrollTop = scrollTop - (e.pageY - container.offsetTop - startY);
  });

  window.addEventListener('mouseup', () => {
    dragging = false;
    container.classList.remove("isvg-container--grabbing");
  });

  // --- Wheel zoom ---
  container.addEventListener('wheel', (e) => {
    e.preventDefault();
    if (e.deltaY < 0) scale = Math.min(scale + step, 10);
    else scale = Math.max(scale - step, 0.1);
    applyZoom();
  });

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
