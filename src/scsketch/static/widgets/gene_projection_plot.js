function render({ model, el }) {
  el.classList.add("gene-projection-plot");

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
    const baseAspect = isHighDPI ? (isSmallScreen ? 0.32 : 0.38)
                                 : (isSmallScreen ? 0.40 : 0.50);

    // --- Calculate adaptive height ---
    let parentH = el.parentElement?.getBoundingClientRect()?.height || 180;
    let targetHeight = Math.min(parentH * 0.95, r.width * baseAspect);
    targetHeight = Math.max(180, targetHeight);

    canvas.style.height = `${targetHeight}px`;
    el.style.height = `${targetHeight + 10}px`;
    el.style.maxHeight = `${targetHeight + 50}px`;

    const cssW = canvas.clientWidth;
    const cssH = canvas.clientHeight;

    // --- Retina-aware internal buffer ---
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);

    // --- Proper scaling transform ---
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, cssW, cssH);

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

    const data = model.get("data") || [];
    const gene = model.get("gene");

    const points = data.map(d => ({ x: d.projection, y: d.expression }));
    const fontScale = window._fontScale || 1.0;

    // === Title with gene name ===
    ctx.fillStyle = "black";
    ctx.font = `${(cssW < 700 ? 13 : 16) * fontScale}px sans-serif`;
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    const titleText = gene
      ? `Gene Expression vs Projection for ${gene}`
      : "Gene Expression vs Projection";

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

    const shrinkFactor = 0.98;
    const marginLeft = 45;
    const marginRight = 25;
    const marginBottom = 40;
    const marginTop = 25;

    const innerWidth = (cssW - marginLeft - marginRight) * shrinkFactor;
    const innerHeight = (cssH - marginTop - marginBottom) * shrinkFactor;

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

    // --- Tick marks and numeric labels ---
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

    // --- Axis labels ---
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
  }

  model.on("change:data", drawPlot);
  model.on("change:gene", drawPlot);
  drawPlot();
}
export default { render };
