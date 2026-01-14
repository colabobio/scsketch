import traitlets
from anywidget import AnyWidget
from traitlets import Dict, List, Unicode


class GeneViolinPlot(AnyWidget):
    """Compact violin-like plot for selected vs background distributions.

    This widget expects binned densities, not raw points, to keep payloads small.
    """

    _esm = """
    function render({ model, el }) {
      el.style.display = "flex";
      el.style.flexDirection = "column";
      el.style.alignItems = "stretch";
      el.style.width = "100%";

      const canvas = document.createElement("canvas");
      canvas.style.width = "100%";
      canvas.style.height = "220px";
      el.appendChild(canvas);

      const ctx = canvas.getContext("2d");

      function sizeCanvas() {
        const dpr = window.devicePixelRatio || 1;
        const cssW = el.getBoundingClientRect().width || 600;
        const cssH = 220;
        canvas.width = Math.round(cssW * dpr);
        canvas.height = Math.round(cssH * dpr);
        canvas.style.height = cssH + "px";
        ctx.setTransform(1, 0, 0, 1, 0, 0);
        ctx.scale(dpr, dpr);
      }

      const ro = new ResizeObserver(() => {
        sizeCanvas();
        draw();
      });
      ro.observe(el);

      function drawViolin(xCenter, yVals, densVals, maxDens, color) {
        if (!yVals || !densVals || yVals.length === 0) return;
        ctx.beginPath();
        for (let i = 0; i < yVals.length; i++) {
          const y = yVals[i];
          const w = (densVals[i] / (maxDens || 1)) * 80; // half-width in px
          ctx.lineTo(xCenter + w, y);
        }
        for (let i = yVals.length - 1; i >= 0; i--) {
          const y = yVals[i];
          const w = (densVals[i] / (maxDens || 1)) * 80;
          ctx.lineTo(xCenter - w, y);
        }
        ctx.closePath();
        ctx.fillStyle = color;
        ctx.globalAlpha = 0.25;
        ctx.fill();
        ctx.globalAlpha = 1.0;
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.5;
        ctx.stroke();
      }

      function draw() {
        const cssW = canvas.clientWidth || 600;
        const cssH = 220;
        ctx.clearRect(0, 0, cssW, cssH);

        const gene = model.get("gene") || "";
        const payload = model.get("data") || {};
        const bins = payload.bins || [];
        const sel = payload.selected || [];
        const bg = payload.background || [];

        // Title
        ctx.fillStyle = "#111";
        ctx.font = "14px sans-serif";
        ctx.textAlign = "left";
        ctx.textBaseline = "top";
        const title = gene ? `Differential distribution: ${gene}` : "Differential distribution";
        ctx.fillText(title, 10, 8);

        if (bins.length < 2 || sel.length === 0 || bg.length === 0) {
          ctx.fillStyle = "#555";
          ctx.font = "12px sans-serif";
          ctx.fillText("No data", 10, 32);
          return;
        }

        const top = 34;
        const bottom = cssH - 22;
        const left = 50;
        const right = cssW - 20;

        // y mapping (expression axis)
        const minY = bins[0];
        const maxY = bins[bins.length - 1];
        const yScale = (v) => {
          const t = (v - minY) / ((maxY - minY) || 1);
          return bottom - t * (bottom - top);
        };

        // Precompute y points at bin centers
        const yCenters = [];
        for (let i = 0; i < bins.length - 1; i++) {
          yCenters.push((bins[i] + bins[i + 1]) / 2);
        }
        const yPix = yCenters.map(yScale);

        const maxSel = Math.max(...sel);
        const maxBg = Math.max(...bg);
        const maxDens = Math.max(maxSel, maxBg, 1e-12);

        // Axes
        ctx.strokeStyle = "#ccc";
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(left, top);
        ctx.lineTo(left, bottom);
        ctx.stroke();

        // y ticks
        ctx.fillStyle = "#444";
        ctx.font = "11px sans-serif";
        ctx.textAlign = "right";
        ctx.textBaseline = "middle";
        const ticks = 4;
        for (let i = 0; i <= ticks; i++) {
          const v = minY + (i / ticks) * (maxY - minY);
          const y = yScale(v);
          ctx.strokeStyle = "#eee";
          ctx.beginPath();
          ctx.moveTo(left, y);
          ctx.lineTo(right, y);
          ctx.stroke();
          ctx.fillStyle = "#444";
          ctx.fillText(v.toFixed(2), left - 6, y);
        }

        // Two violins: selected (left) vs background (right)
        const xSel = left + (right - left) * 0.33;
        const xBg = left + (right - left) * 0.67;
        drawViolin(xSel, yPix, sel, maxDens, "#6a00ff");
        drawViolin(xBg, yPix, bg, maxDens, "#999999");

        // Labels
        ctx.fillStyle = "#111";
        ctx.font = "12px sans-serif";
        ctx.textAlign = "center";
        ctx.textBaseline = "bottom";
        ctx.fillText("Selected", xSel, cssH - 4);
        ctx.fillText("Background", xBg, cssH - 4);
      }

      model.on("change:data", () => draw());
      model.on("change:gene", () => draw());
      sizeCanvas();
      draw();
    }
    export default { render };
    """

    data = Dict(default_value={}).tag(sync=True)
    gene = Unicode("").tag(sync=True)

