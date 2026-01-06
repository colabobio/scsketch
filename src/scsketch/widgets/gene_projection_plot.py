import requests
import traitlets
from anywidget import AnyWidget
from traitlets import Bool, Dict, Int, List, Unicode


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
