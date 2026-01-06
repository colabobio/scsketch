from anywidget import AnyWidget
from traitlets import Unicode


class InteractiveSVG(AnyWidget):
    _esm = """
    export function render({ model, el }) {

        // --- container styling ---
        el.style.position = 'relative';
        el.style.overflow = 'hidden';
        el.style.border = '1px solid #ddd';
        el.style.width = '100%';
        el.style.height = 'auto';
        el.style.minHeight = '0px';

        const container = document.createElement('div');
        container.style.width = '100%';
        container.style.maxHeight = '800px';
        container.style.height = 'auto';
        container.style.overflow = 'auto';
        container.style.cursor = 'grab';
        container.style.position = 'relative';

        const img = document.createElement('img');
        img.style.transformOrigin = 'top left';
        img.style.userSelect = 'none';

        // IMPORTANT: give the SVG real layout size
        img.style.width = '100%';
        img.style.height = 'auto';
        img.style.display = 'block';
        img.style.minHeight = '300px';

        container.appendChild(img);
        el.appendChild(container);

        // --- Zoom buttons ---
        const zoomIn = document.createElement('button');
        zoomIn.textContent = '+';
        zoomIn.style.position = 'absolute';
        zoomIn.style.top = '10px';
        zoomIn.style.right = '50px';
        zoomIn.style.zIndex = '10';

        const zoomOut = document.createElement('button');
        zoomOut.textContent = '−';
        zoomOut.style.position = 'absolute';
        zoomOut.style.top = '10px';
        zoomOut.style.right = '10px';
        zoomOut.style.zIndex = '10';

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
            container.style.cursor = 'grabbing';
        });

        container.addEventListener('mousemove', (e) => {
            if (!dragging) return;
            container.scrollLeft = scrollLeft - (e.pageX - container.offsetLeft - startX);
            container.scrollTop = scrollTop - (e.pageY - container.offsetTop - startY);
        });

        window.addEventListener('mouseup', () => {
            dragging = false;
            container.style.cursor = 'grab';
        });

        // --- Wheel zoom ---
        container.addEventListener('wheel', (e) => {
            e.preventDefault();
            if (e.deltaY < 0) scale = Math.min(scale + step, 10);
            else scale = Math.max(scale - step, 0.1);
            applyZoom();
        });

        // --- Blob URL rendering (Firefox–safe) ---
        let currentURL = null;

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
    """

    svg_content = Unicode("").tag(sync=True)
