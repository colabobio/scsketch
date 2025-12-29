"""Interactive SVG viewer widget with zoom and pan capabilities."""

from anywidget import AnyWidget
from traitlets import Unicode


class InteractiveSVG(AnyWidget):
    """Zoomable and pannable SVG viewer for pathway diagrams."""

    _esm = """
    export function render({ model, el }) {
        el.style.position = 'relative';
        el.style.overflow = 'hidden';
        el.style.border = '1px solid #ddd';
        el.style.width = '100%';
        el.style.height = '800px';

        const container = document.createElement('div');
        container.style.width = '100%';
        container.style.height = '100%';
        container.style.overflow = 'auto';
        container.style.cursor = 'grab';
        container.style.position = 'relative';

        const img = document.createElement('img');
        img.style.transformOrigin = 'top left';
        img.style.width = '100%';
        img.style.height = 'auto';
        img.style.userSelect = 'none';

        let scale = 1;
        const scaleStep = 0.1;
        const minScale = 0.1;
        const maxScale = 10;

        const zoomInBtn = document.createElement('button');
        zoomInBtn.innerHTML = '+';
        zoomInBtn.style.position = 'absolute';
        zoomInBtn.style.top = '10px';
        zoomInBtn.style.right = '50px';
        zoomInBtn.style.zIndex = '10';

        const zoomOutBtn = document.createElement('button');
        zoomOutBtn.innerHTML = '−';
        zoomOutBtn.style.position = 'absolute';
        zoomOutBtn.style.top = '10px';
        zoomOutBtn.style.right = '10px';
        zoomOutBtn.style.zIndex = '10';

        function applyTransform() {
            img.style.transform = `scale(${scale})`;
        }

        zoomInBtn.onclick = () => {
            scale = Math.min(scale + scaleStep, maxScale);
            applyTransform();
        };

        zoomOutBtn.onclick = () => {
            scale = Math.max(scale - scaleStep, minScale);
            applyTransform();
        };

        container.appendChild(img);
        el.appendChild(container);
        el.appendChild(zoomInBtn);
        el.appendChild(zoomOutBtn);

        // Drag-to-pan
        let isDragging = false;
        let startX, startY, scrollLeft, scrollTop;

        container.addEventListener('mousedown', (e) => {
            isDragging = true;
            startX = e.pageX - container.offsetLeft;
            startY = e.pageY - container.offsetTop;
            scrollLeft = container.scrollLeft;
            scrollTop = container.scrollTop;
            container.style.cursor = 'grabbing';
            e.preventDefault();
        });

        container.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            const x = e.pageX - container.offsetLeft;
            const y = e.pageY - container.offsetTop;
            const walkX = x - startX;
            const walkY = y - startY;
            container.scrollLeft = scrollLeft - walkX;
            container.scrollTop = scrollTop - walkY;
        });

        window.addEventListener('mouseup', () => {
            isDragging = false;
            container.style.cursor = 'grab';
        });

        // Mouse wheel for zooming
        container.addEventListener('wheel', (e) => {
            e.preventDefault();
            const oldScale = scale;
            if (e.deltaY < 0) {
                scale = Math.min(scale + scaleStep, maxScale);
            } else {
                scale = Math.max(scale - scaleStep, minScale);
            }

            // Calculate zoom towards the mouse cursor
            const rect = container.getBoundingClientRect();
            const offsetX = (e.clientX - rect.left) + container.scrollLeft;
            const offsetY = (e.clientY - rect.top) + container.scrollTop;
            const dx = (offsetX / oldScale) * (scale - oldScale);
            const dy = (offsetY / oldScale) * (scale - oldScale);

            applyTransform();
            container.scrollLeft += dx;
            container.scrollTop += dy;
        });

        const update = () => {
            const svgContent = model.get('svg_content');
            img.src = `data:image/svg+xml;base64,${svgContent}`;
            scale = 1;
            applyTransform();
        };

        model.on('change:svg_content', update);
        update();
    }
    """

    svg_content = Unicode("").tag(sync=True)
