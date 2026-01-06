from anywidget import AnyWidget
from traitlets import Bool, Dict, Unicode


class Label(AnyWidget):
    _esm = """
    function render({ model, el }) {
      const label = document.createElement("div");
      label.classList.add(
        'jupyter-widgets',
        'jupyter-scatter-label'
      );
      label.tabIndex = 0;
      
      const update = () => {
        label.textContent = model.get('name');

        for (const [key, value] of Object.entries(model.get('style'))) {
          label.style[key] = value;
        }
      }
      
      model.on('change:name', update);
      model.on('change:style', update);

      update();

      const createFocusChanger = (value) => () => {
        model.set('focus', value);
        model.save_changes();
      }

      const focusHandler = createFocusChanger(true);
      const blurHandler = createFocusChanger(false);

      label.addEventListener('focus', focusHandler);
      label.addEventListener('blur', blurHandler);

      el.appendChild(label);

      const updateFocus = () => {
        if (model.get('focus')) {
          label.focus();
        }
      }
      
      model.on('change:focus', updateFocus);

      window.requestAnimationFrame(() => {
        updateFocus();
      });

      return () => {
        label.removeEventListener('focus', focusHandler);
        label.removeEventListener('blur', blurHandler);
      }
    }
    export default { render };
    """

    _css = """
    .jupyter-scatter-label {
      display: flex;
      align-items: center;
      width: 100%;
      height: var(--jp-widgets-inline-height);
      padding: var(--jp-widgets-input-padding) calc(var(--jp-widgets-input-padding)* 2);
      border-top-left-radius: var(--jp-border-radius);
      border-rop-right-radius: 0;
      border-bottom-left-radius: var(--jp-border-radius);
      border-bottom-right-radius: 0;
    }
    .jupyter-scatter-label:focus {
      font-weight: bold;
      outline: 1px solid var(--jp-widgets-input-focus-border-color);
      outline-offset: 1px;
    }
    """

    name = Unicode("").tag(sync=True)
    style = Dict({}).tag(sync=True)
    focus = Bool(False).tag(sync=True)
