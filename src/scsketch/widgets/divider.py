from anywidget import AnyWidget
from traitlets import Dict


class Div(AnyWidget):
    _esm = """
    function render({ model, el }) {
      const div = document.createElement("div");
      div.classList.add(
        'jupyter-widgets',
        'jupyter-scatter-div'
      );
      
      const update = () => {
        for (const [key, value] of Object.entries(model.get('style'))) {
          div.style[key] = value;
        }
      }
      
      model.on('change', update);

      update();

      el.appendChild(div);
    }
    export default { render };
    """

    style = Dict({}).tag(sync=True)
