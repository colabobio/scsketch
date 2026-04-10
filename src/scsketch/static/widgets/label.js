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
