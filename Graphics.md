## About backend
As per suggested, we can switch between backends, mainly GLMakie and CairoMakie. The idea is that, 
- in plotting packages, we use both package. 
- For most plots, GLMakie backend is activated so that it's shown interactively. 
- When printing, use `save("myplot.pdf", fig; backend=CairoMakie)` for Cairo backend.
- For very large prints, do not display interactively.
