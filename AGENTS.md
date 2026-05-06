## Local Notes
- Julia launcher path to try first on this machine if the default doesn't work:
  `C:\Users\ky\AppData\Local\Microsoft\WindowsApps\julia.exe`

## Packaging
for now, some codes are not packaged, just included. Do package them in the future to decouple function calls from the environment by using `module`

## Code Style and Conventions
- Prefer Julia functions with clear, narrow responsibilities. Keep processing logic in `src/` and script-specific orchestration in `script/`.
- Use `snake_case` for variables and functions. Function names usually read as verb + noun, such as `calc_dens_sum`, `crop_center`, `find_peak_position_moving`, or `set_axis_full`.
- Put the kind of quantity first, then attributes: examples include `wh_corner`, `smwh_peak`, `val_t`, `path_plot_peak`, and `dens_full_fmt`.
- Use compact domain prefixes consistently: `wh`/`hw` for image sizes, `smwh` for half-width/half-height crop spans, `xy` for pixel centers, `idx`/`ids` for indices, and `pos` for positions in index-like arrays.
- Keep experimental variation metadata named as `name`, `val`, and variation/count variables such as `n_variation` or `n_dim_vars`. Earlier entries in `name` vary more slowly; later entries vary faster.
- Preserve array structure during processing when memory allows. Generate intermediate processed arrays first, keep variation axes explicit, then build visualizations from those results.
- When reshaping flattened variation data, reshape so the fastest-varying variable is innermost before `permutedims` into the consumer-facing order.
- Prefer typed method signatures for public helpers and data-processing functions, especially `AbstractArray`, `AbstractMatrix`, `AbstractVector`, `Tuple{...}`, and `Real` bounds.
- Validate dimensions and user-facing numeric arguments early with `ArgumentError` or `DimensionMismatch`. Include the offending value and expected size when it helps debugging.
- Use functional data flow and Julia pipes for transformations, especially `|>` with short anonymous functions for `read`, `permutedims`, `mapslices`, `dropdims`, `reshape`, and plotting save steps.
- Prefer standard library or established package functions over custom implementations when they fit the task.
- Use `@view` for non-copying array regions and broadcasting for elementwise operations.
- For plotting, use mutating Makie-style helpers ending in `!` when they alter axes, figures, or layouts, such as `set_panel_solo_essn_2d!` and `draw_solo_essn_2d!`.
- Keep comments sparse and practical: short section comments for analysis stages are useful, while old exploratory plotting blocks may remain commented in scripts when they document active workflows.
- For now, scripts may `include(joinpath(@__DIR__, "..", "src", "..."))`; future packaging should move these relationships behind modules.
