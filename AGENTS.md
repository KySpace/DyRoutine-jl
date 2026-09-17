## Local Notes
- Julia launcher path to try first on this machine if the default doesn't work:
  `C:\Users\ky\AppData\Local\Microsoft\WindowsApps\julia.exe`
- When a task focuses on plotting or a late data-processing stage, prefer loading
  and processing the selected data once in a persistent Julia REPL, then rerun
  only the relevant later script section against the variables already in that
  REPL. Keep that REPL alive across related tests instead of repeating the
  expensive upstream processing.
- For excitation analysis, use `script/anlz_excitation_runner.jl` as the fresh
  analysis entry point, and use `script/anlz_excitation_rerun.jl` as the REPL
  debugging/replotting entry point when loading saved extraction or correlation
  caches. The split stage scripts are meant to be included from those entry
  points; visualization is split by purpose, for example correlation plots and
  extraction/solo plots can be toggled by choosing the corresponding vslz script.
- When developing excitation rerun workflows, new late-stage configurable
  variables may be added directly to `script/anlz_excitation_rerun.jl` for quick
  iteration. Once the change is stable, mirror those variables back into
  `script/anlz_excitation_runner.jl` so fresh analysis and reruns share the same
  configuration surface.

## Packaging
for now, some codes are not packaged, just included. Do package them in the future to decouple function calls from the environment by using `module`

## For Agents
Please summarize completed work in `Agent-Log.md`. Use one entry per finished task, and prefer this style:

```
- Request: Compact but faithful summary of what the user asked for, preserving important constraints and omitting incidental examples that could confuse future work.
  Brief: Compact but specific summary of what changed, major design choices, files/functions touched, and verification status or blockers.
```

While a task is still in progress, mostly edit the latest relevant log entry instead of adding new entries. Add a new entry only when a task is finished or when the user clearly starts a separate task. If you have doubt, ask before writing to `Agent-Log.md`.

## Code Style and Conventions
- Number-evolution analysis uses `multidualmot/anlz_num_evol.jl`; dataset runners
  declare the canonical variable specification, x-axis key, roots, and plotting
  configuration. Shared dual-isotope YAML/MAT readers and plot styles live in
  `multidualmot/dualmotcommons.jl`. Each pair folder's `config.yaml` is a
  top-level processing sequence; every processing entry has a unique `tag` and
  a nonempty `data` sequence, and each rectangular data block explicitly lists
  its `source` MAT filenames. A block's `vars` array includes `rep` and specifies
  the complete acquisition order; `rep` may be `:auto` and need not be
  outermost. Reshape in that order, then permute into the runner's canonical
  order. `anlz_num_evol.jl` handles one rectangular block per processing entry;
  use `anlz_cmpr_loadcfg.jl` for the three-block DCS/SCS loading comparison,
  which calculates each block independently and combines curves only after
  statistics. Its ratio uncertainty assumes independent DCS and SCS samples and
  propagates their sample standard deviations through `DCS/SCS`.
  DCS always contributes to reshaping; lifetime plots omit it, while MOT-loading
  plots include it. Plot colors encode `loadcfg`: DDM `#3a6ac4`, DIS `#b4452a`,
  DCS `#f4d33d`, and SCS `#55ca6f`. Marker shapes encode isotope: 160 diamond,
  161 up-triangle, 162 square, 163 down-triangle, and 164 circle. Use
  `Symbol("162")` (or `Symbol.(string.(values))`) for
  isotope symbols: Julia's `:162` evaluates to an integer, not a Symbol.
  Pair folder names must contain the configured ordered isotope pair as a
  complete numeric token. Processing descriptions such as `n-balanced` belong
  in the YAML `tag` rather than the folder name.
  Keep rejected number-evolution shots as `missing` entries in `num_fmt` so
  acquisition axes remain intact. In `anlz_num_evol.jl`, a shot is invalid when
  either fitted size is nonfinite or outside the runner's inclusive
  `bounds_sigmax_num` / `bounds_sigmay_num`, or when `atomnum` lies outside the
  inclusive range `0:num_max_num`. Compute means, sample deviations, and valid repetition
  counts by condition from the nonmissing repetitions.
  Configure display-unit conversion through the plot spec's `scale_x`; keep
  configured/acquired variable values unchanged and divide only `val_x_plot`.
- Prefer Julia functions with clear, narrow responsibilities. Keep processing logic in `src/` and script-specific orchestration in `script/`.
- Use `snake_case` for variables and functions. Function names usually read as verb + noun, such as `calc_dens_sum`, `crop_center`, `find_peak_position_moving`, or `set_axis_full`.
- Put the kind of quantity first, then attributes: examples include `wh_corner`, `smwh_peak`, `val_t`, `path_plot_peak`, and `dens_full_fmt`.
- Use compact domain prefixes consistently: `wh`/`hw` for image sizes, `smwh` for half-width/half-height crop spans, `xy` for pixel centers, `idx`/`ids` for indices, and `pos` for positions in index-like arrays.
- Keep experimental variation metadata named as `name`, `val_vars`, and variation/count variables such as `n_variation` or `n_dim_vars`. Earlier entries in `name` vary more slowly; later entries vary faster.
- Preserve array structure during processing when memory allows. Generate intermediate processed arrays first, keep variation axes explicit, then build visualizations from those results.
- Use `_fmt` for arrays whose outer dimensions correspond to experimental variables in `val_vars` order. Arrays ending in `_fmt` should generally have the same number of outer dimensions as there are variables, and those axes should align across related `_fmt` arrays. Store image-like or structured per-condition payloads as elements when possible, for example an n-d variable array whose entries are `h x w` images, rather than appending image height/width as extra `_fmt` axes. If a calculation reduces or replaces a variable axis, make that explicit in the variable name, metadata, or axis values.
- When reshaping flattened variation data, reshape so the fastest-varying variable is innermost before `permutedims` into the consumer-facing order.
- Prefer typed method signatures for public helpers and data-processing functions, especially `AbstractArray`, `AbstractMatrix`, `AbstractVector`, `Tuple{...}`, and `Real` bounds.
- Validate dimensions and user-facing numeric arguments early with `ArgumentError` or `DimensionMismatch`. Include the offending value and expected size when it helps debugging.
- Use functional data flow and Julia pipes for transformations, especially `|>` with short anonymous functions for `read`, `permutedims`, `mapslices`, `dropdims`, `reshape`, and plotting save steps. Prefer using @pipe and _ when the data is being used only once in between |>s (for fear of performance issues caused by 2 _ after @pipe)
- Prefer standard library or established package functions over custom implementations when they fit the task.
- Use `@view` for non-copying array regions and broadcasting for elementwise operations.
- For plotting, use mutating Makie-style helpers ending in `!` when they alter axes, figures, or layouts, such as `set_panel_solo_essn_2d!` and `draw_solo_essn_2d!`.
- Keep comments sparse and practical: short section comments for analysis stages are useful, while old exploratory plotting blocks may remain commented in scripts when they document active workflows.
- For now, scripts may `include(joinpath(@__DIR__, "..", "src", "..."))`; future packaging should move these relationships behind modules.
- In scripts, keep one-time workflow steps inline when they are short and their main purpose is to show the analysis flow. This includes simple conversion, packing/unpacking, selection masks, one-off cropping/centering passes, and final visualization steps. Prefer a named function only when the logic is reused, belongs in `src/`, or hiding it improves rather than obscures the script narrative.
