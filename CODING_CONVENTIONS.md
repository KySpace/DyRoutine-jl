# For Agents
Ask before looking for references. Ask before if a new technical stack is required and you cannot access for now (like loading a new library), especially due to permission issues.
Prefer using well-established functions than custom defined ones.
# Naming Conventions

This project uses the following names for experimental variation metadata:

- snake_casing for local variables or function names. For function names: tend to use verb + noun structure.
- When naming, usually put the type of the thing first, then attributes, like `wh_crop`, `val_t_hold`.
- use `idx` or `ids` for index or indices, `pos` for positions in an array of indices (second order or higher order indices).
- Prefer recording experimental variables in a named tuple, for example `vars=(; rep, t_hold, istp)` or `vars=(IB=5.378, rep=1:5, bias=0.1:0.05:0.6, t_hold=6:5:56, istp)`.
- Put commonly identical variable values outside the `runinfos` list and use named-tuple shorthand inside each `runinfo`. If all variable axes are shared, write `vars` once and attach it to each `runinfo` before calling the per-run script.
- `name` or `name_dims`: the names of the varying experimental variables, usually inferred from `propertynames(val)` instead of written separately.
- `val`: the possible values for each variable, preferably as a named tuple whose field order is the axis order.
- `variation`: the full Cartesian variation across all variables.
- use `wh` or `hw` to refer to image sizes, rather than `size` or `sz` to avoid confusion.
- Values named `istp`, `rep`, `t_hold`, etc. may be used as shared axis values. Avoid reusing these names as loop counters in top-level script scope; use names like `idx_istp` or `val_istp` instead.

# Coding style
functional style, prefers piping operations when processing data.

# Data processing
As long as the processing is not too memory consuming:
Generate the data first, keep the structure of the data as much as possible (e.g. per-solo data should have the same dimension as the original `dens` other than the image dimensions), then generate visualizations based on the results of the data.

For arrays ending in `_fmt`:

- The outer dimensions should correspond to experimental variables in `val` order.
- Related `_fmt` arrays should have aligned axes and generally the same number of outer dimensions.
- For image-like payloads, prefer an n-d array over variables whose entries are `h x w` images, rather than appending image height and width as extra `_fmt` axes.
- If a calculation reduces or replaces a variable axis, make that explicit in the variable name, metadata, or axis values. For example, a statistic over `rep` should not be treated as if the `rep` axis were still raw repeats.
- When averaging or transforming one variable axis, prefer clear axis-aware operations such as `mapslices(...; dims=axis)` when they express intent better than manual `CartesianIndex` slicing.

# Current example:

```julia
rep = 1:3
t_hold = 6:2:200
istp = ["162", "164"]
vars = (; rep, t_hold, istp)
val = map(collect, vars)
name_dims = propertynames(val)
n_dim_vars = map(length, val)
n_variation = prod(n_dim_vars)
```

Ordering rule:

- Earlier entries in `name_dims` / `val` vary more slowly.
- Later entries in `name_dims` / `val` vary more quickly.
- In the current example, `rep` is outermost and `istp` is fastest.

When a raw data axis stores the flattened combined variation, reshape it so that
the fastest-varying variable is the innermost variation axis before permuting to
the desired consumer-facing order.

For loaded density data, a typical formatted result is:

```julia
dens_full_fmt  # size == n_dim_vars
first(dens_full_fmt)  # one cropped h x w image
```

When reducing over `istp` while preserving a singleton final axis, a clear pattern is:

```julia
xy_peak_duet = r.dens_full_fmt |>
               ds -> mapslices(f, ds; dims=ndims(ds))
```
