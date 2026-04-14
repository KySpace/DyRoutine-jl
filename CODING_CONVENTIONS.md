# For Agents
Ask before looking for references. Ask before if a new technical stack is required and you cannot access for now (like loading a new library), especially due to permission issues.
# Naming Conventions

This project uses the following names for experimental variation metadata:

- snake_casing for local variables or function names. For function names: tend to use verb + noun structure.
- When naming, usually put the type of the thing first, then attributes, like `wh_crop`, `val_t_hold`.
- use `idx` or `ids` for index or indices, `pos` for positions in an array of indices (second order or higher order indices).
- `name`: the names of the varying experimental variables.
- `val`: the possible values for each variable, aligned by index with `name`.
- `variation`: the full Cartesian variation across all variables.
- use `wh` or `hw` to refer to image sizes, rather than `size` or `sz` to avoid confusion.

# Coding style
functional style, prefers piping operations when processing data.

# Data processing
As long as the processing is not too memory consuming:
Generate the data first, keep the structure of the data as much as possible (e.g. per-shot data should have the same dimension as the original `dens` other than the image dimensions), then generate visualizations based on the results of the data.

# Current example:

```julia
name = ["repeat", "t_hold", "istp"]
val = (
    1:3,
    6:2:200,
    [5, 0],
)
variation = 3 * 98 * 2
```

Ordering rule:

- Earlier entries in `name` vary more slowly.
- Later entries in `name` vary more quickly.
- In the current example, `repeat` is outermost and `istp` is fastest.

When a raw data axis stores the flattened combined variation, reshape it so that
the fastest-varying variable is the innermost variation axis before permuting to
the desired consumer-facing order.
