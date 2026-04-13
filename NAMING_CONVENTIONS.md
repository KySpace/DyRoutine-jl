# Naming Conventions

This project uses the following names for experimental variation metadata:

- `name`: the names of the varying experimental variables.
- `val`: the possible values for each variable, aligned by index with `name`.
- `variation`: the full Cartesian variation across all variables.

Current example:

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
