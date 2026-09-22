# Included by a dataset-specific runner. Sources are aligned by their configured
# variable values and appended along rep; opted-in scans may pad a truncated tail.
blocks_num_evol = map(enumerate(runinfo.data)) do (idx_data, runinfo_data)
    calc_num_evol_block(runinfo_data;
        label="$tag_head data[$idx_data]",
        bounds_sigmax_num,
        bounds_sigmay_num,
        num_max_num,
        allow_partial_rep=plot_num_evol.allow_partial_rep,
    )
end
combined_num_evol = combine_num_evol_blocks(blocks_num_evol; label=tag_head)
vars = combined_num_evol.vars
name, val_vars = keys(vars), values(vars)
n_rep = combined_num_evol.n_rep
val_rep = vars.rep
num_fmt = combined_num_evol.num_fmt

idx_rep_axis = findfirst(==(:rep), name)
num_stat = dropdims(mapslices(num_fmt; dims=idx_rep_axis) do values
    valid = collect(skipmissing(vec(values)))
    isempty(valid) ? NaN : mean(valid)
end; dims=idx_rep_axis)
std_num_stat = dropdims(mapslices(num_fmt; dims=idx_rep_axis) do values
    valid = collect(skipmissing(vec(values)))
    length(valid) < 2 ? NaN : std(valid)
end; dims=idx_rep_axis)
n_rep_stat = dropdims(sum(.!ismissing.(num_fmt); dims=idx_rep_axis); dims=idx_rep_axis)
name_stat = Tuple(key for key in name if key != :rep)
n_rep == 1 && @warn "$tag_head: one repetition; sample standard deviations are undefined (NaN)"
length(blocks_num_evol) > 1 && println("$tag_head: combined $(length(blocks_num_evol)) data blocks " *
    "into $(size(num_fmt)) with $n_rep total repetition slots")

key_x = plot_num_evol.key_x
key_panel = plot_num_evol.key_panel
keys_curve = keys(plot_num_evol.curves)
keys_plot = isnothing(key_panel) ? (key_x, keys_curve...) : (key_x, key_panel, keys_curve...)
Set(keys_plot) == Set(name_stat) && length(keys_plot) == length(name_stat) ||
    throw(ArgumentError("$tag_head: plot axes must cover $(join(name_stat, ", ")) exactly once"))

val_x = getproperty(vars, key_x)
val_x_plot = val_x ./ plot_num_evol.scale_x
val_panel = isnothing(key_panel) ? (nothing,) : getproperty(vars, key_panel)
vals_curve = map(keys_curve) do key
    values_all = getproperty(vars, key)
    selector = getproperty(plot_num_evol.curves, key)
    selector == :all ? values_all : [value for value in values_all if value in selector]
end
conditions_curve = vec([NamedTuple{keys_curve}(reverse(Tuple(values)))
    for values in Iterators.product(reverse(vals_curve)...)])
if :loadcfg in keys_curve
    order_loadcfg = (:SCS, :DCS, :DIS, :DDM)
    sort!(conditions_curve;
        by=condition -> findfirst(==(getproperty(condition, :loadcfg)), order_loadcfg),
        alg=Base.Sort.MergeSort)
end
idx_axis_stat = Dict(key => idx for (idx, key) in enumerate(name_stat))
