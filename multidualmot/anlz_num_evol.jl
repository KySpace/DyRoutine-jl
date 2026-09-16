# Included by a dataset-specific runner. Configured variable order is slowest to fastest.
num_data = read_atomnums(runinfo.files)
len_data = length(num_data)
n_variation = prod(length(getproperty(runinfo.vars, key)) for key in keys(runinfo.vars) if key != :rep)
len_data > 0 && rem(len_data, n_variation) == 0 ||
    throw(DimensionMismatch("$tag_head: $len_data atom numbers must be a positive integer multiple of $n_variation variations"))
n_rep = div(len_data, n_variation)
val_rep = 1:n_rep
vars = merge(runinfo.vars, (; rep=val_rep))
name, val_vars = keys(vars), values(vars)

name_acq = runinfo.var_order
n_dims_acq = map(key -> length(getproperty(vars, key)), name_acq)
num_acq = reshape(num_data, reverse(n_dims_acq)) |>
    data -> permutedims(data, reverse(1:length(n_dims_acq)))
num_fmt = permutedims(num_acq, indexin(collect(name), collect(name_acq)))

idx_rep_axis = findfirst(==(:rep), name)
num_stat = dropdims(mean(num_fmt; dims=idx_rep_axis); dims=idx_rep_axis)
std_num_stat = dropdims(std(num_fmt; dims=idx_rep_axis); dims=idx_rep_axis)
name_stat = Tuple(key for key in name if key != :rep)
n_rep == 1 && @warn "$tag_head: one repetition; sample standard deviations are undefined (NaN)"
println("$tag_head: $len_data samples / $n_variation variations = $n_rep repetitions")

key_x = plot_num_evol.key_x
key_panel = plot_num_evol.key_panel
keys_curve = keys(plot_num_evol.curves)
keys_plot = (key_x, key_panel, keys_curve...)
Set(keys_plot) == Set(name_stat) && length(keys_plot) == length(name_stat) ||
    throw(ArgumentError("$tag_head: plot axes must cover $(join(name_stat, ", ")) exactly once"))

val_x = getproperty(vars, key_x)
val_panel = getproperty(vars, key_panel)
vals_curve = map(keys_curve) do key
    values_all = getproperty(vars, key)
    selector = getproperty(plot_num_evol.curves, key)
    selector == :all ? values_all : [value for value in values_all if value in selector]
end
conditions_curve = [NamedTuple{keys_curve}(reverse(Tuple(values)))
    for values in Iterators.product(reverse(vals_curve)...)]
idx_axis_stat = Dict(key => idx for (idx, key) in enumerate(name_stat))

names_output = [plot_num_evol.filename(scale, panel)
    for panel in val_panel for scale in plot_num_evol.scales]
allunique(names_output) || throw(ArgumentError("$tag_head: plot filenames collide across panel values"))

figs_num_evol = Dict{Tuple{Any, Symbol}, Figure}()
for (idx_panel, panel) in enumerate(val_panel), scale in plot_num_evol.scales
    scale_num = scale == :lin ? plot_num_evol.scale_num_linear : 1.0
    fig = Figure(size=plot_num_evol.size)
    ax = Axis(fig[1, 1]; xlabel=plot_num_evol.xlabel,
        ylabel=scale == :lin ? "$(plot_num_evol.ylabel) (×10⁷)" : plot_num_evol.ylabel,
        title=plot_num_evol.title(tag_head, panel, n_rep),
        yscale=scale == :log ? log10 : identity,
        yminorticks=IntervalsBetween(5), yminorticksvisible=true)

    for condition in conditions_curve
        indices = Any[Colon() for _ in name_stat]
        indices[idx_axis_stat[key_panel]] = idx_panel
        for key in keys_curve
            values_axis = getproperty(vars, key)
            indices[idx_axis_stat[key]] = findfirst(==(getproperty(condition, key)), values_axis)
        end
        nums = vec(@view num_stat[indices...])
        stds = vec(@view std_num_stat[indices...])
        mask = scale == :log ? nums .> 0 : trues(length(nums))
        !all(mask) && @warn "$tag_head: omitting nonpositive log points" panel condition count=count(!, mask)
        nums_plot = ifelse.(mask, nums ./ scale_num, NaN)
        style = plot_num_evol.curve_style(condition)
        scatterlines!(ax, val_x, nums_plot; style...,
            label=plot_num_evol.curve_label(condition))

        mask_error = mask .& isfinite.(stds)
        if scale == :log
            mask_error .&= nums .- stds .> 0
            count_crossing = count(mask .& isfinite.(stds) .& (nums .- stds .<= 0))
            count_crossing > 0 && @warn "$tag_head: log error bars crossing zero omitted" panel condition count_crossing
        end
        errorbars!(ax, val_x[mask_error], nums[mask_error] ./ scale_num,
            stds[mask_error] ./ scale_num;
            color=style.color, whiskerwidth=7, linewidth=1.2)
    end
    axislegend(ax; position=:rt)
    name_output = plot_num_evol.filename(scale, panel)
    for format in plot_num_evol.formats
        save(joinpath(path_output, "$name_output.$format"), fig)
    end
    figs_num_evol[(panel, scale)] = fig
end
