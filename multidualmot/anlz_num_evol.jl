# Included by a dataset-specific runner. Configured variable order is slowest to fastest.
length(runinfo.data) == 1 ||
    throw(ArgumentError("$tag_head: anlz_num_evol.jl requires exactly one rectangular data entry"))
runinfo_data = only(runinfo.data)
shot_data = read_cres_fields(runinfo_data.files, ("atomnum", "sigmax", "sigmay"))
num_data, sigmax_data, sigmay_data = shot_data
all(isfinite, num_data) || throw(ArgumentError("$tag_head: atomnum must contain only finite values"))
len_data = length(num_data)
n_variation = prod(length(getproperty(runinfo_data.vars, key)) for key in keys(runinfo_data.vars) if key != :rep)
len_data > 0 && rem(len_data, n_variation) == 0 ||
    throw(DimensionMismatch("$tag_head: $len_data atom numbers must be a positive integer multiple of $n_variation variations"))
n_rep = div(len_data, n_variation)
val_rep = 1:n_rep
vars = merge(runinfo_data.vars, (; rep=val_rep))
name, val_vars = keys(vars), values(vars)

name_acq = runinfo_data.var_order
n_dims_acq = map(key -> length(getproperty(vars, key)), name_acq)
for (name_bound, bounds) in
    (("sigmax", bounds_sigmax_num), ("sigmay", bounds_sigmay_num))
    lower, upper = bounds
    isfinite(lower) && !isnan(upper) && lower <= upper ||
        throw(ArgumentError("$name_bound bounds must satisfy finite lower <= upper, got $bounds"))
end
mask_sigmax = isfinite.(sigmax_data) .&
    (sigmax_data .>= bounds_sigmax_num[1]) .& (sigmax_data .<= bounds_sigmax_num[2])
mask_sigmay = isfinite.(sigmay_data) .&
    (sigmay_data .>= bounds_sigmay_num[1]) .& (sigmay_data .<= bounds_sigmay_num[2])
mask_size = mask_sigmax .& mask_sigmay
mask_num_low = num_data .>= 0
mask_num_high = num_data .<= num_max_num
mask_num = mask_num_low .& mask_num_high
mask_valid = mask_size .& mask_num
num_data_masked = Vector{Union{Missing, Float64}}(num_data)
num_data_masked[.!mask_valid] .= missing
n_masked_size = count(!, mask_size)
n_masked_sigmax = count(!, mask_sigmax)
n_masked_sigmay = count(!, mask_sigmay)
n_masked_num_low = count(!, mask_num_low)
n_masked_num_high = count(!, mask_num_high)
n_masked_total = count(!, mask_valid)

num_acq = reshape(num_data_masked, reverse(n_dims_acq)) |>
    data -> permutedims(data, reverse(1:length(n_dims_acq)))
num_fmt = permutedims(num_acq, indexin(collect(name), collect(name_acq)))

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
println("$tag_head: $len_data samples / $n_variation variations = $n_rep repetitions; " *
    "$n_masked_total rejected ($n_masked_size by size: $n_masked_sigmax sigmax, " *
    "$n_masked_sigmay sigmay; $n_masked_num_low below zero, " *
    "$n_masked_num_high above $num_max_num)")

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
    # Stable ordering preserves isotope order while keeping DIS/SCS visible at coincident points.
    sort!(conditions_curve;
        by=condition -> getproperty(condition, :loadcfg) in (:DIS, :SCS) ? 1 : 0,
        alg=Base.Sort.MergeSort)
end
idx_axis_stat = Dict(key => idx for (idx, key) in enumerate(name_stat))

names_output = [plot_num_evol.filename(scale, panel)
    for panel in val_panel for scale in plot_num_evol.scales]
allunique(names_output) || throw(ArgumentError("$tag_head: plot filenames collide across panel values"))

figs_num_evol = Dict{Tuple{Any, Symbol}, Figure}()
for (idx_panel, panel) in enumerate(val_panel), scale in plot_num_evol.scales
    scale_num = scale == :lin ? plot_num_evol.scale_num_linear : 1.0
    fig = Figure(size=plot_num_evol.size, fontsize=plot_num_evol.fontsize,
        figure_padding=1)
    indices_panel = Any[Colon() for _ in name_stat]
    isnothing(key_panel) || (indices_panel[idx_axis_stat[key_panel]] = idx_panel)
    reps_panel = vec(@view n_rep_stat[indices_panel...])
    reps_min, reps_max = extrema(reps_panel)
    reps_used = reps_min == reps_max ? string(reps_min) : "$(reps_min)–$(reps_max)"
    ax = Axis(fig[1, 1]; xlabel=plot_num_evol.xlabel(panel),
        ylabel=scale == :lin ? "$(plot_num_evol.ylabel) (×10⁷)" : plot_num_evol.ylabel,
        title=plot_num_evol.title(tag_head, panel, reps_used),
        yscale=scale == :log ? log10 : identity,
        aspect=AxisAspect(4 / 3),
        dualmot_axis_kwargs(; log_y=scale == :log)...)

    curves_plot = map(conditions_curve) do condition
        indices = Any[Colon() for _ in name_stat]
        isnothing(key_panel) || (indices[idx_axis_stat[key_panel]] = idx_panel)
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
        idx_istp = :istp in keys_curve ?
            findfirst(==(condition.istp), vars.istp) : nothing
        val_x_curve = plot_num_evol.transform_x(
            val_x_plot, condition, panel, idx_istp)
        length(val_x_curve) == length(val_x_plot) ||
            throw(DimensionMismatch("$tag_head: transformed x length $(length(val_x_curve)) must equal $(length(val_x_plot))"))
        all(isfinite, val_x_curve) ||
            throw(ArgumentError("$tag_head: transformed x values must be finite for $condition"))
        xautolimits = plot_num_evol.xautolimits(condition, panel)
        xautolimits isa Bool ||
            throw(ArgumentError("$tag_head: xautolimits must return Bool for $condition"))
        yautolimits = plot_num_evol.yautolimits(condition, panel)
        yautolimits isa Bool ||
            throw(ArgumentError("$tag_head: yautolimits must return Bool for $condition"))
        mask_error = mask .& isfinite.(stds)
        if scale == :log
            mask_error .&= nums .- stds .> 0
            count_crossing = count(mask .& isfinite.(stds) .& (nums .- stds .<= 0))
            count_crossing > 0 && @warn "$tag_head: log error bars crossing zero omitted" panel condition count_crossing
        end
        (; condition, nums, stds, nums_plot, val_x_curve, mask_error, style,
            xautolimits, yautolimits)
    end
    for curve in curves_plot
        lines!(ax, curve.val_x_curve, curve.nums_plot;
            curve.style.line_options...,
            xautolimits=curve.xautolimits, yautolimits=curve.yautolimits)
    end
    for curve in curves_plot
        mask_error = curve.mask_error
        mask_marker = isfinite.(curve.nums_plot)
        label = plot_num_evol.curve_label(curve.condition)
        if any(mask_error)
            marker_errorbars!(ax,
                curve.val_x_curve[mask_error], curve.nums_plot[mask_error],
                curve.stds[mask_error] ./ scale_num;
                curve.style.marker_options...,
                curve.style.errorbar_options...,
                marker=curve.style.marker,
                xautolimits=curve.xautolimits,
                yautolimits=curve.yautolimits,
                label)
            mask_marker_only = mask_marker .& .!mask_error
            any(mask_marker_only) && scatter!(ax,
                curve.val_x_curve[mask_marker_only], curve.nums_plot[mask_marker_only];
                curve.style.marker_options...,
                marker=curve.style.marker,
                xautolimits=curve.xautolimits,
                yautolimits=curve.yautolimits,
                label=nothing)
        else
            scatter!(ax, curve.val_x_curve[mask_marker], curve.nums_plot[mask_marker];
                curve.style.marker_options...,
                marker=curve.style.marker,
                xautolimits=curve.xautolimits,
                yautolimits=curve.yautolimits,
                label)
        end
    end
    key_x in (:t_load, :t_hold) && set_time_minor_ticks!(ax)
    axislegend(ax; position=plot_num_evol.legend_position, DUALMOT_LEGEND_OPTIONS...)
    name_output = plot_num_evol.filename(scale, panel)
    for format in plot_num_evol.formats
        save(joinpath(path_output, "$name_output.$format"), fig)
    end
    figs_num_evol[(panel, scale)] = fig
end
