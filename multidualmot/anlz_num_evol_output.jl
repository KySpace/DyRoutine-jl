if isnothing(plot_num_evol_target)
    names_output = [plot_num_evol.filename(scale, panel)
        for panel in val_panel for scale in plot_num_evol.scales]
    allunique(names_output) || throw(ArgumentError("$tag_head: plot filenames collide across panel values"))
end

figs_num_evol = Dict{Tuple{Any, Symbol}, Figure}()
jobs_num_evol = isnothing(plot_num_evol_target) ?
    [(idx_panel, panel, scale) for (idx_panel, panel) in enumerate(val_panel)
        for scale in plot_num_evol.scales] :
    [(findfirst(==(plot_num_evol_target.panel), val_panel),
        plot_num_evol_target.panel, plot_num_evol_target.scale)]
for (idx_panel, panel, scale) in jobs_num_evol
    local reps_min, reps_max, reps_used
    isnothing(idx_panel) && throw(ArgumentError("$tag_head: requested panel $panel is unavailable"))
    scale_num = scale == :lin ? plot_num_evol.scale_num_linear : 1.0
    fig = isnothing(plot_num_evol_target) ?
        Figure(size=plot_num_evol.size, fontsize=plot_num_evol.fontsize,
            figure_padding=1) : plot_num_evol_target.fig
    slot = isnothing(plot_num_evol_target) ? fig[1, 1] : plot_num_evol_target.slot
    frame_options = isnothing(plot_num_evol_target) ?
        (; width=plot_num_evol.frame_size[1], height=plot_num_evol.frame_size[2]) :
        (; aspect=AxisAspect(4 / 3))
    indices_panel = Any[Colon() for _ in name_stat]
    isnothing(key_panel) || (indices_panel[idx_axis_stat[key_panel]] = idx_panel)
    reps_panel = vec(@view n_rep_stat[indices_panel...])
    reps_min, reps_max = extrema(reps_panel)
    reps_used = reps_min == reps_max ? string(reps_min) : "$(reps_min)–$(reps_max)"
    axis_options = dualmot_axis_kwargs(; log_y=scale == :log,
        compact_spacing=!isnothing(plot_num_evol_target))
    !isnothing(plot_num_evol_target) &&
        (axis_options = merge(axis_options, plot_num_evol_target.axis_options))
    ax = Axis(slot; xlabel=plot_num_evol.xlabel(panel),
        ylabel=scale == :lin ? "$(plot_num_evol.ylabel) (10⁷)" : plot_num_evol.ylabel,
        title=isnothing(plot_num_evol_target) ?
            plot_num_evol.title(tag_head, panel, reps_used) : "",
        yscale=scale == :log ? log10 : identity,
        frame_options...,
        axis_options...)

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
    if isnothing(plot_num_evol_target)
        key_x in (:t_load, :t_hold) && set_time_minor_ticks!(ax)
        if scale == :lin && !isnothing(plot_num_evol.limits_linear) && iszero(panel)
            xlims!(ax, plot_num_evol.limits_linear.x...)
            ylims!(ax, plot_num_evol.limits_linear.y...)
        end
        axislegend(ax; position=plot_num_evol.legend_position, DUALMOT_LEGEND_OPTIONS...)
    else
        if !isnothing(plot_num_evol_target.limits)
            xlimits = plot_num_evol_target.limits.x
            ylimits = plot_num_evol_target.limits.y
            isnothing(xlimits) || xlims!(ax, xlimits...)
            isnothing(ylimits) || ylims!(ax, ylimits...)
        end
        plot_num_evol_target.show_legend &&
            axislegend(ax; position=plot_num_evol.legend_position, DUALMOT_LEGEND_OPTIONS...)
    end
    if isnothing(plot_num_evol_target)
        name_output = plot_num_evol.filename(scale, panel)
        for format in plot_num_evol.formats
            save(joinpath(path_output, "$name_output.$format"), fig)
        end
    end
    figs_num_evol[(panel, scale)] = fig
end
