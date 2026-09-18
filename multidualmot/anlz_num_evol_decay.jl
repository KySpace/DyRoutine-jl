# Included after anlz_num_evol.jl; use its statistics and curve selections directly.
fits_num_decay = Dict()
figs_num_decay = Dict()
for (idx_panel, panel) in enumerate(val_panel)
    fig = Figure(size=(plot_num_evol.size[1] + 310, plot_num_evol.size[2]))
    indices_panel = Any[Colon() for _ in name_stat]
    indices_panel[idx_axis_stat[key_panel]] = idx_panel
    reps_min, reps_max = extrema(@view n_rep_stat[indices_panel...])
    reps_used = reps_min == reps_max ? string(reps_min) : "$(reps_min)–$(reps_max)"
    ax = Axis(fig[1, 1]; xlabel=plot_num_evol.xlabel(panel), ylabel=plot_num_evol.ylabel,
        title=plot_num_evol.title(tag_head, panel, reps_used),
        yscale=log10, dualmot_axis_kwargs(; log_y=true)...)
    labels_fit = fig[1, 2] = GridLayout()
    curves_decay = map(enumerate(conditions_curve)) do (idx_condition, condition)
        indices = copy(indices_panel)
        for key in keys_curve
            indices[idx_axis_stat[key]] = findfirst(==(getproperty(condition, key)), getproperty(vars, key))
        end
        nums = vec(@view num_stat[indices...])
        stds = vec(@view std_num_stat[indices...])
        style = plot_num_evol.curve_style(condition)
        label = plot_num_evol.curve_label(condition)
        mask = isfinite.(nums) .& (nums .> 0)
        mask_error = mask .& isfinite.(stds) .& (nums .- stds .> 0)
        result = try
            fit_num_decay(val_x_plot, nums; fit_num_decay_config...)
        catch err
            @warn "$tag_head: decay fit failed" panel condition exception=err
            nothing
        end
        fits_num_decay[(panel, condition)] = result
        text_fit = if isnothing(result)
            "$label\nFit unavailable"
        else
            println("$tag_head / $panel / $label: ", replace(label_num_decay(result), '\n' => "; "))
            "$label\n$(label_num_decay(result))"
        end
        (; idx_condition, condition, nums, stds, style, label, mask, mask_error,
            result, text_fit)
    end
    for curve in curves_decay
        isnothing(curve.result) && continue
        ts = collect(range(minimum(val_x_plot), maximum(val_x_plot); length=400))
        lines!(ax, ts, fit_num_decay_config.model(ts, curve.result.params);
            curve.style.line_options...,
            xautolimits=false, yautolimits=false)
    end
    for curve in curves_decay
        mask_error = curve.mask_error
        errorbars!(ax, val_x_plot[mask_error], curve.nums[mask_error],
            curve.stds[mask_error]; curve.style.errorbar_options...)
        scatter!(ax, val_x_plot[curve.mask], curve.nums[curve.mask];
            curve.style.marker_options...,
            marker=curve.style.marker,
            label=curve.label)
        Label(labels_fit[curve.idx_condition, 1], curve.text_fit; color=curve.style.color,
            fontsize=11, halign=:left, justification=:left)
    end
    Label(labels_fit[length(conditions_curve) + 1, 1],
        "± approximate 1σ fit errors\n* parameter at bound; errors are local";
        fontsize=9, halign=:left, justification=:left)
    key_x in (:t_load, :t_hold) && set_time_minor_ticks!(ax)
    axislegend(ax; position=plot_num_evol.legend_position, DUALMOT_LEGEND_OPTIONS...)
    name_output = replace(plot_num_evol.filename(:log, panel), ".[log]." => ".[fit.log].")
    for format in plot_num_evol.formats
        save(joinpath(path_output, "$name_output.$format"), fig)
    end
    figs_num_decay[panel] = fig
end
