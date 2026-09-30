# Included by the CMOT/MOT lifetime runners after the number and fit figures.
figs_size_evol = Dict{Any,Figure}()
for (idx_panel, panel) in enumerate(val_panel)
    local indices_panel = Any[Colon() for _ in name_stat]
    indices_panel[idx_axis_stat[key_panel]] = idx_panel
    fig = Figure(size=(plot_num_evol.size[1], 255),
        fontsize=plot_num_evol.fontsize, figure_padding=1)
    ax_x = Axis(fig[1, 1];
        ylabel="σₓ (mm)",
        title=plot_num_evol.title(tag_head, panel, "$(minimum(n_rep_stat[indices_panel...]))–$(maximum(n_rep_stat[indices_panel...]))"),
        dualmot_axis_kwargs(; text_size=plot_num_evol.fontsize)...)
    ax_y = Axis(fig[2, 1];
        xlabel=plot_num_evol.xlabel(panel),
        ylabel="σᵧ (mm)",
        dualmot_axis_kwargs(; text_size=plot_num_evol.fontsize)...)
    linkxaxes!(ax_x, ax_y)
    rowgap!(fig.layout, 3)

    headers = Any["$(get(DUALMOT_CONFIG_VAR_NAMES, key_x, string(key_x))) (raw)",
        "$(key_x) plot (s)"]
    curves_size = map(conditions_curve) do condition
        indices = copy(indices_panel)
        for key in keys_curve
            indices[idx_axis_stat[key]] = findfirst(==(getproperty(condition, key)), getproperty(vars, key))
        end
        (; condition,
            sigmax=vec(@view sigmax_stat[indices...]),
            std_sigmax=vec(@view std_sigmax_stat[indices...]),
            sigmay=vec(@view sigmay_stat[indices...]),
            std_sigmay=vec(@view std_sigmay_stat[indices...]),
            n_reps=vec(@view n_rep_stat[indices...]),
            style=plot_num_evol.curve_style(condition),
            label=plot_num_evol.curve_label(condition))
    end
    for curve in curves_size
        append!(headers, ["σₓ mean [$(curve.label)] (mm)",
            "σₓ std [$(curve.label)] (mm)", "n_rep [$(curve.label)]",
            "σᵧ mean [$(curve.label)] (mm)", "σᵧ std [$(curve.label)] (mm)",
            "n_rep [$(curve.label)]", "source [$(curve.label)]"])
    end
    matrix = Matrix{Any}(undef, length(val_x) + 1, length(headers))
    matrix[1, :] .= headers
    for idx_x in eachindex(val_x)
        matrix[idx_x + 1, 1:2] .= (val_x[idx_x], val_x_plot[idx_x])
    end

    for (idx_curve, curve) in enumerate(curves_size)
        x = val_x_plot
        σx = curve.sigmax .* 1e3
        dσx = curve.std_sigmax .* 1e3
        σy = curve.sigmay .* 1e3
        dσy = curve.std_sigmay .* 1e3
        for (ax, values, errors, label) in
            ((ax_x, σx, dσx, curve.label), (ax_y, σy, dσy, nothing))
            valid = isfinite.(values)
            any(valid) || continue
            line_mask = valid .& isfinite.(x)
            count(line_mask) >= 2 && lines!(ax, x[line_mask], values[line_mask];
                curve.style.line_options...)
            error_mask = valid .& isfinite.(errors)
            if any(error_mask)
                marker_errorbars!(ax, x[error_mask], values[error_mask], errors[error_mask];
                    curve.style.marker_options..., curve.style.errorbar_options...,
                    marker=curve.style.marker, label)
                marker_only = valid .& .!error_mask
                any(marker_only) && scatter!(ax, x[marker_only], values[marker_only];
                    curve.style.marker_options..., marker=curve.style.marker)
            else
                scatter!(ax, x[valid], values[valid]; curve.style.marker_options...,
                    marker=curve.style.marker, label)
            end
        end
        col = 3 + 7 * (idx_curve - 1)
        for idx_x in eachindex(val_x)
            source_condition = merge(curve.condition,
                NamedTuple{(key_panel, key_x)}((panel, val_x[idx_x])))
            matrix[idx_x + 1, col:col + 6] .= (
                σx[idx_x], dσx[idx_x],
                curve.n_reps[idx_x] == 0 ? NaN : curve.n_reps[idx_x],
                σy[idx_x], dσy[idx_x],
                curve.n_reps[idx_x] == 0 ? NaN : curve.n_reps[idx_x],
                figure_xlsx_source(runinfo, source_condition))
        end
    end
    key_x in (:t_load, :t_hold) && set_time_minor_ticks!(ax_y)
    axislegend(ax_x; position=plot_num_evol.legend_position,
        DUALMOT_LEGEND_OPTIONS...)
    resize_to_layout!(fig)
    filename = replace(plot_num_evol.filename(:lin, panel),
        ".[lin]." => ".[size].")
    for format in plot_num_evol.formats
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
    sheet_name = figure_xlsx_sheet_name("size β$(@sprintf("%+.2f", Float64(panel))) $(runinfo.tag)")
    push!(get!(figure_data_sheets, runinfo.folder,
        Pair{String,Matrix{Any}}[]), sheet_name => matrix)
    figs_size_evol[panel] = fig
end
