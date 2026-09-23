# Included after anlz_num_evol.jl; use its statistics and curve selections directly.
fits_num_decay = Dict()
figs_num_decay = Dict()
for (idx_panel, panel) in enumerate(val_panel)
    local reps_min, reps_max, reps_used
    fig = Figure(size=plot_num_evol.fit_size, fontsize=plot_num_evol.fit_fontsize,
        figure_padding=1)
    indices_panel = Any[Colon() for _ in name_stat]
    indices_panel[idx_axis_stat[key_panel]] = idx_panel
    reps_min, reps_max = extrema(@view n_rep_stat[indices_panel...])
    reps_used = reps_min == reps_max ? string(reps_min) : "$(reps_min)–$(reps_max)"
    axis_size = (; width=plot_num_evol.frame_size[1],
        height=plot_num_evol.frame_size[2])
    ax_log = Axis(fig[1, 1]; ylabel=plot_num_evol.ylabel,
        title=plot_num_evol.title(tag_head, panel, reps_used),
        yscale=log10, xticklabelsvisible=false,
        axis_size...,
        dualmot_axis_kwargs(; log_y=true, text_size=plot_num_evol.fit_fontsize)...)
    ax_linear = Axis(fig[2, 1]; xlabel=plot_num_evol.xlabel(panel),
        ylabel=plot_num_evol.ylabel,
        axis_size...,
        dualmot_axis_kwargs(; text_size=plot_num_evol.fit_fontsize)...)
    linkxaxes!(ax_log, ax_linear)
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 8)
    labels_fit = fig[1:2, 2] = GridLayout()
    curves_decay = map(enumerate(conditions_curve)) do (idx_condition, condition)
        indices = copy(indices_panel)
        for key in keys_curve
            indices[idx_axis_stat[key]] = findfirst(==(getproperty(condition, key)), getproperty(vars, key))
        end
        nums = vec(@view num_stat[indices...])
        stds = vec(@view std_num_stat[indices...])
        style = plot_num_evol.curve_style(condition)
        label = plot_num_evol.curve_label(condition)
        mask_linear = isfinite.(nums)
        mask_log = mask_linear .& (nums .> 0)
        result = try
            fit_num_decay(val_x_plot, nums; fit_num_decay_config...)
        catch err
            @warn "$tag_head: decay fit failed" panel condition exception=err
            nothing
        end
        fits_num_decay[(panel, condition)] = result
        if !isnothing(result)
            p, e = result.params, result.errors
            tau = result.mode == :kappa ? Inf : p[2]
            std_tau = result.mode == :kappa ? NaN : e[2]
            kappa = result.mode == :tau ? 0.0 : p[result.mode == :full ? 3 : 2]
            std_kappa = result.mode == :tau ? NaN : e[result.mode == :full ? 3 : 2]
            push!(fit_records_num_decay, (
                pair=runinfo.folder,
                tag=runinfo.tag,
                panel=Float64(panel),
                loadcfg=condition.loadcfg,
                istp=condition.istp,
                mode=result.mode,
                n0=p[1],
                std_n0=e[1],
                tau,
                std_tau,
                kappa,
                std_kappa,
                at_bound=result.at_bound,
                n_points=count(result.mask),
            ))
        end
        text_fit = if isnothing(result)
            "$label\nFit unavailable"
        else
            number_unit = plot_num_evol.field_num == "pixsum" ? "pixsum" : "atom"
            fit_label = label_num_decay(result; number_unit)
            println("$tag_head / $panel / $label: ", replace(fit_label, '\n' => "; "))
            "$label\n$fit_label"
        end
        mask_selected = isnothing(result) ? trues(length(nums)) : result.mask_selected
        (; idx_condition, condition, nums, stds, style, label, mask_linear,
            mask_log, mask_selected, result, text_fit)
    end

    # Fitted curves stay underneath every marker and extend over all displayed times.
    for ax in (ax_log, ax_linear), curve in curves_decay
        isnothing(curve.result) && continue
        ts = collect(range(minimum(val_x_plot), maximum(val_x_plot); length=400))
        lines!(ax, ts, curve.result.model(ts, curve.result.params);
            curve.style.line_options...,
            xautolimits=false, yautolimits=false)
    end

    for (ax, is_log) in ((ax_log, true), (ax_linear, false)), curve in curves_decay
        mask_display = is_log ? curve.mask_log : curve.mask_linear
        mask_error = mask_display .& isfinite.(curve.stds)
        is_log && (mask_error .&= curve.nums .- curve.stds .> 0)
        mask_filled = mask_display .& curve.mask_selected
        mask_hollow = mask_display .& .!curve.mask_selected
        mask_filled_error = mask_filled .& mask_error
        mask_filled_only = mask_filled .& .!mask_error
        mask_hollow_error = mask_hollow .& mask_error
        mask_hollow_only = mask_hollow .& .!mask_error
        label_plot = ax === ax_log ? curve.label : nothing
        if any(mask_filled_error)
            marker_errorbars!(ax,
                val_x_plot[mask_filled_error], curve.nums[mask_filled_error],
                curve.stds[mask_filled_error];
                curve.style.marker_options...,
                curve.style.errorbar_options...,
                marker=curve.style.marker,
                label=label_plot)
            label_plot = nothing
        end
        if any(mask_filled_only)
            scatter!(ax, val_x_plot[mask_filled_only], curve.nums[mask_filled_only];
                curve.style.marker_options...,
                marker=curve.style.marker,
                label=label_plot)
            label_plot = nothing
        end
        hollow_options = merge(curve.style.marker_options, (; color=:transparent))
        if any(mask_hollow_error)
            marker_errorbars!(ax,
                val_x_plot[mask_hollow_error], curve.nums[mask_hollow_error],
                curve.stds[mask_hollow_error];
                hollow_options...,
                curve.style.errorbar_options...,
                marker=curve.style.marker,
                label=label_plot)
            label_plot = nothing
        end
        any(mask_hollow_only) && scatter!(ax,
            val_x_plot[mask_hollow_only], curve.nums[mask_hollow_only];
            hollow_options...,
            marker=curve.style.marker,
            label=label_plot)
    end

    for curve in curves_decay
        Label(labels_fit[curve.idx_condition, 1], curve.text_fit;
            color=curve.style.color, fontsize=8,
            halign=:left, justification=:left)
    end
    Label(labels_fit[length(conditions_curve) + 1, 1],
        "± approximate 1σ fit errors\n* parameter at bound; errors are local";
        fontsize=7, halign=:left, justification=:left)
    key_x in (:t_load, :t_hold) && begin
        set_time_minor_ticks!(ax_log)
        set_time_minor_ticks!(ax_linear)
    end
    axislegend(ax_log; position=plot_num_evol.legend_position,
        DUALMOT_LEGEND_OPTIONS...)
    resize_to_layout!(fig)
    name_output = replace(plot_num_evol.filename(:log, panel),
        ".[log]." => ".[fit.log].")
    for format in plot_num_evol.formats
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$name_output.$format"), fig; save_options...)
    end
    figs_num_decay[panel] = fig
end
