# Included after anlz_num_evol.jl and anlz_num_evol_output.jl for an ODT scan.
fits_num_gaussian = Dict()
reps_gaussian = vec(n_rep_stat)
reps_min_gaussian, reps_max_gaussian = extrema(reps_gaussian)
reps_used_gaussian = reps_min_gaussian == reps_max_gaussian ?
    string(reps_min_gaussian) : "$(reps_min_gaussian)–$(reps_max_gaussian)"
fig_gaussian = Figure(size=plot_num_evol.size, fontsize=plot_num_evol.fontsize,
    figure_padding=1)
axis_options_gaussian = merge(
    dualmot_axis_kwargs(; text_size=plot_num_evol.fontsize),
    plot_num_evol.axis_options)
ax_gaussian = Axis(fig_gaussian[1, 1];
    xlabel=plot_num_evol.xlabel(nothing),
    ylabel=rich(plot_num_evol.ylabel, " (10⁷)"),
    title=plot_num_evol.title(tag_head, nothing, reps_used_gaussian),
    axis_options_gaussian...)
scale_num_gaussian = plot_num_evol.scale_num_linear
records_gaussian = NamedTuple[]
rows_gaussian = Vector{Vector{Any}}()

for condition in conditions_curve
    indices = Any[Colon() for _ in name_stat]
    for key in keys_curve
        indices[idx_axis_stat[key]] = findfirst(==(getproperty(condition, key)),
            getproperty(vars, key))
    end
    nums = vec(@view num_stat[indices...])
    stds = vec(@view std_num_stat[indices...])
    n_reps = vec(@view n_rep_stat[indices...])
    idx_istp = :istp in keys_curve ? findfirst(==(condition.istp), vars.istp) : nothing
    x_curve = plot_num_evol.transform_x(val_x_plot, condition, nothing, idx_istp)
    fit_result = fit_num_gaussian(x_curve, nums)
    fits_num_gaussian[condition] = fit_result
    style = plot_num_evol.curve_style(condition)
    fit_x = x_curve[fit_result.mask]
    fit_x_sorted = collect(range(minimum(fit_x), maximum(fit_x); length=400))
    fit_nums = fit_result.model(fit_x_sorted, fit_result.params) ./ scale_num_gaussian
    lines!(ax_gaussian, fit_x_sorted, fit_nums;
        style.line_options..., xautolimits=false, yautolimits=false)

    mask_data = fit_result.mask
    mask_error = mask_data .& isfinite.(stds) .& (stds .>= 0)
    mask_marker = mask_data
    label = plot_num_evol.curve_label(condition)
    if any(mask_error)
        marker_errorbars!(ax_gaussian, x_curve[mask_error],
            nums[mask_error] ./ scale_num_gaussian,
            stds[mask_error] ./ scale_num_gaussian;
            style.marker_options..., style.errorbar_options...,
            marker=style.marker, label)
        mask_marker_only = mask_marker .& .!mask_error
        any(mask_marker_only) && scatter!(ax_gaussian,
            x_curve[mask_marker_only], nums[mask_marker_only] ./ scale_num_gaussian;
            style.marker_options..., marker=style.marker, label=nothing)
    else
        scatter!(ax_gaussian, x_curve[mask_marker], nums[mask_marker] ./ scale_num_gaussian;
            style.marker_options..., marker=style.marker, label)
    end

    baseline, amplitude, center, sigma = fit_result.params
    err_baseline, err_amplitude, err_center, err_sigma = fit_result.errors
    source = figure_xlsx_source(runinfo, condition)
    record = (
        pair=runinfo.folder,
        tag=runinfo.tag,
        loadcfg=condition.loadcfg,
        istp=condition.istp,
        baseline,
        std_baseline=err_baseline,
        amplitude,
        std_amplitude=err_amplitude,
        center_G=center,
        std_center_G=err_center,
        sigma_G=sigma,
        std_sigma_G=err_sigma,
        at_bound=fit_result.at_bound,
        n_points=fit_result.n_points,
        residual_sse=fit_result.residual_sse,
        residual_rms=fit_result.residual_rms,
        source,
    )
    push!(records_gaussian, record)
    push!(rows_gaussian, Any[record.tag, string(record.loadcfg), string(record.istp),
        record.baseline, record.std_baseline, record.amplitude, record.std_amplitude,
        record.center_G, record.std_center_G, record.sigma_G, record.std_sigma_G,
        record.at_bound, record.n_points, record.residual_sse, record.residual_rms,
        record.source])
end

axislegend(ax_gaussian; position=plot_num_evol.legend_position,
    DUALMOT_LEGEND_OPTIONS...)
filename_gaussian = replace(plot_num_evol.filename(:lin, nothing),
    "].[lin]." => "].[fit.lin].")
for format in plot_num_evol.formats
    save_options = format == "png" ? (; px_per_unit=4) : (;)
    save(joinpath(path_output, "$filename_gaussian.$format"), fig_gaussian;
        save_options...)
end

headers_gaussian = Any["tag", "loadcfg", "istp", "baseline", "std baseline",
    "amplitude", "std amplitude", "center (G)", "std center (G)", "sigma (G)",
    "std sigma (G)", "at bound", "n points", "residual SSE", "residual RMS", "source"]
matrix_gaussian = Matrix{Any}(undef, length(rows_gaussian) + 1, length(headers_gaussian))
matrix_gaussian[1, :] .= headers_gaussian
for (idx_row, row) in enumerate(rows_gaussian)
    matrix_gaussian[idx_row + 1, :] .= row
end
sheet_gaussian = num_evol_xlsx_sheet_name(runinfo, nothing, nothing; fit=true)
push!(get!(figure_data_sheets, runinfo.folder,
    Pair{String,Matrix{Any}}[]), sheet_gaussian => matrix_gaussian)
append!(fit_records_num_gaussian, records_gaussian)
