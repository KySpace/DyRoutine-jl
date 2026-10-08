length(runinfo.data) == 3 ||
    throw(ArgumentError("$tag_head: expected exactly three rectangular data blocks"))
stats_data = [calc_num_evol_block(data;
    label="$tag_head data[$idx]", bounds_sigmax_num, bounds_sigmay_num, num_max_num)
    for (idx, data) in enumerate(runinfo.data)]

val_istp = Symbol.(split(runinfo.folder, "-"))
val_t_load = stats_data[1].vars.t_load
all(stats -> stats.vars.t_load == val_t_load, stats_data) ||
    throw(ArgumentError("$tag_head: all data blocks must use the same t_load values"))
isfinite(plot_cmpr_loadcfg.γ_active) && plot_cmpr_loadcfg.γ_active > 0 ||
    throw(ArgumentError("γ_active must be finite and positive, got $(plot_cmpr_loadcfg.γ_active)"))
val_t_load_plot = val_t_load .* plot_cmpr_loadcfg.γ_active
mask_t_load_target = isnothing(plot_cmpr_loadcfg_target) ? trues(length(val_t_load)) :
    plot_cmpr_loadcfg_target.mask_x(val_t_load)
length(mask_t_load_target) == length(val_t_load) ||
    throw(DimensionMismatch("$tag_head: target loading-time mask length $(length(mask_t_load_target)) must equal $(length(val_t_load))"))

curves_num = Dict{Tuple{Symbol, Symbol}, NamedTuple}()
for stats in stats_data
    idx_axis = Dict(key => idx for (idx, key) in enumerate(stats.name_stat))
    for (idx_loadcfg, loadcfg) in enumerate(stats.vars.loadcfg),
        (idx_istp, istp) in enumerate(stats.vars.istp)
        key = (loadcfg, istp)
        haskey(curves_num, key) && throw(ArgumentError("$tag_head: duplicate curve $key"))
        indices = Any[Colon() for _ in stats.name_stat]
        indices[idx_axis[:loadcfg]] = idx_loadcfg
        indices[idx_axis[:istp]] = idx_istp
        curves_num[key] = (
            nums=vec(@view stats.num_stat[indices...]),
            stds=vec(@view stats.std_num_stat[indices...]),
            n_reps=vec(@view stats.n_rep_stat[indices...]),
        )
    end
end
expected_curves = Set((loadcfg, istp) for loadcfg in (:DCS, :SCS) for istp in val_istp)
Set(keys(curves_num)) == expected_curves ||
    throw(ArgumentError("$tag_head: expected DCS and SCS curves for $(collect(val_istp))"))

all_n_reps = vcat((curve.n_reps for curve in values(curves_num))...)
reps_min, reps_max = extrema(all_n_reps)
reps_used = reps_min == reps_max ? string(reps_min) : "$(reps_min)–$(reps_max)"
title_plot = dualmot_title(tag_head, reps_used; bias=0.0)

fig_nums = isnothing(plot_cmpr_loadcfg_target) ?
    Figure(size=plot_cmpr_loadcfg.size, fontsize=plot_cmpr_loadcfg.fontsize,
        figure_padding=1) : plot_cmpr_loadcfg_target.fig
slot_nums = isnothing(plot_cmpr_loadcfg_target) ? fig_nums[1, 1] :
    plot_cmpr_loadcfg_target.slot
frame_options_nums = isnothing(plot_cmpr_loadcfg_target) ?
    (; width=plot_cmpr_loadcfg.frame_size[1],
        height=plot_cmpr_loadcfg.frame_size[2]) :
    get(plot_cmpr_loadcfg_target, :frame_options, (; aspect=AxisAspect(4 / 3)))
axis_options_nums = dualmot_axis_kwargs(;
    compact_spacing=!isnothing(plot_cmpr_loadcfg_target) &&
        get(plot_cmpr_loadcfg_target, :compact_spacing, true))
!isnothing(plot_cmpr_loadcfg_target) &&
    (axis_options_nums = merge(axis_options_nums, plot_cmpr_loadcfg_target.axis_options))
ax_nums = Axis(slot_nums; xlabel=plot_cmpr_loadcfg.xlabel,
    ylabel=rich(plot_cmpr_loadcfg.ylabel_num, " (10", superscript("7"), ")"),
    title=isnothing(plot_cmpr_loadcfg_target) ? title_plot : "",
    frame_options_nums...,
    axis_options_nums...)
curves_nums_plot = [begin
    curve = curves_num[(loadcfg, istp)]
    style = dualmot_curve_style((; loadcfg, istp))
    nums_plot = ifelse.(mask_t_load_target,
        curve.nums ./ plot_cmpr_loadcfg.scale_num, NaN)
    mask_error = mask_t_load_target .& isfinite.(curve.nums) .& isfinite.(curve.stds)
    (; loadcfg, istp, curve, style, nums_plot, mask_error)
end for loadcfg in (:SCS, :DCS) for istp in val_istp]
for curve_plot in curves_nums_plot
    lines!(ax_nums, val_t_load_plot, curve_plot.nums_plot;
        curve_plot.style.line_options...)
end
for curve_plot in curves_nums_plot
    mask_error = curve_plot.mask_error
    mask_marker = isfinite.(curve_plot.nums_plot)
    label = "$(curve_plot.istp) $(curve_plot.loadcfg)"
    if any(mask_error)
        marker_errorbars!(ax_nums,
            val_t_load_plot[mask_error], curve_plot.nums_plot[mask_error],
            curve_plot.curve.stds[mask_error] ./ plot_cmpr_loadcfg.scale_num;
            curve_plot.style.marker_options...,
            curve_plot.style.errorbar_options...,
            marker=curve_plot.style.marker,
            label)
        mask_marker_only = mask_marker .& .!mask_error
        any(mask_marker_only) && scatter!(ax_nums,
            val_t_load_plot[mask_marker_only], curve_plot.nums_plot[mask_marker_only];
            curve_plot.style.marker_options...,
            marker=curve_plot.style.marker,
            label=nothing)
    else
        scatter!(ax_nums, val_t_load_plot[mask_marker], curve_plot.nums_plot[mask_marker];
            curve_plot.style.marker_options...,
            marker=curve_plot.style.marker,
            label)
    end
end
if isnothing(plot_cmpr_loadcfg_target)
    set_time_minor_ticks!(ax_nums)
    !isnothing(limits_linear_cmpr) && begin
        xlims!(ax_nums, limits_linear_cmpr.x...)
        ylims!(ax_nums, limits_linear_cmpr.y...)
    end
else
    if get(plot_cmpr_loadcfg_target, :show_legend, true)
        for loadcfg in (:DIS, :DDM), istp in val_istp
            style = dualmot_curve_style((; loadcfg, istp))
            scatter!(ax_nums, Float64[], Float64[];
                style.marker_options..., marker=style.marker,
                label="$(istp) $loadcfg")
        end
    end
    xlimits = plot_cmpr_loadcfg_target.limits.x
    ylimits = plot_cmpr_loadcfg_target.limits.y
    isnothing(xlimits) || xlims!(ax_nums, xlimits...)
    isnothing(ylimits) || ylims!(ax_nums, ylimits...)
end
if isnothing(plot_cmpr_loadcfg_target) ||
    get(plot_cmpr_loadcfg_target, :show_legend, true)
    axislegend(ax_nums; position=:rb, DUALMOT_LEGEND_OPTIONS...)
end

curves_ratio = Dict{Symbol, NamedTuple}()
if isnothing(plot_cmpr_loadcfg_target)
fig_ratio = Figure(size=plot_cmpr_loadcfg.size, fontsize=plot_cmpr_loadcfg.fontsize,
    figure_padding=1)
ax_ratio = Axis(fig_ratio[1, 1]; xlabel=plot_cmpr_loadcfg.xlabel,
    ylabel="DCS / SCS number", title=title_plot,
    width=plot_cmpr_loadcfg.frame_size[1],
    height=plot_cmpr_loadcfg.frame_size[2],
    dualmot_axis_kwargs()...)
hlines!(ax_ratio, [0.9, 1.1]; color=RGBAf(0, 0, 0, 0.45),
    linestyle=:dash, linewidth=0.5)
curves_ratio_plot = NamedTuple[]
for istp in val_istp
    curve_dcs = curves_num[(:DCS, istp)]
    curve_scs = curves_num[(:SCS, istp)]
    denom = curve_scs.nums
    mask_ratio = isfinite.(curve_dcs.nums) .& isfinite.(denom) .& (denom .> 0)
    ratios = fill(NaN, length(val_t_load))
    ratios[mask_ratio] .= curve_dcs.nums[mask_ratio] ./ denom[mask_ratio]
    stds = fill(NaN, length(val_t_load))
    mask_error = mask_ratio .& isfinite.(curve_dcs.stds) .& isfinite.(curve_scs.stds)
    stds[mask_error] .= sqrt.((curve_dcs.stds[mask_error] ./ denom[mask_error]).^2 .+
        (curve_dcs.nums[mask_error] .* curve_scs.stds[mask_error] ./ denom[mask_error].^2).^2)
    curves_ratio[istp] = (; ratios, stds)

    style = merge(dualmot_ratio_style(istp), (; marker=:diamond))
    push!(curves_ratio_plot, (; istp, ratios, stds, mask_error, style))
end
for curve_plot in curves_ratio_plot
    lines!(ax_ratio, val_t_load_plot, curve_plot.ratios;
        curve_plot.style.line_options...)
end
for curve_plot in curves_ratio_plot
    mask_error = curve_plot.mask_error
    mask_marker = isfinite.(curve_plot.ratios)
    label = string(curve_plot.istp)
    if any(mask_error)
        marker_errorbars!(ax_ratio,
            val_t_load_plot[mask_error], curve_plot.ratios[mask_error],
            curve_plot.stds[mask_error];
            curve_plot.style.marker_options...,
            curve_plot.style.errorbar_options...,
            marker=curve_plot.style.marker,
            label)
        mask_marker_only = mask_marker .& .!mask_error
        any(mask_marker_only) && scatter!(ax_ratio,
            val_t_load_plot[mask_marker_only], curve_plot.ratios[mask_marker_only];
            curve_plot.style.marker_options...,
            marker=curve_plot.style.marker,
            label=nothing)
    else
        scatter!(ax_ratio, val_t_load_plot[mask_marker], curve_plot.ratios[mask_marker];
            curve_plot.style.marker_options...,
            marker=curve_plot.style.marker,
            label)
    end
end
set_time_minor_ticks!(ax_ratio)
axislegend(ax_ratio; position=:rb, DUALMOT_LEGEND_OPTIONS...)

if isdefined(@__MODULE__, :figure_data_sheets)
    number_headers = Any["t_load (raw)", "effective t_load (s)"]
    for loadcfg in (:SCS, :DCS), istp in val_istp
        append!(number_headers, ["mean [$loadcfg $istp]", "std [$loadcfg $istp]",
            "n_rep [$loadcfg $istp]", "source [$loadcfg $istp]"])
    end
    number_matrix = Matrix{Any}(undef, length(val_t_load) + 1, length(number_headers))
    number_matrix[1, :] .= number_headers
    for idx_t in eachindex(val_t_load)
        number_matrix[idx_t + 1, 1:2] .= (val_t_load[idx_t], val_t_load_plot[idx_t])
        col = 3
        for loadcfg in (:SCS, :DCS), istp in val_istp
            curve = curves_num[(loadcfg, istp)]
            number_matrix[idx_t + 1, col:col + 3] .= (
                curve.nums[idx_t], curve.stds[idx_t],
                curve.n_reps[idx_t] == 0 ? NaN : curve.n_reps[idx_t],
                figure_xlsx_source(runinfo, (; t_load=val_t_load[idx_t], loadcfg, istp)))
            col += 4
        end
    end

    ratio_headers = Any["t_load (raw)", "effective t_load (s)"]
    for istp in val_istp
        append!(ratio_headers, ["ratio DCS/SCS [$istp]", "std [$istp]",
            "n_rep DCS [$istp]", "n_rep SCS [$istp]",
            "source DCS [$istp]", "source SCS [$istp]"])
    end
    ratio_matrix = Matrix{Any}(undef, length(val_t_load) + 1, length(ratio_headers))
    ratio_matrix[1, :] .= ratio_headers
    for idx_t in eachindex(val_t_load)
        ratio_matrix[idx_t + 1, 1:2] .= (val_t_load[idx_t], val_t_load_plot[idx_t])
        col = 3
        for istp in val_istp
            dcs, scs = curves_num[(:DCS, istp)], curves_num[(:SCS, istp)]
            ratio = curves_ratio[istp]
            ratio_matrix[idx_t + 1, col:col + 5] .= (
                ratio.ratios[idx_t], ratio.stds[idx_t],
                isfinite(ratio.ratios[idx_t]) && dcs.n_reps[idx_t] > 0 ? dcs.n_reps[idx_t] : NaN,
                isfinite(ratio.ratios[idx_t]) && scs.n_reps[idx_t] > 0 ? scs.n_reps[idx_t] : NaN,
                figure_xlsx_source(runinfo, (; t_load=val_t_load[idx_t],
                    loadcfg=:DCS, istp)),
                figure_xlsx_source(runinfo, (; t_load=val_t_load[idx_t],
                    loadcfg=:SCS, istp)))
            col += 6
        end
    end
    sheets = get!(figure_data_sheets, runinfo.folder, Pair{String,Matrix{Any}}[])
    push!(sheets, "$(runinfo.tag) numbers" => number_matrix)
    push!(sheets, "$(runinfo.tag) ratio" => ratio_matrix)
end

for format in plot_cmpr_loadcfg.formats
    save_options = format == "png" ? (; px_per_unit=4) : (;)
    save(joinpath(path_output, "[$(plot_cmpr_loadcfg.file_head)].[$(runinfo.tag)].[nums].$format"), fig_nums; save_options...)
    save(joinpath(path_output, "[$(plot_cmpr_loadcfg.file_head)].[$(runinfo.tag)].[ratio].$format"), fig_ratio; save_options...)
end
end
