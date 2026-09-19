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
title_plot = "$tag_head · reps = $reps_used"

fig_nums = isnothing(plot_cmpr_loadcfg_target) ?
    Figure(size=plot_cmpr_loadcfg.size, fontsize=plot_cmpr_loadcfg.fontsize,
        figure_padding=1) : plot_cmpr_loadcfg_target.fig
slot_nums = isnothing(plot_cmpr_loadcfg_target) ? fig_nums[1, 1] :
    plot_cmpr_loadcfg_target.slot
frame_options_nums = isnothing(plot_cmpr_loadcfg_target) ?
    (; width=plot_cmpr_loadcfg.frame_size[1],
        height=plot_cmpr_loadcfg.frame_size[2]) :
    (; aspect=AxisAspect(4 / 3))
ax_nums = Axis(slot_nums; xlabel=plot_cmpr_loadcfg.xlabel,
    ylabel="$(plot_cmpr_loadcfg.ylabel_num) (×10⁷)",
    title=isnothing(plot_cmpr_loadcfg_target) ? title_plot : "",
    frame_options_nums...,
    dualmot_axis_kwargs(; compact_spacing=!isnothing(plot_cmpr_loadcfg_target))...)
curves_nums_plot = [begin
    curve = curves_num[(loadcfg, istp)]
    style = dualmot_curve_style((; loadcfg, istp))
    nums_plot = curve.nums ./ plot_cmpr_loadcfg.scale_num
    mask_error = isfinite.(curve.nums) .& isfinite.(curve.stds)
    (; loadcfg, istp, curve, style, nums_plot, mask_error)
end for loadcfg in (:DCS, :SCS) for istp in val_istp]
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
set_time_minor_ticks!(ax_nums)
axislegend(ax_nums; position=:rb, DUALMOT_LEGEND_OPTIONS...)

curves_ratio = Dict{Symbol, NamedTuple}()
if isnothing(plot_cmpr_loadcfg_target)
fig_ratio = Figure(size=plot_cmpr_loadcfg.size, fontsize=plot_cmpr_loadcfg.fontsize,
    figure_padding=1)
ax_ratio = Axis(fig_ratio[1, 1]; xlabel=plot_cmpr_loadcfg.xlabel,
    ylabel="DCS / SCS number", title=title_plot,
    width=plot_cmpr_loadcfg.frame_size[1],
    height=plot_cmpr_loadcfg.frame_size[2],
    dualmot_axis_kwargs()...)
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

for format in plot_cmpr_loadcfg.formats
    save(joinpath(path_output, "[$(plot_cmpr_loadcfg.file_head)].[$(runinfo.tag)].[nums].$format"), fig_nums)
    save(joinpath(path_output, "[$(plot_cmpr_loadcfg.file_head)].[$(runinfo.tag)].[ratio].$format"), fig_ratio)
end
end
