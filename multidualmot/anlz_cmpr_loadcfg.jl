length(runinfo.data) == 3 ||
    throw(ArgumentError("$tag_head: expected exactly three rectangular data blocks"))
stats_data = [calc_num_evol_block(data;
    label="$tag_head data[$idx]", bounds_sigmax_num, bounds_sigmay_num, num_max_num)
    for (idx, data) in enumerate(runinfo.data)]

val_istp = Symbol.(split(runinfo.folder, "-"))
val_t_load = stats_data[1].vars.t_load
all(stats -> stats.vars.t_load == val_t_load, stats_data) ||
    throw(ArgumentError("$tag_head: all data blocks must use the same t_load values"))

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

fig_nums = Figure(size=plot_cmpr_loadcfg.size)
ax_nums = Axis(fig_nums[1, 1]; xlabel=plot_cmpr_loadcfg.xlabel,
    ylabel="$(plot_cmpr_loadcfg.ylabel_num) (×10⁷)", title=title_plot,
    yminorticks=IntervalsBetween(5), yminorticksvisible=true)
for loadcfg in (:DCS, :SCS), istp in val_istp
    curve = curves_num[(loadcfg, istp)]
    style = dualmot_curve_style((; loadcfg, istp))
    nums_plot = curve.nums ./ plot_cmpr_loadcfg.scale_num
    scatterlines!(ax_nums, val_t_load, nums_plot; style..., label="$istp $loadcfg")
    mask_error = isfinite.(curve.nums) .& isfinite.(curve.stds)
    errorbars!(ax_nums, val_t_load[mask_error], nums_plot[mask_error],
        curve.stds[mask_error] ./ plot_cmpr_loadcfg.scale_num;
        color=style.color, whiskerwidth=7, linewidth=1.2)
end
dualmot_axislegend(ax_nums; position=:rb)

curves_ratio = Dict{Symbol, NamedTuple}()
fig_ratio = Figure(size=plot_cmpr_loadcfg.size)
ax_ratio = Axis(fig_ratio[1, 1]; xlabel=plot_cmpr_loadcfg.xlabel,
    ylabel="DCS / SCS number", title=title_plot,
    yminorticks=IntervalsBetween(5), yminorticksvisible=true)
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

    style = dualmot_ratio_style(istp)
    scatterlines!(ax_ratio, val_t_load, ratios; style..., label=string(istp))
    errorbars!(ax_ratio, val_t_load[mask_error], ratios[mask_error], stds[mask_error];
        color=style.color, whiskerwidth=7, linewidth=1.2)
end
dualmot_axislegend(ax_ratio; position=:rb)

for format in plot_cmpr_loadcfg.formats
    save(joinpath(path_output, "[$(plot_cmpr_loadcfg.file_head)].[$(runinfo.tag)].[nums].$format"), fig_nums)
    save(joinpath(path_output, "[$(plot_cmpr_loadcfg.file_head)].[$(runinfo.tag)].[ratio].$format"), fig_ratio)
end
