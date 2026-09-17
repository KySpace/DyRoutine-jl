function calc_num_evol_block(runinfo_data::NamedTuple;
    label::AbstractString,
    size_min_num::Real,
    num_max_num::Real,
)
    shot_data = read_cres_fields(runinfo_data.files, ("atomnum", "sigmax", "sigmay"))
    num_data, sigmax_data, sigmay_data = shot_data
    all(isfinite, num_data) || throw(ArgumentError("$label: atomnum must contain only finite values"))
    len_data = length(num_data)
    n_variation = prod(length(getproperty(runinfo_data.vars, key))
        for key in keys(runinfo_data.vars) if key != :rep)
    len_data > 0 && rem(len_data, n_variation) == 0 ||
        throw(DimensionMismatch("$label: $len_data atom numbers must be a positive integer multiple of $n_variation variations"))
    n_rep = div(len_data, n_variation)
    vars = merge(runinfo_data.vars, (; rep=1:n_rep))
    name = keys(vars)

    name_acq = runinfo_data.var_order
    n_dims_acq = map(key -> length(getproperty(vars, key)), name_acq)
    mask_size = isfinite.(sigmax_data) .& isfinite.(sigmay_data) .&
        (sigmax_data .>= size_min_num) .& (sigmay_data .>= size_min_num)
    mask_num_low = num_data .>= 0
    mask_num_high = num_data .<= num_max_num
    mask_valid = mask_size .& mask_num_low .& mask_num_high
    num_data_masked = Vector{Union{Missing, Float64}}(num_data)
    num_data_masked[.!mask_valid] .= missing

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

    n_masked_size = count(!, mask_size)
    n_masked_num_low = count(!, mask_num_low)
    n_masked_num_high = count(!, mask_num_high)
    n_masked_total = count(!, mask_valid)
    println("$label: $len_data samples / $n_variation variations = $n_rep repetitions; " *
        "$n_masked_total rejected ($n_masked_size by size, $n_masked_num_low below zero, " *
        "$n_masked_num_high above $num_max_num)")
    (; vars, name_stat, num_fmt, num_stat, std_num_stat, n_rep_stat, n_rep,
        mask_valid, n_masked_total)
end

length(runinfo.data) == 3 ||
    throw(ArgumentError("$tag_head: expected exactly three rectangular data blocks"))
stats_data = [calc_num_evol_block(data;
    label="$tag_head data[$idx]", size_min_num, num_max_num)
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
