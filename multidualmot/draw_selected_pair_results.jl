"""Draw the isotope-pair comparison figure for the selected Results folder."""

path_result_pairs = joinpath(path_root, "Result")
mkpath(path_result_pairs)
axis_label_fontsize = 7 / 0.75
fig_selected_pair_comparison = Figure(size=(640, 495), fontsize=8,
    figure_padding=3, rowgap=4, colgap=6)

function result_pair_axis(fig::Figure, slot; ylabel, yticks=nothing,
    log_y::Bool=false, ratio::Bool=false)
    kwargs = dualmot_axis_kwargs(; log_y, text_size=8)
    isnothing(yticks) || (kwargs = merge(kwargs, (; yticks)))
    ratio && (kwargs = merge(kwargs, (; yminorticks=IntervalsBetween(2),
        yminorticksvisible=true)))
    kwargs = merge(kwargs, (; xlabelsize=axis_label_fontsize,
        ylabelsize=axis_label_fontsize))
    Axis(fig[slot...]; xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair", ylabel,
        yscale=log_y ? log10 : identity, kwargs...)
end

# MOT loading 421, 30-second t-balanced DDM/DIS atom numbers.
ax_result_421_numbers = result_pair_axis(fig_selected_pair_comparison,
    (1, 1); ylabel="CMOT N", log_y=true)
draw_pair_spans!(ax_result_421_numbers)
for (idx_pair, pair) in enumerate(val_pair),
    loadcfg in (:DDM, ), istp in Symbol.(split(pair, "-"))
    haskey(points_421[:t_balanced], pair) || continue
    point = get(points_421[:t_balanced][pair], (loadcfg, istp), nothing)
    isnothing(point) && continue
    isfinite(point.num) && point.num > 0 || continue
    style = dualmot_curve_style((; loadcfg, istp))
    options = marker_style(style; markersize=6)
    x = pair_istp_position(idx_pair, pair, istp)
    if isfinite(point.std) && point.std >= 0 && point.num > 1e6
        marker_errorbars_log!(ax_result_421_numbers, [x], [point.num], [point.std];
            floor=1e6, options..., errorlinewidth=0.75)
    else
        scatter!(ax_result_421_numbers, [x], [point.num]; options...)
    end
end
ylims!(ax_result_421_numbers, 0.8e6, 1.2e8)

# MOT loading 421 DDM/DIS number ratio, t-balanced.
ax_result_421_ratio = result_pair_axis(fig_selected_pair_comparison,
    (1, 2); ylabel="CMOT N DDM—DIS ratio", yticks=0:0.2:1.2, ratio=true)
draw_pair_spans!(ax_result_421_ratio)
hlines!(ax_result_421_ratio, [0.9, 1.1]; color=RGBAf(0, 0, 0, 0.45),
    linestyle=:dash, linewidth=0.5)
for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-"))
    haskey(points_421[:t_balanced], pair) || continue
    pair_points = points_421[:t_balanced][pair]
    haskey(pair_points, (:DDM, istp)) && haskey(pair_points, (:DIS, istp)) || continue
    top = pair_points[(:DDM, istp)]
    bottom = pair_points[(:DIS, istp)]
    isfinite(top.num) && isfinite(bottom.num) && bottom.num > 0 || continue
    ratio = top.num / bottom.num
    error = isfinite(top.std) && isfinite(bottom.std) ?
        sqrt((top.std / bottom.num)^2 + (top.num * bottom.std / bottom.num^2)^2) : NaN
    style = dualmot_curve_style((; loadcfg=:DDM, istp))
    options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
    x = pair_istp_position(idx_pair, pair, istp)
    if isfinite(error)
        marker_errorbars!(ax_result_421_ratio, [x], [ratio], [error];
            options..., errorlinewidth=0.75)
    else
        scatter!(ax_result_421_ratio, [x], [ratio]; options...)
    end
end
ylims!(ax_result_421_ratio, 0.35, 1.25)

# CMOT decay inverse-kappa ratio, with references at 1.0
ax_result_cmot_ratio = result_pair_axis(fig_selected_pair_comparison,
    (2, 1); ylabel="κ DDM—DIS ratio", yticks=0:0.2:1.2,
    ratio=true)
draw_pair_spans!(ax_result_cmot_ratio)
hlines!(ax_result_cmot_ratio, [1.0]; color=RGBAf(0, 0, 0, 0.45),
    linestyle=:dash, linewidth=0.5)
points_cmot_ratio = comparison_points_lifetime[:n_balanced].cmot
for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-"))
    haskey(points_cmot_ratio, pair) || continue
    pair_points = points_cmot_ratio[pair]
    haskey(pair_points, (:DDM, istp)) && haskey(pair_points, (:DIS, istp)) || continue
    top = pair_points[(:DDM, istp)]
    bottom = pair_points[(:DIS, istp)]
    isfinite(top.value) && top.value > 0 && isfinite(bottom.value) && bottom.value > 0 || continue
    ratio = top.value / bottom.value
    error = isfinite(top.std) && isfinite(bottom.std) ?
        sqrt((top.std / bottom.value)^2 + (top.value * bottom.std / bottom.value^2)^2) : NaN
    style = dualmot_curve_style((; loadcfg=:DDM, istp))
    options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
    x = pair_istp_position(idx_pair, pair, istp)
    if isfinite(error)
        marker_errorbars!(ax_result_cmot_ratio, [x], [ratio], [error];
            options..., errorlinewidth=0.75)
    else
        scatter!(ax_result_cmot_ratio, [x], [ratio]; options...)
    end
end
ylims!(ax_result_cmot_ratio, 0.35, 1.25)

# DDM/DIS ratio of ODT-to-CMOT efficiencies; x/z directions are offset slightly.
ax_result_odt_efficiency = result_pair_axis(fig_selected_pair_comparison,
    (2, 2); ylabel="ODT loading efficiency \n DDM—DIS ratio", yticks=0:0.2:1.2,
    ratio=true)
draw_pair_spans!(ax_result_odt_efficiency)
hlines!(ax_result_odt_efficiency, [0.9, 1.1]; color=RGBAf(0, 0, 0, 0.45),
    linestyle=:dash, linewidth=0.5)
for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-"))
    direction_results = filter(result -> !isnothing(result), [
        odt_cmot_efficiency_ratio(direction, pair, istp) for direction in (:x, :z)])
    isempty(direction_results) && continue
    result = direction_results[argmax(getproperty.(direction_results, :value))]
    style = dualmot_curve_style((; loadcfg=:DDM, istp))
    options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
    x = pair_istp_position(idx_pair, pair, istp)
    if isfinite(result.error)
        marker_errorbars!(ax_result_odt_efficiency, [x], [result.value], [result.error];
            options..., errorlinewidth=0.75)
    else
        scatter!(ax_result_odt_efficiency, [x], [result.value]; options...)
    end
end
ylims!(ax_result_odt_efficiency, 0.35, 1.25)

# MOT loading 626 DCS/SCS CMOT-number ratio, at 30 seconds.
ax_result_626_ratio = result_pair_axis(fig_selected_pair_comparison,
    (3, 1); ylabel="CMOT N DCS—SCS ratio", yticks=0:0.2:1.2,
    ratio=true)
draw_pair_spans!(ax_result_626_ratio)
hlines!(ax_result_626_ratio, [0.9, 1.1]; color=RGBAf(0, 0, 0, 0.45),
    linestyle=:dash, linewidth=0.5)
for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-"))
    haskey(points_626, pair) || continue
    pair_points = points_626[pair]
    haskey(pair_points, (:DCS, istp)) && haskey(pair_points, (:SCS, istp)) || continue
    top = pair_points[(:DCS, istp)]
    bottom = pair_points[(:SCS, istp)]
    isfinite(top.num) && top.num >= 0 && isfinite(bottom.num) && bottom.num > 0 || continue
    ratio = top.num / bottom.num
    error = isfinite(top.std) && isfinite(bottom.std) ?
        sqrt((top.std / bottom.num)^2 + (top.num * bottom.std / bottom.num^2)^2) : NaN
    style = dualmot_curve_style((; loadcfg=:DCS, istp))
    options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
    x = pair_istp_position(idx_pair, pair, istp)
    if isfinite(error)
        marker_errorbars!(ax_result_626_ratio, [x], [ratio], [error];
            options..., errorlinewidth=0.75)
    else
        scatter!(ax_result_626_ratio, [x], [ratio]; options...)
    end
end
ylims!(ax_result_626_ratio, 0, 1.25)

# MOT lifetime τ DDM/DIS ratio, using the n-balanced fit results.
ax_result_mot_lifetime_ratio = result_pair_axis(fig_selected_pair_comparison,
    (3, 2); ylabel=rich("τ", subscript("DDM"), " / τ", subscript("DIS")),
    yticks=0:0.2:1.2, ratio=true)
draw_pair_spans!(ax_result_mot_lifetime_ratio)
hlines!(ax_result_mot_lifetime_ratio, [0.9, 1.1];
    color=RGBAf(0, 0, 0, 0.45), linestyle=:dash, linewidth=0.5)
points_mot_lifetime_ratio = comparison_points_lifetime[:n_balanced].mot
for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-"))
    haskey(points_mot_lifetime_ratio, pair) || continue
    pair_points = points_mot_lifetime_ratio[pair]
    haskey(pair_points, (:DDM, istp)) && haskey(pair_points, (:DIS, istp)) || continue
    top = pair_points[(:DDM, istp)]
    bottom = pair_points[(:DIS, istp)]
    isfinite(top.value) && top.value > 0 && isfinite(bottom.value) && bottom.value > 0 || continue
    ratio = top.value / bottom.value
    error = isfinite(top.std) && isfinite(bottom.std) ?
        sqrt((top.std / bottom.value)^2 + (top.value * bottom.std / bottom.value^2)^2) : NaN
    style = dualmot_curve_style((; loadcfg=:DDM, istp))
    options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
    x = pair_istp_position(idx_pair, pair, istp)
    if isfinite(error)
        marker_errorbars!(ax_result_mot_lifetime_ratio, [x], [ratio], [error];
            options..., errorlinewidth=0.75)
    else
        scatter!(ax_result_mot_lifetime_ratio, [x], [ratio]; options...)
    end
end
ylims!(ax_result_mot_lifetime_ratio, 0, 1.25)

for format in ("svg", "png", "pdf")
    save_options = format == "png" ? (; px_per_unit=4) : (;)
    save(joinpath(path_result_pairs, "selected_isotope_pair_comparison.$format"),
        fig_selected_pair_comparison; save_options...)
end
println("Saved selected isotope-pair comparison figure to $path_result_pairs")
