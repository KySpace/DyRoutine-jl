"""Draw the isotope-pair comparison figure for the selected Results folder."""

path_result_pairs = joinpath(path_root, "Result")
mkpath(path_result_pairs)
axis_ticklabel_fontsize = 7 / 0.75
axis_label_fontsize = 8 / 0.75
result_axis_font = "Times New Roman"
result_pair_axis_width = round(Int, 75 / 25.4 * 96)
result_pair_axis_height = round(Int, result_pair_axis_width / (320 / 149))
result_pair_tick_labels = replace.(val_pair, "-" => "–")
result_number_label() = rich(result_math_italic("N"),
    subscript(result_math_italic("i")), superscript("CMOT"))
result_ratio_label(quantity::AbstractString, qualifier::AbstractString) =
    rich(result_math_italic("R"),
        subscript(result_math_italic("i," * quantity)), superscript(qualifier))
fig_selected_pair_comparison = Figure(size=(fig_width, 322), fontsize=8 / 0.75,
    figure_padding=1, rowgap=4, colgap=2)

function result_pair_axis(fig::Figure, slot; ylabel, yticks=nothing,
    log_y::Bool=false, ratio::Bool=false)
    kwargs = dualmot_axis_kwargs(; log_y, text_size=axis_ticklabel_fontsize)
    isnothing(yticks) || (kwargs = merge(kwargs, (; yticks)))
    ratio && (kwargs = merge(kwargs, (; yminorticks=IntervalsBetween(2),
        yminorticksvisible=true)))
    kwargs = merge(kwargs, (; xlabelsize=axis_label_fontsize,
        ylabelsize=axis_label_fontsize,
        xticklabelsize=axis_ticklabel_fontsize,
        yticklabelsize=axis_ticklabel_fontsize,
        xlabelfont=result_axis_font,
        ylabelfont=result_axis_font,
        xticklabelfont=result_axis_font,
        yticklabelfont=result_axis_font))
    Axis(fig[slot...]; xticks=(eachindex(val_pair), result_pair_tick_labels),
        xlabel="Isotope pair", ylabel,
        yscale=log_y ? log10 : identity,
        width=result_pair_axis_width, height=result_pair_axis_height,
        aspect=AxisAspect(320 / 149), kwargs...)
end

# MOT loading 421, 30-second t-balanced DDM/DIS atom numbers.
ax_result_421_numbers = result_pair_axis(fig_selected_pair_comparison,
    (1, 1); ylabel=result_number_label(), log_y=true)
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
    (1, 2); ylabel=result_ratio_label("N", "load"),
    yticks=0:0.2:1.2, ratio=true)
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
    (2, 1); ylabel=result_ratio_label("κ", "CMOT"), yticks=0:0.2:1.2,
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
    (2, 2); ylabel=result_ratio_label("η", "ODT"), yticks=0:0.2:1.2,
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

for format in ("svg", "png", "pdf")
    save_options = format == "png" ? (; px_per_unit=4) : (;)
    save(joinpath(path_result_pairs, "selected_isotope_pair_comparison.$format"),
        fig_selected_pair_comparison; save_options...)
end
println("Saved selected isotope-pair comparison figure to $path_result_pairs")
