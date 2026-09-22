function decay_parameter_points(latest::NamedTuple, expected_mode::Symbol,
    parameter::Symbol)
    latest.fit_mode == expected_mode ||
        throw(ArgumentError("expected $expected_mode fits in $(latest.path), got $(latest.fit_mode)"))
    points = Dict{String, Dict{Tuple{Symbol,Symbol},NamedTuple}}()
    for pair in val_pair
        records_pair = filter(record -> record.pair == pair, latest.records)
        isempty(records_pair) && throw(ArgumentError("$(latest.path): no records for $pair"))
        panels = unique(record.panel for record in records_pair)
        panel = 0.0 in panels ? 0.0 : only(panels)
        records_panel = filter(record -> record.panel == panel, records_pair)
        val_istp = Symbol.(split(pair, "-"))
        expected = Set((loadcfg, istp) for loadcfg in (:DIS, :DDM) for istp in val_istp)
        actual = Set((record.loadcfg, record.istp) for record in records_panel)
        length(records_panel) == length(expected) && actual == expected || throw(ArgumentError(
            "$(latest.path): $pair panel $panel has conditions $actual, expected $expected"))
        points[pair] = Dict((record.loadcfg, record.istp) => (
            value=Float64(getproperty(record, parameter)),
            std=Float64(getproperty(record, Symbol("std_", parameter))),
            panel,
        ) for record in records_panel)
    end
    points
end

function inverse_decay_points(points_by_pair::AbstractDict)
    Dict(pair => Dict(key => begin
        isfinite(point.value) && point.value > 0 ||
            throw(ArgumentError("$pair $key: κ must be finite and positive"))
        (value=inv(point.value),
            std=isfinite(point.std) ? point.std / point.value^2 : NaN,
            panel=point.panel)
    end for (key, point) in points) for (pair, points) in points_by_pair)
end

function draw_pair_decay_values(points_by_pair::AbstractDict, parameter::Symbol;
    ylabel, title::AbstractString, filename::AbstractString, scale_y::Real=1.0, limits_y::Tuple)
    fig = Figure(size=(320, 164), fontsize=8, figure_padding=1)
    ax = Axis(fig[1, 1];
        xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair",
        ylabel,
        title,
        yscale=log10,
        dualmot_axis_kwargs(; log_y=true, text_size=8)...,
    )
    draw_pair_spans!(ax)
    lower = Inf
    upper = 0.0
    for (idx_pair, pair) in enumerate(val_pair), loadcfg in (:DIS, :DDM),
        istp in Symbol.(split(pair, "-"))
        point = points_by_pair[pair][(loadcfg, istp)]
        isfinite(point.value) && point.value > 0 ||
            throw(ArgumentError("$pair $loadcfg $istp: $parameter must be finite and positive"))
        value = point.value / scale_y
        std = point.std / scale_y
        style = dualmot_curve_style((; loadcfg, istp))
        marker_options = marker_style(style; markersize=6)
        if isfinite(std) && std >= 0
            value - std > 0 || throw(ArgumentError(
                "$pair $loadcfg $istp: log-scale $parameter error bar crosses zero"))
            marker_errorbars!(ax, [idx_pair], [value], [std];
                marker_options..., errorlinewidth=0.75)
            lower = min(lower, value - std)
            upper = max(upper, value + std)
        else
            scatter!(ax, [idx_pair], [value]; marker_options...)
            lower = min(lower, value)
            upper = max(upper, value)
        end
    end
    lower_limit, upper_limit = lower / 1.15, 1.15 * upper
    has_decade_tick = any(lower_limit <= 10.0^exponent <= upper_limit for exponent in -20:20)
    if !has_decade_tick
        lower_limit = 10.0^floor(log10(lower)) / 1.02
        upper_limit = 10.0^ceil(log10(upper)) * 1.02
    end
    ylims!(ax, limits_y...)
    draw_number_style_key!(fig, (:DIS, :DDM))
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
    fig
end

function draw_pair_decay_ratio(points_by_pair::AbstractDict,
    numerator::Symbol, denominator::Symbol;
    ylabel, title::AbstractString, filename::AbstractString)
    fig = Figure(size=(320, 149), fontsize=8, figure_padding=1)
    ax = Axis(fig[1, 1];
        xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair",
        yticks=0:0.2:1.2,
        yminorticks=IntervalsBetween(2),
        yminorticksvisible=true,
        ylabel,
        title,
        dualmot_axis_kwargs(; text_size=8)...,
    )
    draw_pair_spans!(ax)
    hlines!(ax, [0.9, 1.1]; color=RGBAf(0, 0, 0, 0.45),
        linestyle=:dash, linewidth=0.5)
    upper = 1.1
    for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-"))
        point_num = points_by_pair[pair][(numerator, istp)]
        point_den = points_by_pair[pair][(denominator, istp)]
        isfinite(point_num.value) && point_num.value > 0 &&
            isfinite(point_den.value) && point_den.value > 0 ||
            throw(ArgumentError("$pair $istp: invalid $numerator/$denominator fit ratio"))
        ratio = point_num.value / point_den.value
        std_ratio = if isfinite(point_num.std) && isfinite(point_den.std)
            sqrt((point_num.std / point_den.value)^2 +
                (point_num.value * point_den.std / point_den.value^2)^2)
        else
            NaN
        end
        style = dualmot_curve_style((; loadcfg=numerator, istp))
        marker_options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
        if isfinite(std_ratio)
            marker_errorbars!(ax, [idx_pair], [ratio], [std_ratio];
                marker_options..., errorlinewidth=0.75)
            upper = max(upper, ratio + std_ratio)
        else
            scatter!(ax, [idx_pair], [ratio]; marker_options...)
            upper = max(upper, ratio)
        end
    end
    ylims!(ax, 0, max(1.25, 1.08 * upper))
    ax.yminorticks = IntervalsBetween(2)
    draw_ratio_style_key!(fig)
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
    fig
end

latest_cmot_decay = load_latest_num_decay_results(
    joinpath(dirname(path_root_421), "CMOT lifetime"), "CMOT")
latest_mot_decay = load_latest_num_decay_results(
    joinpath(dirname(path_root_421), "MOT lifetime"), "MOT")
println("Using CMOT decay fits: $(latest_cmot_decay.path)")
println("Using MOT decay fits: $(latest_mot_decay.path)")

points_cmot_kappa = decay_parameter_points(latest_cmot_decay, :kappa, :kappa)
points_mot_tau = decay_parameter_points(latest_mot_decay, :tau, :tau)
points_cmot_inverse_kappa = inverse_decay_points(points_cmot_kappa)

fig_cmot_kappa = draw_pair_decay_values(points_cmot_inverse_kappa, :inverse_kappa;
    ylabel=rich("1 / κ (atom s)"),
    title="CMOT decay · κ-only fit",
    filename="[CMOT.decay.pairs].[kappa].[values]",
    limits_y = (0.5e6, 1.0e7),
)
fig_cmot_kappa_ratio = draw_pair_decay_ratio(points_cmot_inverse_kappa, :DDM, :DIS;
    ylabel=rich("(1 / κ)", subscript("DDM"), " / (1 / κ)", subscript("DIS")),
    title="CMOT decay · κ-only fit",
    filename="[CMOT.decay.pairs].[DIS-DDM].[ratio]",
)
fig_mot_tau = draw_pair_decay_values(points_mot_tau, :tau;
    ylabel="τ (s)",
    title="MOT decay · τ-only fit",
    filename="[MOT.decay.pairs].[tau].[values]",
    limits_y = (0.9e0, 5.5e1),
)
fig_mot_tau_ratio = draw_pair_decay_ratio(points_mot_tau, :DDM, :DIS;
    ylabel=rich("τ", subscript("DDM"), " / τ", subscript("DIS")),
    title="MOT decay · τ-only fit",
    filename="[MOT.decay.pairs].[DDM-DIS].[ratio]",
)
