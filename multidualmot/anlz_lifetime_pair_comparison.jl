if !isdefined(@__MODULE__, :pair_istp_offset)
    global pair_istp_offset = 0.05
end
isdefined(@__MODULE__, :pair_istp_position) ||
    include(joinpath(@__DIR__, "anlz_pair_comparison_plot_helpers.jl"))

function decay_parameter_points(latest::NamedTuple, expected_mode::Symbol,
    parameter::Symbol, biases::AbstractDict{String,<:Real})
    latest.fit_mode == expected_mode ||
        throw(ArgumentError("expected $expected_mode fits in $(latest.path), got $(latest.fit_mode)"))
    points = Dict{String, Dict{Tuple{Symbol,Symbol},NamedTuple}}()
    for pair in val_pair
        records_pair = filter(record -> record.pair == pair, latest.records)
        isempty(records_pair) && continue
        panel = Float64(biases[pair])
        panel in unique(record.panel for record in records_pair) || continue
        records_panel = filter(record -> record.panel == panel, records_pair)
        val_istp = Symbol.(split(pair, "-"))
        expected = Set((loadcfg, istp) for loadcfg in (:DIS, :DDM) for istp in val_istp)
        actual = Set((record.loadcfg, record.istp) for record in records_panel)
        issubset(actual, expected) && length(records_panel) == length(actual) ||
            throw(ArgumentError("$(latest.path): $pair panel $panel has unexpected or duplicate conditions $actual"))
        points[pair] = Dict((record.loadcfg, record.istp) => (
            value=Float64(getproperty(record, parameter)),
            std=Float64(getproperty(record, Symbol("std_", parameter))),
            n0=Float64(record.n0),
            std_n0=Float64(record.std_n0),
            panel,
            n_rep=Float64(get(record, :n_points, NaN)),
            source=String(get(record, :sources, basename(latest.path))),
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
            panel=point.panel,
            n_rep=point.n_rep,
            source=point.source,
            n0=point.n0,
            std_n0=point.std_n0)
    end for (key, point) in points) for (pair, points) in points_by_pair)
end

function draw_pair_lifetime_values(points_by_pair::AbstractDict, value_key::Symbol,
    error_key::Symbol; ylabel, title::AbstractString, filename::AbstractString,
    limits_y::Union{Nothing,Tuple}=nothing)
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
        haskey(points_by_pair, pair) || continue
        haskey(points_by_pair[pair], (loadcfg, istp)) || continue
        point = points_by_pair[pair][(loadcfg, istp)]
        value = Float64(getproperty(point, value_key))
        std = Float64(getproperty(point, error_key))
        isfinite(value) && value > 0 ||
            throw(ArgumentError("$pair $loadcfg $istp: $value_key must be finite and positive"))
        x = pair_istp_position(idx_pair, pair, istp)
        style = dualmot_curve_style((; loadcfg, istp))
        marker_options = marker_style(style; markersize=6)
        if isfinite(std) && std >= 0
            if value - std > 0
                marker_errorbars!(ax, [x], [value], [std];
                    marker_options..., errorlinewidth=0.75)
                lower = min(lower, value - std)
                upper = max(upper, value + std)
            else
                scatter!(ax, [x], [value]; marker_options...)
                lower = min(lower, value)
                upper = max(upper, value)
            end
        else
            scatter!(ax, [x], [value]; marker_options...)
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
    if isnothing(limits_y)
        ylims!(ax, lower_limit, upper_limit)
    else
        ylims!(ax, limits_y...)
    end
    draw_number_style_key!(fig, (:DIS, :DDM))
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
    fig
end

draw_pair_decay_values(points_by_pair::AbstractDict, parameter::Symbol;
    kwargs...) = draw_pair_lifetime_values(points_by_pair, :value, :std; kwargs...)

draw_pair_n0_values(points_by_pair::AbstractDict;
    title::AbstractString, filename::AbstractString) =
    draw_pair_lifetime_values(points_by_pair, :n0, :std_n0;
        ylabel="N₀ (atom)", title, filename)

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
        haskey(points_by_pair, pair) || continue
        haskey(points_by_pair[pair], (numerator, istp)) || continue
        haskey(points_by_pair[pair], (denominator, istp)) || continue
        point_num = points_by_pair[pair][(numerator, istp)]
        point_den = points_by_pair[pair][(denominator, istp)]
        isfinite(point_num.value) && point_num.value > 0 &&
            isfinite(point_den.value) && point_den.value > 0 ||
            throw(ArgumentError("$pair $istp: invalid $numerator/$denominator fit ratio"))
        ratio = point_num.value / point_den.value
        x = pair_istp_position(idx_pair, pair, istp)
        std_ratio = if isfinite(point_num.std) && isfinite(point_den.std)
            sqrt((point_num.std / point_den.value)^2 +
                (point_num.value * point_den.std / point_den.value^2)^2)
        else
            NaN
        end
        style = dualmot_curve_style((; loadcfg=numerator, istp))
        marker_options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
        if isfinite(std_ratio)
            marker_errorbars!(ax, [x], [ratio], [std_ratio];
                marker_options..., errorlinewidth=0.75)
            upper = max(upper, ratio + std_ratio)
        else
            scatter!(ax, [x], [ratio]; marker_options...)
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

figs_lifetime_pairs = Dict{Symbol,NamedTuple}()
comparison_points_lifetime = Dict{Symbol,NamedTuple}()
for (variant, biases) in ((:t_balanced, balance_biases[:t_balanced]),
    (:n_balanced, balance_biases[:n_balanced]))
    label = variant == :t_balanced ? "t-balanced" : "n-balanced"
    points_cmot_kappa = decay_parameter_points(latest_cmot_decay, :kappa, :kappa, biases)
    points_mot_tau = decay_parameter_points(latest_mot_decay, :tau, :tau, biases)
    points_cmot_inverse_kappa = inverse_decay_points(points_cmot_kappa)
    comparison_points_lifetime[variant] = (; cmot=points_cmot_inverse_kappa,
        mot=points_mot_tau)
    figs_lifetime_pairs[variant] = (
        cmot_values=draw_pair_decay_values(points_cmot_inverse_kappa, :inverse_kappa;
            ylabel=rich("1 / κ (atom s)"),
            title="CMOT decay · κ-only fit · $label",
            filename="[CMOT.decay.pairs].[kappa].[values.$label]",
            limits_y=(0.5e6, 1.0e7)),
        cmot_n0=draw_pair_n0_values(points_cmot_inverse_kappa;
            title="CMOT decay · fitted N₀ · $label",
            filename="[CMOT.decay.pairs].[n0].[values.$label]"),
        cmot_ratio=draw_pair_decay_ratio(points_cmot_inverse_kappa, :DDM, :DIS;
            ylabel=rich("(1 / κ)", subscript("DDM"), " / (1 / κ)", subscript("DIS")),
            title="CMOT decay · κ-only fit · $label",
            filename="[CMOT.decay.pairs].[DIS-DDM].[ratio.$label]"),
        mot_values=draw_pair_decay_values(points_mot_tau, :tau;
            ylabel="τ (s)",
            title="MOT decay · τ-only fit · $label",
            filename="[MOT.decay.pairs].[tau].[values.$label]",
            limits_y=(0.9e0, 5.5e1)),
        mot_n0=draw_pair_n0_values(points_mot_tau;
            title="MOT decay · fitted N₀ · $label",
            filename="[MOT.decay.pairs].[n0].[values.$label]"),
        mot_ratio=draw_pair_decay_ratio(points_mot_tau, :DDM, :DIS;
            ylabel=rich("τ", subscript("DDM"), " / τ", subscript("DIS")),
            title="MOT decay · τ-only fit · $label",
            filename="[MOT.decay.pairs].[DDM-DIS].[ratio.$label]"),
    )
end

function lifetime_pair_matrix(points_by_pair::AbstractDict;
    biases::AbstractDict, parameter::AbstractString, ratio=false, n0=false,
    numerator=:DDM, denominator=:DIS, quantity::AbstractString)
    headers = ratio ? Any["pair", "istp", "loadcfg ratio", "CMOT/MOT value",
        "tbiasmot", "ratio", "std", "n_rep numerator", "n_rep denominator",
        "source numerator", "source denominator"] :
        Any["pair", "istp", "loadcfg", "CMOT/MOT value", "tbiasmot",
            n0 ? "N₀" : "value", "std", "n_rep", "source", "source denominator"]
    rows = Vector{Vector{Any}}()
    for pair in val_pair, istp in Symbol.(split(pair, "-"))
        if ratio
            point_num = get(get(points_by_pair, pair, Dict()), (numerator, istp), nothing)
            point_den = get(get(points_by_pair, pair, Dict()), (denominator, istp), nothing)
            valid = !isnothing(point_num) && !isnothing(point_den) &&
                isfinite(point_num.value) && point_num.value > 0 &&
                isfinite(point_den.value) && point_den.value > 0
            value = valid ? point_num.value / point_den.value : NaN
            err = valid && isfinite(point_num.std) && isfinite(point_den.std) ?
                sqrt((point_num.std / point_den.value)^2 +
                    (point_num.value * point_den.std / point_den.value^2)^2) : NaN
            push!(rows, Any[pair, string(istp), "$(numerator)/$(denominator)",
                quantity, get(biases, pair, NaN), value, err,
                isnothing(point_num) ? NaN : point_num.n_rep,
                isnothing(point_den) ? NaN : point_den.n_rep,
                isnothing(point_num) ? "" : point_num.source,
                isnothing(point_den) ? "" : point_den.source])
        else
            for loadcfg in (:DDM, :DIS)
                point = get(get(points_by_pair, pair, Dict()), (loadcfg, istp), nothing)
                value = isnothing(point) ? NaN : getproperty(point, n0 ? :n0 : :value)
                error = isnothing(point) ? NaN : getproperty(point, n0 ? :std_n0 : :std)
                push!(rows, Any[pair, string(istp), string(loadcfg), quantity,
                    isnothing(point) ? get(biases, pair, NaN) : point.panel,
                    value,
                    error,
                    isnothing(point) ? NaN : point.n_rep,
                    isnothing(point) ? "" : point.source, ""])
            end
        end
    end
    matrix = Matrix{Any}(undef, length(rows) + 1, length(headers))
    matrix[1, :] .= headers
    for (idx, row) in enumerate(rows)
        matrix[idx + 1, :] .= row
    end
    matrix
end

for (dataset, parameter, field) in (("CMOT lifetime", "1/κ (s)", :cmot),
    ("MOT lifetime", "τ (s)", :mot))
    local sheets = Pair{String,Matrix{Any}}[]
    for variant in (:t_balanced, :n_balanced)
        label = variant == :t_balanced ? "t" : "n"
        points = comparison_points_lifetime[variant][field]
        biases = balance_biases[variant]
        push!(sheets, "$label values" => lifetime_pair_matrix(points;
            biases, parameter, quantity=parameter))
        push!(sheets, "$label N₀" => lifetime_pair_matrix(points;
            biases, parameter, n0=true, quantity="N₀ (atom)"))
        push!(sheets, "$label DDM-DIS ratio" => lifetime_pair_matrix(points;
            biases, parameter, ratio=true, numerator=:DDM, denominator=:DIS,
            quantity="$dataset $parameter ratio"))
    end
    write_figure_workbook(joinpath(path_output, "$dataset.xlsx"), sheets)
end
