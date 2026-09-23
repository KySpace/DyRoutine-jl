include(joinpath(@__DIR__, "dualmotcommons.jl"))
using CSV

path_root_421 = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading 421"
path_root_626 = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading 626"
path_output = joinpath(dirname(path_root_421), "Isotope pair comparison")
val_pair = ["160-162", "162-164", "161-162", "161-164", "162-163",
    "163-164", "161-163"]
pair_istp_offset = 0.05 # isotope marker displacement as a fraction of the pair interval
isfinite(pair_istp_offset) && 0 <= pair_istp_offset < 0.5 ||
    throw(ArgumentError("pair_istp_offset must be finite and in [0, 0.5), got $pair_istp_offset"))

function pair_istp_position(idx_pair::Integer, pair::AbstractString, istp::Symbol)
    isotopes = Symbol.(split(pair, "-"))
    length(isotopes) == 2 || throw(ArgumentError("invalid isotope pair: $pair"))
    direction = istp == isotopes[1] ? -1 : istp == isotopes[2] ? 1 :
        throw(ArgumentError("$istp is not part of isotope pair $pair"))
    idx_pair + direction * pair_istp_offset
end
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8
formats_output = ("svg", "png")

function collect_loading_curves(runinfo::NamedTuple;
    bias::Union{Nothing,Real}=nothing,
    bounds_sigmax_num::Tuple{<:Real,<:Real},
    bounds_sigmay_num::Tuple{<:Real,<:Real},
    num_max_num::Real,
)
    stats_data = [calc_num_evol_block(data;
        label="$(runinfo.folder) $(runinfo.tag) data[$idx]",
        bounds_sigmax_num, bounds_sigmay_num, num_max_num)
        for (idx, data) in enumerate(runinfo.data)]
    val_t_load = only(unique([stats.vars.t_load for stats in stats_data]))
    issorted(val_t_load) ||
        throw(ArgumentError("$(runinfo.folder): t_load values must be sorted"))

    curves = Dict{Tuple{Symbol, Symbol}, NamedTuple}()
    for stats in stats_data
        idx_axis = Dict(key => idx for (idx, key) in enumerate(stats.name_stat))
        idx_bias = if :β_MOT in stats.name_stat
            isnothing(bias) && throw(ArgumentError("$(runinfo.folder): β_MOT selection is required"))
            idx = findfirst(==(bias), stats.vars.β_MOT)
            isnothing(idx) &&
                throw(ArgumentError("$(runinfo.folder): β_MOT = $bias is unavailable"))
            idx
        else
            isnothing(bias) || throw(ArgumentError("$(runinfo.folder): unexpected β_MOT selection"))
            nothing
        end
        for (idx_loadcfg, loadcfg) in enumerate(stats.vars.loadcfg),
            (idx_istp, istp) in enumerate(stats.vars.istp)
            key = (loadcfg, istp)
            haskey(curves, key) &&
                throw(ArgumentError("$(runinfo.folder): duplicate curve $key"))
            indices = Any[Colon() for _ in stats.name_stat]
            isnothing(idx_bias) || (indices[idx_axis[:β_MOT]] = idx_bias)
            indices[idx_axis[:loadcfg]] = idx_loadcfg
            indices[idx_axis[:istp]] = idx_istp
            curves[key] = (
                nums=vec(@view stats.num_stat[indices...]),
                stds=vec(@view stats.std_num_stat[indices...]),
                n_reps=vec(@view stats.n_rep_stat[indices...]),
            )
        end
    end
    (; val_t_load, curves)
end

function final_loading_points(loading_curves::NamedTuple,
    loadcfgs::Tuple, val_istp::AbstractVector{Symbol}, label::AbstractString;
    idx_final = lastindex(loading_curves.val_t_load)
    )
    points = Dict{Tuple{Symbol, Symbol}, NamedTuple}()
    for loadcfg in loadcfgs, istp in val_istp
        key = (loadcfg, istp)
        haskey(loading_curves.curves, key) ||
            throw(ArgumentError("$label: missing curve $key"))
        curve = loading_curves.curves[key]
        points[key] = (
            num=curve.nums[idx_final],
            std=curve.stds[idx_final],
            n_rep=curve.n_reps[idx_final],
            t_load=loading_curves.val_t_load[idx_final],
        )
    end
    points
end

runinfos_421 = Dict(pair => read_num_evol_runinfos(path_root_421, pair;
    var_specs=DUALMOT_LOADING_VAR_SPECS,
    validate_vars=validate_dualmot_vars,
) for pair in val_pair)
runinfos_626 = Dict(pair => read_num_evol_runinfos(path_root_626, pair;
    var_specs=DUALMOT_LOADCFG_VAR_SPECS,
    validate_vars=validate_loadcfg_comparison_vars,
) for pair in val_pair)

points_421 = Dict{Symbol, Dict{String, Dict{Tuple{Symbol, Symbol}, NamedTuple}}}(
    variant => Dict{String, Dict{Tuple{Symbol, Symbol}, NamedTuple}}()
    for variant in (:t_balanced, :n_balanced))
points_626 = Dict{String, Dict{Tuple{Symbol, Symbol}, NamedTuple}}()
t_load_query = 30.0 # 30 seconds loading time
rows_balance = collect(CSV.File(joinpath(dirname(path_root_421), "MOT loading balance", "balance.csv");
    header=["Pair", "tbiasmot"], skipto=2, stripwhitespace=true))
canonical_comparison_pair_name(pair) = join(sort(parse.(Int, split(strip(string(pair)), "-"))), "-")
balance_by_pair = Dict(canonical_comparison_pair_name(row.Pair) => Float64(row.tbiasmot) for row in rows_balance)
Set(keys(balance_by_pair)) == Set(val_pair) ||
    throw(ArgumentError("MOT loading balance/balance.csv must contain exactly the configured isotope pairs"))
balance_biases = Dict(
    :t_balanced => Dict(pair => 0.0 for pair in val_pair),
    :n_balanced => balance_by_pair,
)
for pair in val_pair
    runinfo_626 = only(runinfos_626[pair])
    val_istp = Symbol.(split(pair, "-"))

    curves_626 = collect_loading_curves(runinfo_626;
        bounds_sigmax_num, bounds_sigmay_num, num_max_num)
    points_626[pair] = final_loading_points(
        curves_626, (:DCS, :SCS), val_istp, "$pair MOT loading 626"; idx_final=only(indexin(t_load_query, curves_626.val_t_load)))
    for variant in (:t_balanced, :n_balanced)
        bias = balance_biases[variant][pair]
        candidates_421 = filter(runinfos_421[pair]) do runinfo
            any(data -> :β_MOT in keys(data.vars) && bias in data.vars.β_MOT, runinfo.data)
        end
        if isempty(candidates_421) && variant == :n_balanced
            @warn "$pair: no MOT loading 421 data at n-balanced β_MOT = $bias; leaving this comparison point blank"
            continue
        end
        length(candidates_421) == 1 || throw(ArgumentError(
            "$pair $(variant): expected one 421 processing entry containing β_MOT = $bias"))
        runinfo_421 = only(candidates_421)
        curves_421 = collect_loading_curves(runinfo_421;
            bias, bounds_sigmax_num, bounds_sigmay_num, num_max_num)
        points_421[variant][pair] = final_loading_points(
            curves_421, (:DDM, :DIS), val_istp, "$pair MOT loading 421 $(variant)";
            idx_final=only(indexin(t_load_query, curves_421.val_t_load)))
    end
end

marker_style(style::NamedTuple; markersize::Real=style.markersize) = (
    color=style.markercolor,
    marker=style.marker,
    markersize,
    strokecolor=style.strokecolor,
    strokewidth=style.strokewidth,
)

function draw_pair_spans!(ax::Axis)
    for idx_pair in eachindex(val_pair)
        vspan!(ax, idx_pair - 0.5, idx_pair + 0.5; color=:white) |>
            span -> translate!(span, 0, 0, -100)
    end
    nothing
end

function style_key_axis(fig::Figure; width::Real, limits)
    border_color = RGBAf(0.25, 0.25, 0.25, 0.45)
    Axis(fig[1, 1];
        width,
        height=58,
        halign=0.05,
        valign=0.10,
        tellwidth=false,
        tellheight=false,
        limits,
        xticksvisible=false,
        xticklabelsvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
        xgridvisible=false,
        ygridvisible=false,
        bottomspinecolor=border_color,
        leftspinecolor=border_color,
        topspinecolor=border_color,
        rightspinecolor=border_color,
        spinewidth=0.75,
        backgroundcolor=RGBAf(1, 1, 1, 0.86),
    )
end

function draw_number_style_key!(fig::Figure, loadcfgs::Tuple)
    val_istp_key = Symbol.(string.(160:164))
    pos_y = collect(5:-1:1)
    pos_x = 1.10 .+ 0.34 .* (0:length(loadcfgs)-1)
    ax_key = style_key_axis(fig;
        width=56,
        limits=(0.65, 1.70, 0.5, 6.2),
    )
    text!(ax_key, mean(pos_x), 5.8;
        text=join(string.(loadcfgs), "  "),
        align=(:center, :center), fontsize=22 / 3)
    for (idx_istp, istp) in enumerate(val_istp_key)
        text!(ax_key, 0.72, pos_y[idx_istp]; text=string(istp),
            align=(:left, :center), fontsize=22 / 3)
    end
    for (idx_loadcfg, loadcfg) in enumerate(loadcfgs),
        (idx_istp, istp) in enumerate(val_istp_key)
        style = dualmot_curve_style((; loadcfg, istp))
        scatter!(ax_key, [pos_x[idx_loadcfg]], [pos_y[idx_istp]];
            marker_style(style; markersize=6)...)
    end
    ax_key
end

function draw_ratio_style_key!(fig::Figure)
    val_istp_key = Symbol.(string.(160:164))
    pos_y = collect(5:-1:1)
    ax_key = style_key_axis(fig;
        width=28,
        limits=(0.65, 1.25, 0.5, 5.5),
    )
    for (idx_istp, istp) in enumerate(val_istp_key)
        text!(ax_key, 0.7, pos_y[idx_istp]; text=string(istp),
            align=(:left, :center), fontsize=22 / 3)
        style = dualmot_curve_style((; loadcfg=:DDM, istp))
        scatter!(ax_key, [1.07], [pos_y[idx_istp]];
            merge(marker_style(style; markersize=6), (; marker=:hexagon))...)
    end
    ax_key
end

function draw_pair_ratio(points_by_pair::AbstractDict, numerator::Symbol, denominator::Symbol;
    ylabel, title, filename::AbstractString)
    fig = Figure(size=(320, 149), fontsize=8, figure_padding=1)
    ax = Axis(fig[1, 1];
        xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair",
        yticks=(0:0.2:1.2),
        yminorticks=IntervalsBetween(2),
        ylabel,
        title,
        dualmot_axis_kwargs(; text_size=8)...,
    )
    draw_pair_spans!(ax)
    hlines!(ax, [0.9, 1.1]; color=RGBAf(0, 0, 0, 0.45),
        linestyle=:dash, linewidth=0.5)
    for (idx_pair, pair) in enumerate(val_pair),
        istp in Symbol.(split(pair, "-"))
        haskey(points_by_pair, pair) || continue
        point_num = points_by_pair[pair][(numerator, istp)]
        point_den = points_by_pair[pair][(denominator, istp)]
        isfinite(point_num.num) && isfinite(point_den.num) && point_den.num > 0 ||
            throw(ArgumentError("$pair $istp: invalid $numerator/$denominator ratio"))
        ratio = point_num.num / point_den.num
        std_ratio = if isfinite(point_num.std) && isfinite(point_den.std)
            sqrt((point_num.std / point_den.num)^2 +
                (point_num.num * point_den.std / point_den.num^2)^2)
        else
            NaN
        end
        x = pair_istp_position(idx_pair, pair, istp)
        style = dualmot_curve_style((; loadcfg=numerator, istp))
        marker_options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
        if isfinite(std_ratio)
            marker_errorbars!(ax, [x], [ratio], [std_ratio];
                marker_options..., errorlinewidth=0.75)
        else
            scatter!(ax, [x], [ratio]; marker_options...)
        end
    end
    ylims!(ax, 0, 1.25)
    ax.yminorticks = IntervalsBetween(2)
    draw_ratio_style_key!(fig)
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
    fig
end

function draw_pair_numbers(points_by_pair::AbstractDict, loadcfgs::Tuple;
    ylabel, title, filename::AbstractString, limits_y::Tuple)
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
    for (idx_pair, pair) in enumerate(val_pair)
        haskey(points_by_pair, pair) || continue
        val_istp = Symbol.(split(pair, "-"))
        for loadcfg in loadcfgs, istp in val_istp
            key = (loadcfg, istp)
            point = points_by_pair[pair][key]
            isfinite(point.num) && point.num > 0 ||
                throw(ArgumentError("$pair $key: number must be finite and positive"))
            x = pair_istp_position(idx_pair, pair, istp)
            style = dualmot_curve_style((; loadcfg, istp))
            marker_options = marker_style(style; markersize=6)
            if isfinite(point.std) && point.num - point.std > 0
                marker_errorbars!(ax, [x], [point.num], [point.std];
                    marker_options..., errorlinewidth=0.75)
            else
                scatter!(ax, [x], [point.num]; marker_options...)
            end
        end
    end
    ylims!(ax, limits_y...)
    draw_number_style_key!(fig, loadcfgs)
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
    fig
end

mkpath(path_output)
limits_y_numbers = (0.3e6, 1.6e8)
figs_loading_pairs = Dict{Symbol,NamedTuple}()
for variant in (:t_balanced, :n_balanced)
    label = variant == :t_balanced ? "t-balanced" : "n-balanced"
    points = points_421[variant]
    bias_text = variant == :t_balanced ? "β_MOT = 0" : "pair balance β_MOT"
    figs_loading_pairs[variant] = (
        ratio_421=draw_pair_ratio(points, :DDM, :DIS;
            ylabel=rich("N", subscript("DDM"), " / N", subscript("DIS")),
            title="MOT loading 421 · 30 sec · $bias_text · $label",
            filename="[MOT.loading.pairs].[DDM-DIS].[ratio.$label]"),
        ratio_626=draw_pair_ratio(points_626, :DCS, :SCS;
            ylabel=rich("N", subscript("DCS"), " / N", subscript("SCS")),
            title="MOT loading 626 · 30 sec · $label",
            filename="[MOT.loading.pairs].[DCS-SCS].[ratio.$label]"),
        numbers_421=draw_pair_numbers(points, (:DIS, :DDM);
            limits_y=limits_y_numbers,
            ylabel=rich("N", subscript("DDM"), ", N", subscript("DIS")),
            title="MOT loading 421 · 30 sec · $bias_text · $label",
            filename="[MOT.loading.pairs].[DDM-DIS].[nums.$label]"),
        numbers_626=draw_pair_numbers(points_626, (:SCS, :DCS);
            limits_y=limits_y_numbers,
            ylabel=rich("N", subscript("DCS"), ", N", subscript("SCS")),
            title="MOT loading 626 · 30 sec · $label",
            filename="[MOT.loading.pairs].[DCS-SCS].[nums.$label]"),
    )
end

include(joinpath(@__DIR__, "anlz_lifetime_pair_comparison.jl"))
