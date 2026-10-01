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
include(joinpath(@__DIR__, "anlz_pair_comparison_plot_helpers.jl"))
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
    combined = combine_num_evol_blocks(stats_data;
        label="$(runinfo.folder) $(runinfo.tag)")
    vars = combined.vars
    idx_rep_axis = findfirst(==(:rep), keys(vars))
    num_stat = dropdims(mapslices(combined.num_fmt; dims=idx_rep_axis) do values
        valid = collect(skipmissing(vec(values)))
        isempty(valid) ? NaN : mean(valid)
    end; dims=idx_rep_axis)
    std_num_stat = dropdims(mapslices(combined.num_fmt; dims=idx_rep_axis) do values
        valid = collect(skipmissing(vec(values)))
        length(valid) < 2 ? NaN : std(valid)
    end; dims=idx_rep_axis)
    n_rep_stat = dropdims(sum(.!ismissing.(combined.num_fmt); dims=idx_rep_axis);
        dims=idx_rep_axis)
    name_stat = Tuple(key for key in keys(vars) if key != :rep)
    stats_combined = (; vars, name_stat, num_stat, std_num_stat, n_rep_stat)
    val_t_load = vars.t_load

    curves = Dict{Tuple{Symbol, Symbol}, NamedTuple}()
    for stats in (stats_combined,)
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
    runinfo::NamedTuple, bias::Union{Nothing,Real}=nothing,
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
            source=figure_xlsx_source(runinfo, merge((; loadcfg, istp,
                t_load=loading_curves.val_t_load[idx_final]),
                isnothing(bias) ? NamedTuple() : (; β_MOT=bias))),
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
    runinfo_421_comparison = select_for_load421_runinfo(runinfos_421[pair], pair)
    runinfo_626 = latest_num_evol_runinfo(runinfos_626[pair];
        label="$pair MOT loading 626 comparison")
    local val_istp = Symbol.(split(pair, "-"))

    curves_626 = collect_loading_curves(runinfo_626;
        bounds_sigmax_num, bounds_sigmay_num, num_max_num)
    points_626[pair] = final_loading_points(
        curves_626, (:DCS, :SCS), val_istp, "$pair MOT loading 626";
        runinfo=runinfo_626, idx_final=only(indexin(t_load_query, curves_626.val_t_load)))
    for variant in (:t_balanced, :n_balanced)
        bias = balance_biases[variant][pair]
        candidates_421 = filter([runinfo_421_comparison]) do runinfo
            any(data -> :β_MOT in keys(data.vars) && bias in data.vars.β_MOT, runinfo.data)
        end
        if isempty(candidates_421) && variant == :n_balanced
            @warn "$pair: no MOT loading 421 data at n-balanced β_MOT = $bias; leaving this comparison point blank"
            continue
        end
        runinfo_421 = latest_num_evol_runinfo(candidates_421;
            label="$pair MOT loading 421 $(variant) comparison")
        curves_421 = collect_loading_curves(runinfo_421;
            bias, bounds_sigmax_num, bounds_sigmay_num, num_max_num)
        points_421[variant][pair] = final_loading_points(
            curves_421, (:DDM, :DIS), val_istp, "$pair MOT loading 421 $(variant)";
            runinfo=runinfo_421, bias,
            idx_final=only(indexin(t_load_query, curves_421.val_t_load)))
        if variant == :t_balanced
            candidates_dcs = filter([runinfo_421_comparison]) do runinfo
                any(data -> :β_MOT in keys(data.vars) && bias in data.vars.β_MOT &&
                    :DCS in data.vars.loadcfg && t_load_query in data.vars.t_load,
                    runinfo.data)
            end
            if isempty(candidates_dcs)
                @warn "$pair: no t-balanced MOT loading 421 DCS data at $(t_load_query) sec; omitting DCS points"
            else
                runinfo_dcs = latest_num_evol_runinfo(candidates_dcs;
                    label="$pair t-balanced MOT loading 421 DCS comparison")
                curves_dcs = collect_loading_curves(runinfo_dcs;
                    bias, bounds_sigmax_num, bounds_sigmay_num, num_max_num)
                points_dcs = final_loading_points(curves_dcs, (:DCS,), val_istp,
                    "$pair t-balanced MOT loading 421 DCS";
                    runinfo=runinfo_dcs, bias,
                    idx_final=only(indexin(t_load_query, curves_dcs.val_t_load)))
                merge!(points_421[variant][pair], points_dcs)
            end
        end
    end
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
        haskey(points_by_pair[pair], (numerator, istp)) || continue
        haskey(points_by_pair[pair], (denominator, istp)) || continue
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
            haskey(points_by_pair[pair], key) || continue
            point = points_by_pair[pair][key]
            isfinite(point.num) && point.num > 0 ||
                throw(ArgumentError("$pair $key: number must be finite and positive"))
            x = pair_istp_position(idx_pair, pair, istp)
            style = dualmot_curve_style((; loadcfg, istp))
            marker_options = marker_style(style; markersize=6)
            if isfinite(point.std) && point.std >= 0 && point.num > limits_y[1]
                marker_errorbars_log!(ax, [x], [point.num], [point.std];
                    floor=limits_y[1], marker_options..., errorlinewidth=0.75)
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
        ratio_421_dis_dcs=variant == :t_balanced ? draw_pair_ratio(points, :DIS, :DCS;
            ylabel=rich("N", subscript("DIS"), " / N", subscript("DCS")),
            title="MOT loading 421 · 30 sec · β_MOT = 0 · t-balanced",
            filename="[MOT.loading.pairs].[DIS-DCS].[ratio.t-balanced]") : nothing,
        ratio_626=draw_pair_ratio(points_626, :DCS, :SCS;
            ylabel=rich("N", subscript("DCS"), " / N", subscript("SCS")),
            title="MOT loading 626 · 30 sec · $label",
            filename="[MOT.loading.pairs].[DCS-SCS].[ratio.$label]"),
        numbers_421=draw_pair_numbers(points,
            variant == :t_balanced ? (:DIS, :DDM, :DCS) : (:DIS, :DDM);
            limits_y=limits_y_numbers,
            ylabel=variant == :t_balanced ?
                rich("N", subscript("DDM"), ", N", subscript("DIS"), ", N", subscript("DCS")) :
                rich("N", subscript("DDM"), ", N", subscript("DIS")),
            title="MOT loading 421 · 30 sec · $bias_text · $label",
            filename=variant == :t_balanced ?
                "[MOT.loading.pairs].[DDM-DIS-DCS].[nums.$label]" :
                "[MOT.loading.pairs].[DDM-DIS].[nums.$label]"),
        numbers_626=draw_pair_numbers(points_626, (:SCS, :DCS);
            limits_y=limits_y_numbers,
            ylabel=rich("N", subscript("DCS"), ", N", subscript("SCS")),
            title="MOT loading 626 · 30 sec · $label",
            filename="[MOT.loading.pairs].[DCS-SCS].[nums.$label]"),
    )
end

function loading_pair_table(points_by_pair::AbstractDict, loadcfgs::Tuple;
    biases::AbstractDict, ratio=false, numerator=nothing, denominator=nothing,
    quantity::AbstractString)
    if ratio
        headers = Any["pair", "istp", "loadcfg ratio", "quantity", "tbiasmot",
            "t_load (s)", "ratio", "std", "n_rep numerator", "n_rep denominator",
            "source numerator", "source denominator"]
    else
        headers = Any["pair", "istp", "loadcfg", "quantity", "tbiasmot",
            "t_load (s)", "mean", "std", "n_rep", "source", "source denominator"]
    end
    rows = Vector{Vector{Any}}()
    for pair in val_pair, istp in Symbol.(split(pair, "-"))
        loadcfg_iter = ratio ? (nothing,) : loadcfgs
        for loadcfg in loadcfg_iter
            if ratio
                point_num = get(get(points_by_pair, pair, Dict()), (numerator, istp), nothing)
                point_den = get(get(points_by_pair, pair, Dict()), (denominator, istp), nothing)
                valid = !isnothing(point_num) && !isnothing(point_den) &&
                    isfinite(point_num.num) && isfinite(point_den.num) && point_den.num > 0
                value = valid ? point_num.num / point_den.num : NaN
                err = valid && isfinite(point_num.std) && isfinite(point_den.std) ?
                    sqrt((point_num.std / point_den.num)^2 +
                        (point_num.num * point_den.std / point_den.num^2)^2) : NaN
                push!(rows, Any[pair, string(istp), "$(numerator)/$(denominator)",
                    quantity, get(biases, pair, NaN), ismissing(point_num) ? NaN :
                    isnothing(point_num) ? NaN : point_num.t_load, value, err,
                    isnothing(point_num) || !isfinite(point_num.num) || point_num.n_rep == 0 ? NaN : point_num.n_rep,
                    isnothing(point_den) || !isfinite(point_den.num) || point_den.n_rep == 0 ? NaN : point_den.n_rep,
                    isnothing(point_num) ? "" : point_num.source,
                    isnothing(point_den) ? "" : point_den.source])
            else
                point = get(get(points_by_pair, pair, Dict()), (loadcfg, istp), nothing)
                push!(rows, Any[pair, string(istp), string(loadcfg), quantity,
                    get(biases, pair, NaN), isnothing(point) ? t_load_query : point.t_load,
                    isnothing(point) ? NaN : point.num,
                    isnothing(point) ? NaN : point.std,
                    isnothing(point) || point.n_rep == 0 ? NaN : point.n_rep,
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

comparison_sheets_421 = Pair{String,Matrix{Any}}[]
for variant in (:t_balanced, :n_balanced)
    label = variant == :t_balanced ? "t" : "n"
    points = points_421[variant]
    push!(comparison_sheets_421, "$label numbers" => loading_pair_table(points,
        variant == :t_balanced ? (:DDM, :DIS, :DCS) : (:DDM, :DIS);
        biases=balance_biases[variant], quantity="MOT loading 421"))
    push!(comparison_sheets_421, "$label DDM-DIS ratio" => loading_pair_table(points,
        (); biases=balance_biases[variant], ratio=true, numerator=:DDM,
        denominator=:DIS, quantity="MOT loading 421 ratio"))
    if variant == :t_balanced
        push!(comparison_sheets_421, "t DIS-DCS ratio" => loading_pair_table(points,
            (); biases=balance_biases[variant], ratio=true, numerator=:DIS,
            denominator=:DCS, quantity="MOT loading 421 ratio"))
    end
end
write_figure_workbook(joinpath(path_output, "MOT loading 421.xlsx"), comparison_sheets_421)
comparison_sheets_626 = Pair{String,Matrix{Any}}[
    "t numbers" => loading_pair_table(points_626, (:DCS, :SCS);
        biases=Dict(pair => NaN for pair in val_pair), quantity="MOT loading 626"),
    "t DCS-SCS ratio" => loading_pair_table(points_626, ();
        biases=Dict(pair => NaN for pair in val_pair), ratio=true,
        numerator=:DCS, denominator=:SCS, quantity="MOT loading 626 ratio"),
]
write_figure_workbook(joinpath(path_output, "MOT loading 626.xlsx"), comparison_sheets_626)
