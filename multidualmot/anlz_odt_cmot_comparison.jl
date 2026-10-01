include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root_odt = joinpath(path_root, "ODT BField")
path_root_421 = joinpath(path_root, "MOT loading 421")
path_output = joinpath(dirname(path_root_odt), "Isotope pair comparison")
val_pair = ["160-162", "162-164", "161-162", "161-164", "162-163", "163-164", "161-163"]
bounds_sigmax_num = (2e-4, 10e-4)
bounds_sigmay_num = (1e-4, 5e-4)
num_max_num = 2e8
mot_bounds_sigmax_num = (4e-4, Inf)
mot_bounds_sigmay_num = (4e-4, Inf)
formats_output = ("svg", "png")
mkpath(path_output)

function comparison_stats(runinfo; bounds_sigmax_num, bounds_sigmay_num, num_max_num)
    blocks = [calc_num_evol_block(data;
        label="$(runinfo.folder) $(runinfo.tag) data[$idx]",
        bounds_sigmax_num, bounds_sigmay_num, num_max_num)
        for (idx, data) in enumerate(runinfo.data)]
    combined = combine_num_evol_blocks(blocks; label="$(runinfo.folder) $(runinfo.tag)")
    num_fmt = combined.num_fmt
    idx_rep = findfirst(==(:rep), keys(combined.vars))
    num_stat = dropdims(mapslices(num_fmt; dims=idx_rep) do values
        valid = collect(skipmissing(vec(values)))
        isempty(valid) ? NaN : mean(valid)
    end; dims=idx_rep)
    std_stat = dropdims(mapslices(num_fmt; dims=idx_rep) do values
        valid = collect(skipmissing(vec(values)))
        length(valid) < 2 ? NaN : std(valid)
    end; dims=idx_rep)
    n_rep_stat = dropdims(sum(.!ismissing.(num_fmt); dims=idx_rep); dims=idx_rep)
    (; vars=combined.vars, name=Tuple(k for k in keys(combined.vars) if k != :rep),
        num_stat, std_stat, n_rep_stat, runinfo)
end

function stats_point(stats, condition::NamedTuple)
    idx = Any[Colon() for _ in stats.name]
    for key in keys(condition)
        axis_idx = findfirst(==(key), stats.name)
        isnothing(axis_idx) && throw(ArgumentError("missing axis $key in $(stats.name)"))
        idx[axis_idx] = findfirst(==(getproperty(condition, key)), getproperty(stats.vars, key))
        isnothing(idx[axis_idx]) && throw(ArgumentError("condition $key=$(getproperty(condition,key)) unavailable"))
    end
    (; num=only(stats.num_stat[idx...]), std=only(stats.std_stat[idx...]),
        n_rep=only(stats.n_rep_stat[idx...]),
        source=figure_xlsx_source(stats.runinfo, condition))
end

function draw_marker_point!(ax, x, y, err; marker_options, label=nothing,
    log_floor::Union{Nothing,Real}=nothing)
    if isfinite(err) && err >= 0
        if isnothing(log_floor)
            marker_errorbars!(ax, [x], [y], [err]; marker_options..., errorlinewidth=0.8, label)
        elseif y > log_floor
            marker_errorbars_log!(ax, [x], [y], [err]; floor=log_floor,
                marker_options..., errorlinewidth=0.8, label)
        end
    else
        scatter!(ax, [x], [y]; marker_options..., label)
    end
    nothing
end

odt_runinfos = Dict(pair => read_num_evol_runinfos(path_root_odt, pair;
    var_specs=DUALMOT_ODT_BFIELD_VAR_SPECS, validate_vars=validate_odt_bfield_vars) for pair in val_pair)
mot_runinfos = Dict(pair => read_num_evol_runinfos(path_root_421, pair;
    var_specs=DUALMOT_LOADING_VAR_SPECS, validate_vars=validate_dualmot_vars) for pair in val_pair)
source_lines = readlines(joinpath(path_root_odt, "source table.txt"))

points_odt = Dict{String, Dict{Tuple{Symbol,Symbol,Symbol},NamedTuple}}()
points_cmot = Dict{String, Dict{Tuple{Symbol,Symbol},NamedTuple}}()
conditions = Dict{String,NamedTuple}()
for pair in val_pair
    source_pair_lines = filter(line -> startswith(strip(line), pair * " "), source_lines)
    isempty(source_pair_lines) && throw(ArgumentError("$pair: no source table rows"))
    source_conditions = map(source_pair_lines) do line
        fields = split(strip(line))
        length(fields) >= 6 || throw(ArgumentError("$pair: malformed ODT source table row"))
        parse(Float64, fields[5]), begin
            matches = collect(eachmatch(r"\[([-+0-9.eE]+)\]", line))
            length(matches) >= 1 || throw(ArgumentError("$pair: source table row is missing tbiasmot"))
            parse(Float64, only(matches[1].captures))
        end
    end
    all(==(first(source_conditions)), source_conditions) ||
        throw(ArgumentError("$pair: ODT source rows disagree on tmotload/tbiasmot: $source_conditions"))
    tmotload, bias = first(source_conditions)
    conditions[pair] = (; bias, tmotload)
    points_odt[pair] = Dict{Tuple{Symbol,Symbol,Symbol},NamedTuple}()
    for runinfo in odt_runinfos[pair]
        stats = comparison_stats(runinfo; bounds_sigmax_num, bounds_sigmay_num, num_max_num)
        for direction in (Symbol(runinfo.tag),), loadcfg in (:DDM, :DIS), istp in Symbol.(split(pair, "-"))
            # Each processing entry is independent and has one direction; take its largest mean.
            candidates = [stats_point(stats, (; ib, loadcfg, istp)) for ib in stats.vars.ib]
            idx_max = findmax(getproperty.(candidates, :num))[2]
            key = (direction, loadcfg, istp)
            points_odt[pair][key] = candidates[idx_max]
        end
    end
    # Keep one 421 processing entry with this pair's bias, then select matching loading time.
    bias_candidates = filter(mot_runinfos[pair]) do runinfo
        any(data -> :β_MOT in keys(data.vars) && bias in data.vars.β_MOT, runinfo.data)
    end
    config_421 = YAML.load_file(joinpath(path_root_421, pair, "config.yaml"))
    selected_tag = if length(config_421) == 1
        String(only(config_421)["tag"])
    else
        tags_for_odt = [String(entry["tag"]) for entry in config_421
            if get(entry, "for_odt", false) === true]
        length(tags_for_odt) == 1 || throw(ArgumentError(
            "$pair: multiple 421 processing tags require exactly one for_odt: true, got $tags_for_odt"))
        only(tags_for_odt)
    end
    matching_tags = filter(runinfo -> runinfo.tag == selected_tag, mot_runinfos[pair])
    runinfo_mot = only(matching_tags)
    runinfo_mot in bias_candidates || throw(ArgumentError(
        "$pair: selected 421 tag '$selected_tag' has no data at tbiasmot=$bias"))
    mot_stats = comparison_stats(runinfo_mot;
        bounds_sigmax_num=mot_bounds_sigmax_num,
        bounds_sigmay_num=mot_bounds_sigmay_num, num_max_num)
    tmotload in mot_stats.vars.t_load || throw(ArgumentError("$pair: t_load=$tmotload unavailable"))
    points_cmot[pair] = Dict((loadcfg, istp) => stats_point(mot_stats,
        (; β_MOT=bias, t_load=tmotload, loadcfg, istp)) for loadcfg in (:DDM,:DIS), istp in Symbol.(split(pair,"-")))
end

# Plot two cross-pair figures. ODT x and z maxima are both shown when available.
function draw_comparison(; ratio::Bool)
    fig = Figure(size=(690, 340), fontsize=11, figure_padding=4)
    ax = Axis(fig[1,1]; xticks=(eachindex(val_pair), ["$pair\n$(conditions[pair].bias >= 0 ? "+" : "")$(round(conditions[pair].bias,digits=3)) $(Int(conditions[pair].tmotload))s" for pair in val_pair]),
        xlabel="Isotope pair · tbiasmot / tmotload", ylabel=ratio ? "ODT / CMOT number" : "Atom number",
        title=ratio ? "ODT / CMOT" : "ODT (filled) · CMOT (open)",
        yscale=ratio ? identity : log10, dualmot_axis_kwargs(; log_y=!ratio, text_size=9)...)
    xlims!(ax, 0.5, length(val_pair)+0.5)
    ratio && hlines!(ax, [1]; color=(:black,0.45), linestyle=:dash, linewidth=0.8)
    log_floor = nothing
    if !ratio
        log_values = Float64[]
        log_errors = Float64[]
        for pair in val_pair
            append!(log_values, (pt.num for pt in values(points_odt[pair])))
            append!(log_errors, (pt.std for pt in values(points_odt[pair])))
            append!(log_values, (pt.num for pt in values(points_cmot[pair])))
            append!(log_errors, (pt.std for pt in values(points_cmot[pair])))
        end
        log_limits = log_plot_limits(log_values, log_errors)
        ylims!(ax, log_limits...)
        log_floor = log_limits[1]
    end
    for (ip,pair) in enumerate(val_pair), loadcfg in (:DDM,:DIS), istp in Symbol.(split(pair,"-"))
        style = dualmot_curve_style((; loadcfg, istp))
        pos = ip + (istp == Symbol(first(split(pair,"-"))) ? -0.07 : 0.07) + (loadcfg == :DDM ? -0.025 : 0.025)
        mot = points_cmot[pair][(loadcfg,istp)]
        if ratio
            odt_pts = [v for ((_,lc,it),v) in points_odt[pair] if lc == loadcfg && it == istp]
            for pt in odt_pts
                value = pt.num / mot.num
                err = sqrt((pt.std/mot.num)^2 + (pt.num*mot.std/mot.num^2)^2)
                draw_marker_point!(ax, pos, value, err; marker_options=(;
                    color=style.markercolor, marker=style.marker, markersize=7,
                    strokecolor=style.strokecolor, strokewidth=style.strokewidth))
            end
        else
            for ((dir,lc,it),pt) in points_odt[pair]
                lc == loadcfg && it == istp || continue
                shift = dir == :x ? -0.025 : 0.025
                draw_marker_point!(ax, pos+shift, pt.num, pt.std; marker_options=(;
                    color=style.markercolor, marker=style.marker, markersize=7,
                    strokecolor=style.strokecolor, strokewidth=style.strokewidth),
                    log_floor)
            end
            hollow = (; color=:transparent, marker=style.marker, markersize=7,
                strokecolor=style.strokecolor, strokewidth=style.strokewidth)
            draw_marker_point!(ax, pos, mot.num, mot.std; marker_options=hollow,
                log_floor)
        end
    end
    filename = ratio ? "[ODT.CMOT.comparison].[lin].[ratios]" : "[ODT.CMOT.comparison].[log].[numbers]"
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=3) : (;)
        save(joinpath(path_output,"$filename.$format"), fig; save_options...)
    end
    fig
end

function odt_cmot_efficiency_ratio(direction::Symbol, pair::AbstractString,
    istp::Symbol)
    odt_ddm = get(points_odt[pair], (direction, :DDM, istp), nothing)
    odt_dis = get(points_odt[pair], (direction, :DIS, istp), nothing)
    cmot_ddm = points_cmot[pair][(:DDM, istp)]
    cmot_dis = points_cmot[pair][(:DIS, istp)]
    all(point -> !isnothing(point) && isfinite(point.num) && point.num > 0,
        (odt_ddm, odt_dis, cmot_ddm, cmot_dis)) || return nothing
    value = (odt_ddm.num / cmot_ddm.num) / (odt_dis.num / cmot_dis.num)
    relvar = sum((point.std / point.num)^2 for point in
        (odt_ddm, cmot_ddm, odt_dis, cmot_dis) if isfinite(point.std))
    n_valid_errors = count(point -> isfinite(point.std),
        (odt_ddm, cmot_ddm, odt_dis, cmot_dis))
    error = n_valid_errors == 4 ? value * sqrt(relvar) : NaN
    (; value, error, odt_ddm, cmot_ddm, odt_dis, cmot_dis)
end

function draw_efficiency_ratio()
    labels = ["$pair\n$(conditions[pair].bias >= 0 ? "+" : "")$(round(conditions[pair].bias, digits=3)) $(Int(conditions[pair].tmotload))s"
        for pair in val_pair]
    fig = Figure(size=(690, 340), fontsize=11, figure_padding=4)
    ax = Axis(fig[1, 1]; xticks=(eachindex(val_pair), labels),
        xlabel="Isotope pair · tbiasmot / tmotload",
        ylabel="(ODT / CMOT) DDM / (ODT / CMOT) DIS",
        title="ODT loading efficiency · DDM / DIS",
        yticks=0:0.2:1.2, yminorticks=IntervalsBetween(2),
        dualmot_axis_kwargs(; text_size=9)...)
    ax.yminorticks = IntervalsBetween(2)
    xlims!(ax, 0.5, length(val_pair) + 0.5)
    hlines!(ax, [0.9, 1.1]; color=(:black, 0.45), linestyle=:dash, linewidth=0.8)
    for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-")),
        (idx_direction, direction) in enumerate((:x, :z))
        result = odt_cmot_efficiency_ratio(direction, pair, istp)
        isnothing(result) && continue
        pos = idx_pair + (istp == Symbol(first(split(pair, "-"))) ? -0.07 : 0.07) +
            (idx_direction == 1 ? -0.025 : 0.025)
        style = dualmot_curve_style((; loadcfg=:DDM, istp))
        options = merge(marker_style(style; markersize=7), (; marker=:hexagon))
        draw_marker_point!(ax, pos, result.value, result.error;
            marker_options=options)
    end
    ylims!(ax, 0, 1.25)
    filename = "[ODT.CMOT.comparison].[lin].[DDM-DIS-efficiency-ratio]"
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=3) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
    fig
end

fig_numbers = draw_comparison(; ratio=false)
fig_ratios = draw_comparison(; ratio=true)
fig_efficiency_ratios = draw_efficiency_ratio()

number_headers = Any["pair", "istp", "loadcfg", "ODT or CMOT", "tbiasmot",
    "tmotload (s)", "mean", "std", "n_rep", "source", "source denominator"]
ratio_headers = Any["pair", "istp", "loadcfg", "ODT/CMOT", "tbiasmot",
    "tmotload (s)", "ratio", "std", "n_rep numerator", "n_rep denominator",
    "source numerator", "source denominator"]
number_rows = Vector{Vector{Any}}()
ratio_rows = Vector{Vector{Any}}()
efficiency_rows = Vector{Vector{Any}}()
for pair in val_pair, loadcfg in (:DDM, :DIS), istp in Symbol.(split(pair, "-"))
    condition = conditions[pair]
    mot = points_cmot[pair][(loadcfg, istp)]
    for direction in (:x, :z)
        odt = get(points_odt[pair], (direction, loadcfg, istp), nothing)
        push!(number_rows, Any[pair, string(istp), string(loadcfg), "ODT $direction",
            condition.bias, condition.tmotload,
            isnothing(odt) ? NaN : odt.num, isnothing(odt) ? NaN : odt.std,
            isnothing(odt) || odt.n_rep == 0 ? NaN : odt.n_rep,
            isnothing(odt) ? "" : odt.source, ""])
        valid_ratio = !isnothing(odt) && isfinite(odt.num) &&
            isfinite(mot.num) && mot.num > 0
        ratio = valid_ratio ? odt.num / mot.num : NaN
        err = valid_ratio && isfinite(odt.std) && isfinite(mot.std) ?
            sqrt((odt.std / mot.num)^2 + (odt.num * mot.std / mot.num^2)^2) : NaN
        push!(ratio_rows, Any[pair, string(istp), string(loadcfg), "ODT/CMOT",
            condition.bias, condition.tmotload, ratio, err,
            isnothing(odt) || odt.n_rep == 0 ? NaN : odt.n_rep,
            mot.n_rep == 0 ? NaN : mot.n_rep,
            isnothing(odt) ? "" : odt.source, mot.source])
    end
    for direction in (:x, :z)
        efficiency = odt_cmot_efficiency_ratio(direction, pair, istp)
        push!(efficiency_rows, Any[pair, string(istp), string(direction),
            "(ODT/CMOT) DDM/DIS", condition.bias, condition.tmotload,
            isnothing(efficiency) ? NaN : efficiency.value,
            isnothing(efficiency) ? NaN : efficiency.error,
            isnothing(efficiency) ? "" : efficiency.odt_ddm.source,
            isnothing(efficiency) ? "" : efficiency.cmot_ddm.source,
            isnothing(efficiency) ? "" : efficiency.odt_dis.source,
            isnothing(efficiency) ? "" : efficiency.cmot_dis.source])
    end
    push!(number_rows, Any[pair, string(istp), string(loadcfg), "CMOT",
        condition.bias, condition.tmotload, mot.num, mot.std, mot.n_rep, mot.source, ""])
end
number_matrix = Matrix{Any}(undef, length(number_rows) + 1, length(number_headers))
number_matrix[1, :] .= number_headers
for (idx, row) in enumerate(number_rows)
    number_matrix[idx + 1, :] .= row
end
ratio_matrix = Matrix{Any}(undef, length(ratio_rows) + 1, length(ratio_headers))
ratio_matrix[1, :] .= ratio_headers
for (idx, row) in enumerate(ratio_rows)
    ratio_matrix[idx + 1, :] .= row
end
write_figure_workbook(joinpath(path_output, "ODT.xlsx"), [
    "numbers" => number_matrix,
    "ODT-CMOT ratio" => ratio_matrix,
    "DDM-DIS efficiency ratio" => begin
        headers = Any["pair", "istp", "ODT direction", "quantity", "tbiasmot",
            "tmotload (s)", "ratio", "std", "source ODT DDM", "source CMOT DDM",
            "source ODT DIS", "source CMOT DIS"]
        matrix = Matrix{Any}(undef, length(efficiency_rows) + 1, length(headers))
        matrix[1, :] .= headers
        for (idx, row) in enumerate(efficiency_rows)
            matrix[idx + 1, :] .= row
        end
        matrix
    end,
])
println("Saved ODT/CMOT comparison figures to $path_output")
