include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root_odt = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\ODT BField"
path_root_421 = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading 421"
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
    (; vars=combined.vars, name=Tuple(k for k in keys(combined.vars) if k != :rep), num_stat, std_stat)
end

function stats_point(stats, condition::NamedTuple)
    idx = Any[Colon() for _ in stats.name]
    for key in keys(condition)
        axis_idx = findfirst(==(key), stats.name)
        isnothing(axis_idx) && throw(ArgumentError("missing axis $key in $(stats.name)"))
        idx[axis_idx] = findfirst(==(getproperty(condition, key)), getproperty(stats.vars, key))
        isnothing(idx[axis_idx]) && throw(ArgumentError("condition $key=$(getproperty(condition,key)) unavailable"))
    end
    (; num=only(stats.num_stat[idx...]), std=only(stats.std_stat[idx...]))
end

function draw_marker_point!(ax, x, y, err; marker_options, label=nothing)
    if isfinite(err) && err >= 0
        marker_errorbars!(ax, [x], [y], [err]; marker_options..., errorlinewidth=0.8, label)
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
                    strokecolor=style.strokecolor, strokewidth=style.strokewidth))
            end
            hollow = (; color=:transparent, marker=style.marker, markersize=7,
                strokecolor=style.strokecolor, strokewidth=style.strokewidth)
            draw_marker_point!(ax, pos, mot.num, mot.std; marker_options=hollow)
        end
    end
    filename = ratio ? "[ODT.CMOT.comparison].[lin].[ratios]" : "[ODT.CMOT.comparison].[log].[numbers]"
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=3) : (;)
        save(joinpath(path_output,"$filename.$format"), fig; save_options...)
    end
    fig
end

fig_numbers = draw_comparison(; ratio=false)
fig_ratios = draw_comparison(; ratio=true)
println("Saved ODT/CMOT comparison figures to $path_output")
