include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root_421 = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading 421"
path_root_626 = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading 626"
path_output = joinpath(dirname(path_root_421), "MOT loading pair comparison")
val_pair = ["162-164", "161-162", "162-163", "161-164", "163-164",
    "161-163", "160-162"]
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
    loadcfgs::Tuple, val_istp::AbstractVector{Symbol}, label::AbstractString)
    idx_final = lastindex(loading_curves.val_t_load)
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

points_421 = Dict{String, Dict{Tuple{Symbol, Symbol}, NamedTuple}}()
points_626 = Dict{String, Dict{Tuple{Symbol, Symbol}, NamedTuple}}()
for pair in val_pair
    candidates_421 = filter(runinfos_421[pair]) do runinfo
        length(runinfo.data) == 1 && 0.0 in only(runinfo.data).vars.β_MOT
    end
    length(candidates_421) == 1 ||
        throw(ArgumentError("$pair: expected exactly one 421 processing entry containing β_MOT = 0"))
    runinfo_421 = only(candidates_421)
    runinfo_626 = only(runinfos_626[pair])
    val_istp = Symbol.(split(pair, "-"))

    curves_421 = collect_loading_curves(runinfo_421;
        bias=0.0, bounds_sigmax_num, bounds_sigmay_num, num_max_num)
    curves_626 = collect_loading_curves(runinfo_626;
        bounds_sigmax_num, bounds_sigmay_num, num_max_num)
    points_421[pair] = final_loading_points(
        curves_421, (:DDM, :DIS), val_istp, "$pair MOT loading 421")
    points_626[pair] = final_loading_points(
        curves_626, (:DCS, :SCS), val_istp, "$pair MOT loading 626")
end

marker_style(style::NamedTuple; markersize::Real=style.markersize) = (
    color=style.markercolor,
    marker=style.marker,
    markersize,
    strokecolor=style.strokecolor,
    strokewidth=style.strokewidth,
)

function draw_style_key!(fig::Figure, loadcfgs::Tuple)
    val_istp_key = Symbol.(string.(160:164))
    pos_y = collect(5:-1:1)
    ax_key = Axis(fig[1, 1];
        width=125,
        height=155,
        halign=0.04,
        valign=:bottom,
        tellwidth=false,
        tellheight=false,
        limits=(0.5, 2.5, 0.5, 5.5),
        xticks=(Float64.(1:length(loadcfgs)), string.(collect(loadcfgs))),
        yticks=(pos_y, string.(val_istp_key)),
        xaxisposition=:top,
        xgridvisible=false,
        ygridvisible=false,
        xticksize=0,
        yticksize=0,
        xticklabelsize=10,
        yticklabelsize=10,
        xticklabelpad=3,
        yticklabelpad=3,
        backgroundcolor=:white,
    )
    for (idx_loadcfg, loadcfg) in enumerate(loadcfgs),
        (idx_istp, istp) in enumerate(val_istp_key)
        style = dualmot_curve_style((; loadcfg, istp))
        scatter!(ax_key, [idx_loadcfg], [pos_y[idx_istp]];
            marker_style(style; markersize=10)...)
    end
    ax_key
end

function draw_pair_ratio(points_by_pair::AbstractDict, numerator::Symbol, denominator::Symbol;
    ylabel::AbstractString, title::AbstractString, filename::AbstractString)
    fig = Figure(size=(900, 420))
    ax = Axis(fig[1, 1];
        xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair",
        ylabel,
        title,
        yminorticks=IntervalsBetween(5),
        yminorticksvisible=true,
    )
    for (idx_pair, pair) in enumerate(val_pair),
        (idx_istp, istp) in enumerate(Symbol.(split(pair, "-")))
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
        x = idx_pair
        style = dualmot_curve_style((; loadcfg=numerator, istp))
        scatter!(ax, [x], [ratio]; marker_style(style; markersize=17)...)
        isfinite(std_ratio) && errorbars!(ax, [x], [ratio], [std_ratio];
            color=style.color, whiskerwidth=7, linewidth=1.2)
    end
    ylims!(ax, 0, nothing)
    draw_style_key!(fig, (numerator, denominator))
    for format in formats_output
        save(joinpath(path_output, "$filename.$format"), fig)
    end
    fig
end

function draw_pair_numbers(points_by_pair::AbstractDict, loadcfgs::Tuple;
    ylabel::AbstractString, title::AbstractString, filename::AbstractString)
    fig = Figure(size=(980, 420))
    ax = Axis(fig[1, 1];
        xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair",
        ylabel,
        title,
        yscale=log10,
        yminorticks=IntervalsBetween(5),
        yminorticksvisible=true,
    )
    # Draw filled markers first and open markers second so coincident loadcfgs remain visible.
    loadcfgs_draw = sort(collect(loadcfgs); by=loadcfg -> loadcfg in (:DDM, :SCS) ? 0 : 1)
    for (idx_pair, pair) in enumerate(val_pair)
        val_istp = Symbol.(split(pair, "-"))
        for loadcfg in loadcfgs_draw, istp in val_istp
            key = (loadcfg, istp)
            point = points_by_pair[pair][key]
            isfinite(point.num) && point.num > 0 ||
                throw(ArgumentError("$pair $key: number must be finite and positive"))
            x = idx_pair
            style = dualmot_curve_style((; loadcfg, istp))
            scatter!(ax, [x], [point.num]; marker_style(style; markersize=17)...)
            isfinite(point.std) && point.num - point.std > 0 &&
                errorbars!(ax, [x], [point.num], [point.std];
                    color=style.color, whiskerwidth=7, linewidth=1.2)
        end
    end
    draw_style_key!(fig, loadcfgs)
    for format in formats_output
        save(joinpath(path_output, "$filename.$format"), fig)
    end
    fig
end

mkpath(path_output)
fig_ratio_ddm_dis = draw_pair_ratio(points_421, :DDM, :DIS;
    ylabel="N_DDM / N_DIS",
    title="MOT loading 421 · final acquired point · β_MOT = 0",
    filename="[MOT.loading.pairs].[DDM-DIS].[ratio]",
)
fig_ratio_dcs_scs = draw_pair_ratio(points_626, :DCS, :SCS;
    ylabel="N_DCS / N_SCS",
    title="MOT loading 626 · final acquired point",
    filename="[MOT.loading.pairs].[DCS-SCS].[ratio]",
)
fig_nums_ddm_dis = draw_pair_numbers(points_421, (:DDM, :DIS);
    ylabel="MOT number",
    title="MOT loading 421 · final acquired point · β_MOT = 0",
    filename="[MOT.loading.pairs].[DDM-DIS].[nums]",
)
fig_nums_dcs_scs = draw_pair_numbers(points_626, (:DCS, :SCS);
    ylabel="MOT number",
    title="MOT loading 626 · final acquired point",
    filename="[MOT.loading.pairs].[DCS-SCS].[nums]",
)
