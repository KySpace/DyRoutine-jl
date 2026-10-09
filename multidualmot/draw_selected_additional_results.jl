"""Draw supplementary selected isotope-pair and isotope-curve Result figures."""

path_result_additional = joinpath(path_root, "Result")
mkpath(path_result_additional)
result_additional_axis_font = "NewComputerModern Math"
result_additional_axis_label_size = 7 / 0.75
result_additional_formats = ("svg", "png", "pdf")

function save_selected_result_formats(fig::Figure, basename::AbstractString)
    for format in result_additional_formats
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_result_additional, "$basename.$format"), fig;
            save_options...)
    end
    fig
end

function selected_result_pair_axis(fig::Figure, slot; ylabel,
    log_y::Bool=false, ratio::Bool=false, yticks=nothing)
    kwargs = dualmot_axis_kwargs(; log_y, text_size=8)
    isnothing(yticks) || (kwargs = merge(kwargs, (; yticks)))
    ratio && (kwargs = merge(kwargs, (; yminorticks=IntervalsBetween(2),
        yminorticksvisible=true)))
    kwargs = merge(kwargs, (; xlabelsize=result_additional_axis_label_size,
        ylabelsize=result_additional_axis_label_size,
        xlabelfont=result_additional_axis_font,
        ylabelfont=result_additional_axis_font,
        xticklabelfont=result_additional_axis_font,
        yticklabelfont=result_additional_axis_font))
    Axis(fig[slot...]; xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair", ylabel, yscale=log_y ? log10 : identity,
        aspect=AxisAspect(320 / 149), kwargs...)
end

function plot_pair_ratio!(ax::Axis, points_by_pair::AbstractDict,
    numerator::Symbol, denominator::Symbol; value_key::Symbol, error_key::Symbol,
    marker_loadcfg::Symbol, reference_lines=(0.9, 1.1))
    draw_pair_spans!(ax)
    hlines!(ax, collect(reference_lines); color=RGBAf(0, 0, 0, 0.45),
        linestyle=:dash, linewidth=0.5)
    for (idx_pair, pair) in enumerate(val_pair), istp in Symbol.(split(pair, "-"))
        haskey(points_by_pair, pair) || continue
        pair_points = points_by_pair[pair]
        haskey(pair_points, (numerator, istp)) &&
            haskey(pair_points, (denominator, istp)) || continue
        top = pair_points[(numerator, istp)]
        bottom = pair_points[(denominator, istp)]
        top_value = Float64(getproperty(top, value_key))
        bottom_value = Float64(getproperty(bottom, value_key))
        isfinite(top_value) && top_value >= 0 &&
            isfinite(bottom_value) && bottom_value > 0 || continue
        ratio_value = top_value / bottom_value
        top_error = Float64(getproperty(top, error_key))
        bottom_error = Float64(getproperty(bottom, error_key))
        error = isfinite(top_error) && isfinite(bottom_error) ?
            sqrt((top_error / bottom_value)^2 +
                (top_value * bottom_error / bottom_value^2)^2) : NaN
        style = dualmot_curve_style((; loadcfg=marker_loadcfg, istp))
        options = merge(marker_style(style; markersize=6), (; marker=:hexagon))
        x = pair_istp_position(idx_pair, pair, istp)
        if isfinite(error)
            marker_errorbars!(ax, [x], [ratio_value], [error]; options...,
                errorlinewidth=0.75)
        else
            scatter!(ax, [x], [ratio_value]; options...)
        end
    end
    ylims!(ax, 0, 1.25)
    ax
end

# MOT 626 DCS/SCS ratio and 421 interrupt DIS/DCS ratio.
fig_selected_number_ratios = Figure(size=(fig_width, 192), fontsize=8,
    figure_padding=3, colgap=6)
ax_selected_626_ratio = selected_result_pair_axis(fig_selected_number_ratios,
    (1, 1); ylabel=result_ratio_label("𝑁", "626"), yticks=0:0.2:1.2,
    ratio=true)
plot_pair_ratio!(ax_selected_626_ratio, points_626, :DCS, :SCS;
    value_key=:num, error_key=:std, marker_loadcfg=:DCS)
ax_selected_intr_ratio = selected_result_pair_axis(fig_selected_number_ratios,
    (1, 2); ylabel=result_ratio_label("𝑁", "intr"), yticks=0:0.2:1.2,
    ratio=true)
plot_pair_ratio!(ax_selected_intr_ratio, points_421[:t_balanced], :DIS, :DCS;
    value_key=:num, error_key=:std, marker_loadcfg=:DIS)
save_selected_result_formats(fig_selected_number_ratios,
    "selected_isotope_pair_number_ratios")

# N-balanced MOT lifetime values and DDM/DIS ratios.
fig_selected_mot_lifetime = Figure(size=(fig_width, 192), fontsize=8,
    figure_padding=3, colgap=6)
points_selected_mot_lifetime = comparison_points_lifetime[:n_balanced].mot
ax_selected_mot_lifetime = selected_result_pair_axis(fig_selected_mot_lifetime,
    (1, 1); ylabel=rich("𝜏", subscript("𝑖"), superscript("CMOT"), " (s)"),
    log_y=true)
draw_pair_spans!(ax_selected_mot_lifetime)
ylims!(ax_selected_mot_lifetime, 0.9, 55)
for (idx_pair, pair) in enumerate(val_pair), loadcfg in (:DIS, :DDM),
    istp in Symbol.(split(pair, "-"))
    haskey(points_selected_mot_lifetime, pair) || continue
    point = get(points_selected_mot_lifetime[pair], (loadcfg, istp), nothing)
    isnothing(point) && continue
    value, error = Float64(point.value), Float64(point.std)
    isfinite(value) && value > 0 || continue
    style = dualmot_curve_style((; loadcfg, istp))
    options = marker_style(style; markersize=6)
    x = pair_istp_position(idx_pair, pair, istp)
    if isfinite(error) && error >= 0
        marker_errorbars_log!(ax_selected_mot_lifetime, [x], [value], [error];
            floor=0.9, options..., errorlinewidth=0.75)
    else
        scatter!(ax_selected_mot_lifetime, [x], [value]; options...)
    end
end
ax_selected_mot_lifetime_ratio = selected_result_pair_axis(fig_selected_mot_lifetime,
    (1, 2); ylabel=result_ratio_label("𝜏", "CMOT"), yticks=0:0.2:1.2,
    ratio=true)
plot_pair_ratio!(ax_selected_mot_lifetime_ratio, points_selected_mot_lifetime,
    :DDM, :DIS; value_key=:value, error_key=:std, marker_loadcfg=:DDM)
save_selected_result_formats(fig_selected_mot_lifetime,
    "selected_isotope_pair_mot_lifetime")

# ODT/CMOT atom-number efficiency, keeping the better of the x/z scans.
fig_selected_odt_efficiency = Figure(size=(fig_width, 340), fontsize=8,
    figure_padding=4)
odt_efficiency_axis_width = round(Int, 120 / 25.4 * 96)
odt_efficiency_axis_height = round(Int, odt_efficiency_axis_width / (320 / 149))
ax_selected_odt_efficiency = Axis(fig_selected_odt_efficiency[1, 1];
    xticks=(eachindex(val_pair), val_pair), xlabel="Isotope pair",
    ylabel=rich("𝜂"), yticks=0:0.1:0.5,
    xlabelsize=result_additional_axis_label_size,
    ylabelsize=result_additional_axis_label_size,
    xlabelfont=result_additional_axis_font,
    ylabelfont=result_additional_axis_font,
    xticklabelfont=result_additional_axis_font,
    yticklabelfont=result_additional_axis_font,
    width=odt_efficiency_axis_width, height=odt_efficiency_axis_height,
    halign=:center, valign=:center, tellwidth=false, tellheight=false,
    dualmot_axis_kwargs(; text_size=8)...)
draw_pair_spans!(ax_selected_odt_efficiency)
for (idx_pair, pair) in enumerate(val_pair), loadcfg in (:DDM, :DIS),
    istp in Symbol.(split(pair, "-"))
    mot_point = get(points_cmot[pair], (loadcfg, istp), nothing)
    isnothing(mot_point) && continue
    mot_point.num > 0 && isfinite(mot_point.num) || continue
    direction_points = filter(point -> !isnothing(point) &&
        isfinite(point.num) && point.num >= 0, [
        get(points_odt[pair], (direction, loadcfg, istp), nothing)
        for direction in (:x, :z)])
    isempty(direction_points) && continue
    ratios = [point.num / mot_point.num for point in direction_points]
    idx_best = findmax(ratios)[2]
    odt_point = direction_points[idx_best]
    value = ratios[idx_best]
    isfinite(value) || continue
    error = isfinite(odt_point.std) && isfinite(mot_point.std) ?
        sqrt((odt_point.std / mot_point.num)^2 +
            (odt_point.num * mot_point.std / mot_point.num^2)^2) : NaN
    x = pair_istp_position(idx_pair, pair, istp) +
        (loadcfg == :DDM ? -0.025 : 0.025)
    style = dualmot_curve_style((; loadcfg, istp))
    options = marker_style(style; markersize=6)
    if isfinite(error)
        marker_errorbars!(ax_selected_odt_efficiency, [x], [value], [error];
            options..., errorlinewidth=0.75)
    else
        scatter!(ax_selected_odt_efficiency, [x], [value]; options...)
    end
end
ylims!(ax_selected_odt_efficiency, 0, 0.56)
save_selected_result_formats(fig_selected_odt_efficiency,
    "selected_isotope_pair_odt_efficiency")

function collect_result_curves(runinfo::NamedTuple, x_key::Symbol;
    bias::Union{Nothing,Real}=nothing,
    bounds_sigmax_num::Tuple{<:Real,<:Real},
    bounds_sigmay_num::Tuple{<:Real,<:Real},
    num_max_num::Real)
    stats_blocks = [calc_num_evol_block(data;
        label="$(runinfo.folder) $(runinfo.tag) data[$idx]",
        bounds_sigmax_num, bounds_sigmay_num, num_max_num)
        for (idx, data) in enumerate(runinfo.data)]
    combined = combine_num_evol_blocks(stats_blocks;
        label="$(runinfo.folder) $(runinfo.tag)")
    idx_rep = findfirst(==(:rep), keys(combined.vars))
    nums = combined.num_fmt
    num_stat = dropdims(mapslices(nums; dims=idx_rep) do values
        valid = collect(skipmissing(vec(values)))
        isempty(valid) ? NaN : mean(valid)
    end; dims=idx_rep)
    std_stat = dropdims(mapslices(nums; dims=idx_rep) do values
        valid = collect(skipmissing(vec(values)))
        length(valid) < 2 ? NaN : std(valid)
    end; dims=idx_rep)
    name_stat = Tuple(key for key in keys(combined.vars) if key != :rep)
    idx_axis = Dict(key => idx for (idx, key) in enumerate(name_stat))
    idx_bias = if :β_MOT in name_stat
        isnothing(bias) && throw(ArgumentError("$(runinfo.folder): β_MOT selection is required"))
        findfirst(==(bias), combined.vars.β_MOT)
    else
        isnothing(bias) || throw(ArgumentError("$(runinfo.folder): unexpected β_MOT selection"))
        nothing
    end
    curves = Dict{Tuple{Symbol,Symbol},NamedTuple}()
    for (idx_cfg, loadcfg) in enumerate(combined.vars.loadcfg),
        (idx_istp, istp) in enumerate(combined.vars.istp)
        indices = Any[Colon() for _ in name_stat]
        isnothing(idx_bias) || (indices[idx_axis[:β_MOT]] = idx_bias)
        indices[idx_axis[:loadcfg]] = idx_cfg
        indices[idx_axis[:istp]] = idx_istp
        curves[(loadcfg, istp)] = (
            nums=vec(@view num_stat[indices...]),
            stds=vec(@view std_stat[indices...]),
        )
    end
    (; x=Float64.(getproperty(combined.vars, x_key)), curves)
end

function validate_result_monofreq(vars::NamedTuple, folder::AbstractString)
    pair_isotopes = Symbol.(split(folder, "-"))
    length(pair_isotopes) == 2 || throw(ArgumentError("$folder: expected an isotope pair"))
    all(in(pair_isotopes), vars.istp) || throw(ArgumentError("$folder: invalid isotope values"))
    all(in((:SDS, :SIS)), vars.loadcfg) || throw(ArgumentError("$folder: expected SDS/SIS"))
    nothing
end

function latest_monofreq_isotope_run(pair::AbstractString, isotope::Symbol)
    infos = read_num_evol_runinfos(joinpath(path_root, "MOT loading 421 monofreq"), pair;
        var_specs=DUALMOT_LOADING_VAR_SPECS,
        validate_vars=validate_result_monofreq)
    candidates = filter(info -> any(data -> isotope in data.vars.istp,
        info.data), infos)
    runinfo = latest_num_evol_runinfo(candidates;
        label="$pair MOT loading 421 monofreq for isotope $isotope")
    biases = sort!(unique(vcat([data.vars.β_MOT for data in runinfo.data
        if isotope in data.vars.istp]...)))
    bias = if length(biases) == 1
        only(biases)
    elseif 0.0 in biases
        0.0
    else
        throw(ArgumentError(
            "$pair $isotope monofrequency curve has multiple bias values: $biases"))
    end
    (; runinfo, bias)
end

function result_curve_axis(fig::Figure, slot; ylabel, log_y=false,
    xlabel, xlim, xticks, xminorticks=nothing, yticks=nothing,
    aspect=AxisAspect(4 / 3))
    options = dualmot_axis_kwargs(; log_y, text_size=8)
    options = merge(options, (; xlabelsize=result_additional_axis_label_size,
        ylabelsize=result_additional_axis_label_size,
        xlabelfont=result_additional_axis_font,
        ylabelfont=result_additional_axis_font,
        xticklabelfont=result_additional_axis_font,
        yticklabelfont=result_additional_axis_font,
        xticks, xminorticks, xminorticksvisible=!isnothing(xminorticks)))
    isnothing(yticks) || (options = merge(options, (; yticks)))
    Axis(fig[slot...]; xlabel, ylabel, yscale=log_y ? log10 : identity,
        limits=(xlim, nothing), aspect, options...)
end

function plot_loading_result_panel!(ax::Axis, curve_data::NamedTuple,
    pair_label::AbstractString, configs, isotopes; ylabel::Bool=false)
    ax.title = pair_label
    ax.titlesize = result_additional_axis_label_size
    x_raw = hasproperty(curve_data, :x) ? curve_data.x : curve_data.val_t_load
    x_values = x_raw .* 0.43
    plotted_values = Float64[]
    for loadcfg in configs, isotope in isotopes
        curve = get(curve_data.curves, (loadcfg, isotope), nothing)
        isnothing(curve) && continue
        y_values = curve.nums ./ 1e7
        valid = isfinite.(x_values) .& isfinite.(y_values)
        any(valid) || continue
        style = dualmot_curve_style((; loadcfg, istp=isotope))
        lines!(ax, x_values[valid], y_values[valid]; style.line_options...)
        scatter!(ax, x_values[valid], y_values[valid];
            marker_style(style; markersize=4)...)
        append!(plotted_values, y_values[valid])
    end
    isempty(plotted_values) && throw(ArgumentError("$pair_label: no finite loading data"))
    ylims!(ax, 0, max(1.0, 1.08 * maximum(plotted_values)))
    ylabel && (ax.ylabel = rich("𝑁", subscript("𝑖"), superscript("CMOT"),
        " (10", superscript("7"), ")"))
    ax
end

# Four loading curves, with the fourth panel centered in the second row.
fig_selected_loading_pairs = Figure(size=(fig_width, 424), fontsize=8,
    figure_padding=4, rowgap=10, colgap=8)
selected_loading_axes = Axis[]
loading_axis_options = (; xlabel="Effective loading time (s)", xlim=(0.0, 20.0),
    xticks=0:5:20, xminorticks=setdiff(collect(0:1:20), collect(0:5:20)))
for (slot, pair_label) in zip(((1, 1), (1, 2), (1, 3), (2, 2)),
    ("162–163", "161–163", "163–164", "161–163"))
    ax = result_curve_axis(fig_selected_loading_pairs, slot;
        ylabel="", loading_axis_options...)
    push!(selected_loading_axes, ax)
    if length(selected_loading_axes) == 1
        selection = latest_monofreq_isotope_run("162-163", Symbol("163"))
        mono_data = collect_result_curves(selection.runinfo, :t_load;
            bias=selection.bias,
            bounds_sigmax_num=(4e-4, Inf), bounds_sigmay_num=(4e-4, Inf),
            num_max_num=2e8)
        plot_loading_result_panel!(ax, mono_data, pair_label, (:SIS, :SDS),
            (Symbol("163"),); ylabel=true)
    elseif length(selected_loading_axes) == 2
        mono_curves = Dict{Tuple{Symbol,Symbol},NamedTuple}()
        mono_x = nothing
        for isotope in (Symbol("161"), Symbol("163"))
            selection = latest_monofreq_isotope_run("161-163", isotope)
            mono_data = collect_result_curves(selection.runinfo, :t_load;
                bias=selection.bias,
                bounds_sigmax_num=(4e-4, Inf), bounds_sigmay_num=(4e-4, Inf),
                num_max_num=2e8)
            isnothing(mono_x) && (mono_x = mono_data.x)
            mono_x == mono_data.x || throw(ArgumentError(
                "161-163 monofrequency isotope curves have incompatible loading-time axes"))
            for loadcfg in (:SIS, :SDS)
                mono_curves[(loadcfg, isotope)] = mono_data.curves[(loadcfg, isotope)]
            end
        end
        mono_data = (; x=mono_x, curves=mono_curves)
        plot_loading_result_panel!(ax, mono_data, pair_label, (:SIS, :SDS),
            (Symbol("161"), Symbol("163")))
    elseif length(selected_loading_axes) == 3
        runinfo_626 = latest_num_evol_runinfo(runinfos_626["163-164"];
            label="163-164 MOT loading 626 Result")
        data_626 = collect_loading_curves(runinfo_626;
            bounds_sigmax_num=(4e-4, Inf), bounds_sigmay_num=(4e-4, Inf),
            num_max_num=2e8)
        plot_loading_result_panel!(ax, data_626, pair_label, (:SCS, :DCS),
            (Symbol("163"), Symbol("164")))
    else
        runinfo_626 = latest_num_evol_runinfo(runinfos_626["161-163"];
            label="161-163 MOT loading 626 Result")
        data_626 = collect_loading_curves(runinfo_626;
            bounds_sigmax_num=(4e-4, Inf), bounds_sigmay_num=(4e-4, Inf),
            num_max_num=2e8)
        plot_loading_result_panel!(ax, data_626, pair_label, (:SCS, :DCS),
            (Symbol("161"), Symbol("163")))
    end
end
save_selected_result_formats(fig_selected_loading_pairs,
    "selected_isotope_pair_loading_curves")

function result_decay_fit_record(pair::AbstractString, bias::Real,
    loadcfg::Symbol, isotope::Symbol)
    selected = filter(latest_cmot_decay.records) do record
        record.pair == pair && record.panel == Float64(bias) &&
            record.loadcfg == loadcfg && record.istp == isotope
    end
    length(selected) == 1 || throw(ArgumentError(
        "$pair $loadcfg $isotope at β_MOT=$bias: expected one CMOT fit record, found $(length(selected))"))
    only(selected)
end

function result_cmot_runinfo(pair::AbstractString, tag::AbstractString)
    infos = read_num_evol_runinfos(joinpath(path_root, "CMOT lifetime"), pair;
        var_specs=DUALMOT_VAR_SPECS, validate_vars=validate_dualmot_vars)
    matches = filter(info -> info.tag == tag, infos)
    length(matches) == 1 || throw(ArgumentError("$pair: expected one CMOT run with tag $tag"))
    only(matches)
end

fig_selected_cmot_curves = Figure(size=(fig_width, 240), fontsize=8,
    figure_padding=4, colgap=8)
selected_cmot_pairs = ("162-163", "163-164", "161-163")
selected_cmot_ylims = Dict(
    "162-163" => (1e6, 5e7),
    "163-164" => (1e6, 4e6),
    "161-163" => (7e5, 1e7),
)
for (idx_pair, pair) in enumerate(selected_cmot_pairs)
    isotope_values = Symbol.(split(pair, "-"))
    fit_bias = balance_biases[:n_balanced][pair]
    fit_tags = unique([result_decay_fit_record(pair, fit_bias, :DDM, isotope).tag
        for isotope in isotope_values])
    length(fit_tags) == 1 || throw(ArgumentError("$pair: CMOT fit tags disagree: $fit_tags"))
    runinfo_cmot = result_cmot_runinfo(pair, only(fit_tags))
    data_cmot = collect_result_curves(runinfo_cmot, :t_hold; bias=fit_bias,
        bounds_sigmax_num=(4e-4, Inf), bounds_sigmay_num=(4e-4, Inf),
        num_max_num=2e8)
    ax = result_curve_axis(fig_selected_cmot_curves, (1, idx_pair);
        ylabel=idx_pair == 1 ? rich("𝑁", subscript("𝑖"), superscript("CMOT")) : "",
        xlabel="CMOT holding time (s)", xlim=(0.0, 1.05),
        xticks=0:0.2:1.0, log_y=true)
    ax.title = replace(pair, "-" => "–")
    ax.titlesize = result_additional_axis_label_size
    for loadcfg in (:DIS, :DDM), isotope in isotope_values
        curve = get(data_cmot.curves, (loadcfg, isotope), nothing)
        isnothing(curve) && continue
        x_values = data_cmot.x ./ 1000
        mask = (x_values .>= 0) .& (x_values .<= 1.05) .&
            isfinite.(curve.nums) .& (curve.nums .> 0)
        any(mask) || continue
        style = dualmot_curve_style((; loadcfg, istp=isotope))
        scatter!(ax, x_values[mask], curve.nums[mask];
            marker_style(style; markersize=4)...)
        fit_record = result_decay_fit_record(pair, fit_bias, loadcfg, isotope)
        fit_x = range(0.0, 1.05; length=250)
        fit_y = model_num_decay_kappa(collect(fit_x), [fit_record.n0, fit_record.kappa])
        lines!(ax, fit_x, fit_y; style.line_options...)
    end
    y_limits = selected_cmot_ylims[pair]
    ylims!(ax, y_limits...)
    ax.yticks = LogTicks(floor(Int, log10(first(y_limits))):
        floor(Int, log10(last(y_limits))))
end
save_selected_result_formats(fig_selected_cmot_curves,
    "selected_isotope_pair_cmot_decay_curves")
