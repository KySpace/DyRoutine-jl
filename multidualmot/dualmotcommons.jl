using YAML
using MAT
using Statistics: mean, std
using CairoMakie
using Colors: Oklch, RGB
using Dates
using JLD2

isdefined(@__MODULE__, :marker_errorbars!) ||
    include(joinpath(@__DIR__, "..", "snippets", "marker_errorbars.jl"))

const DUALMOT_VAR_SPECS = (
    β_MOT=(config="tbiasmot", convert=values -> Float64.(values)),
    t_hold=(config="t_hold", convert=values -> Float64.(values)),
    loadcfg=(config="loadcfg", convert=values -> Symbol.(string.(values))),
    istp=(config="istp", convert=values -> Symbol.(string.(values))),
)

const DUALMOT_LOADING_VAR_SPECS = (
    β_MOT=(config="tbiasmot", convert=values -> Float64.(values)),
    t_load=(config="t_load", convert=values -> Float64.(values)),
    loadcfg=(config="loadcfg", convert=values -> Symbol.(string.(values))),
    istp=(config="istp", convert=values -> Symbol.(string.(values))),
)

const DUALMOT_LOADCFG_VAR_SPECS = (
    t_load=(config="t_load", convert=values -> Float64.(values)),
    loadcfg=(config="loadcfg", convert=values -> Symbol.(string.(values))),
    istp=(config="istp", convert=values -> Symbol.(string.(values))),
)

const DUALMOT_BALANCE_VAR_SPECS = (
    β_MOT=(config="tbiasmot", convert=values -> Float64.(values)),
    loadcfg=(config="loadcfg", convert=values -> Symbol.(string.(values))),
    istp=(config="istp", convert=values -> Symbol.(string.(values))),
)

const DUALMOT_ODT_BFIELD_VAR_SPECS = (
    ib=(config="ib", convert=values -> Float64.(values)),
    loadcfg=(config="loadcfg", convert=values -> Symbol.(string.(values))),
    istp=(config="istp", convert=values -> Symbol.(string.(values))),
)

const ODT_BFIELD_CALIBRATION = Dict(
    :x => (Q=1.48, B0=0.28),
    :z => (Q=1.72, B0=-0.35),
)

function odt_bfield_values(current::AbstractVector{<:Real}, direction::Symbol)
    calibration = get(ODT_BFIELD_CALIBRATION, direction, nothing)
    isnothing(calibration) &&
        throw(ArgumentError("ODT B-field direction must be :x or :z, got $direction"))
    calibration.Q .* current .- calibration.B0
end

function odt_bfield_xlabel(direction::Symbol)
    haskey(ODT_BFIELD_CALIBRATION, direction) ||
        throw(ArgumentError("ODT B-field direction must be :x or :z, got $direction"))
    rich(rich("B", font=:italic),
        subscript(rich(string(direction), font=:italic)), " (A)")
end

function odt_bfield_axis_options(current::AbstractVector{<:Real}, direction::Symbol)
    values = odt_bfield_values(current, direction)
    step_major = 0.05
    idx_min = floor(Int, minimum(values) / step_major)
    idx_max = ceil(Int, maximum(values) / step_major)
    (
        xticks=collect(idx_min:idx_max) .* step_major,
        xminorticks=IntervalsBetween(5),
        xminorticksvisible=true,
    )
end

const HUE_ISTP = Dict(Symbol(string(i)) => h for (i, h) in
    ((160, 195), (161, 306), (162, 21), (163, 90), (164, 259)))
const LIGHTNESS_STROKE, CHROMA_STROKE = 0.45, 0.10
const LIGHTNESS_FACE, CHROMA_FACE = 0.85, 0.06
const LIGHTNESS_DIS_LINE, CHROMA_DIS_LINE = 0.65, 0.08
const MARKER_LOADCFG = Dict(
    :DDM => :rect, :DIS => :utriangle, :DCS => :circle, :SCS => :diamond)
const DUALMOT_LOG_MAJOR_TICKS = LogTicks(-20:20)
const DUALMOT_LOG_MINOR_TICKS = IntervalsBetween(10)
const DUALMOT_FONT = "Helvetica World"
set_theme!(fonts=(;
    regular=DUALMOT_FONT,
    bold=DUALMOT_FONT,
    italic=DUALMOT_FONT,
    bold_italic=DUALMOT_FONT,
))
#        stroke      face
# 160    #006566     #a0dbda
# 161    #634581     #d7c4ee
# 162    #843b3d     #f3bfbd
# 163    #6b5200     #ddcda1
# 164    #31558c     #b7cff6

function validate_dualmot_vars(vars::NamedTuple, tag_head::AbstractString)
    pair = join(string.(vars.istp), "-")
    occursin(Regex("(?<!\\d)" * pair * "(?!\\d)"), tag_head) ||
        throw(ArgumentError("$tag_head: folder tag must contain the configured isotope pair $pair"))
    length(vars.istp) == 2 || throw(ArgumentError("$pair: expected two isotopes"))
    all(in(vars.loadcfg), (:DDM, :DIS)) ||
        throw(ArgumentError("$pair: loadcfg must include DDM and DIS"))
    for key in (:β_MOT, :loadcfg, :istp)
        values = getproperty(vars, key)
        allunique(values) || throw(ArgumentError("$pair: duplicate values in $key"))
    end
    nothing
end

function validate_odt_bfield_vars(vars::NamedTuple, folder::AbstractString)
    pair = Symbol.(split(folder, "-"))
    length(pair) == 2 || throw(ArgumentError("$folder: expected an isotope-pair folder"))
    vars.istp == pair ||
        throw(ArgumentError("$folder: isotope order must be $(collect(pair))"))
    allunique(vars.ib) || throw(ArgumentError("$folder: duplicate ib values"))
    all(isfinite, vars.ib) || throw(ArgumentError("$folder: ib values must be finite"))
    allunique(vars.loadcfg) || throw(ArgumentError("$folder: duplicate loadcfg values"))
    all(in((:DDM, :DIS, :DCS)), vars.loadcfg) ||
        throw(ArgumentError("$folder: loadcfg values must be DDM, DIS, or DCS"))
    all(in(vars.loadcfg), (:DDM, :DIS)) ||
        throw(ArgumentError("$folder: loadcfg must include DDM and DIS"))
    nothing
end

function validate_balance_vars(vars::NamedTuple, folder::AbstractString)
    pair = Symbol.(split(folder, "-"))
    length(pair) == 2 || throw(ArgumentError("$folder: expected an isotope-pair folder"))
    vars.istp == pair ||
        throw(ArgumentError("$folder: isotope order must be $(collect(pair))"))
    allunique(vars.β_MOT) || throw(ArgumentError("$folder: duplicate β_MOT values"))
    all(isfinite, vars.β_MOT) || throw(ArgumentError("$folder: β_MOT values must be finite"))
    allunique(vars.loadcfg) || throw(ArgumentError("$folder: duplicate loadcfg values"))
    !isempty(vars.loadcfg) && all(in((:DDM, :DIS)), vars.loadcfg) ||
        throw(ArgumentError("$folder: loadcfg values must be DDM or DIS"))
    nothing
end

function validate_loadcfg_comparison_vars(vars::NamedTuple, folder::AbstractString)
    pair = Symbol.(split(folder, "-"))
    length(pair) == 2 || throw(ArgumentError("$folder: expected an isotope-pair folder"))
    allunique(vars.t_load) || throw(ArgumentError("$folder: duplicate t_load values"))
    length(vars.loadcfg) == 1 && only(vars.loadcfg) in (:DCS, :SCS) ||
        throw(ArgumentError("$folder: each data block must contain exactly one of DCS or SCS"))
    allunique(vars.istp) || throw(ArgumentError("$folder: duplicate isotope values"))
    if only(vars.loadcfg) == :DCS
        vars.istp == pair ||
            throw(ArgumentError("$folder: DCS isotope order must be $(collect(pair))"))
    else
        length(vars.istp) == 1 && only(vars.istp) in pair ||
            throw(ArgumentError("$folder: each SCS block must contain one isotope from $(collect(pair))"))
    end
    nothing
end

function read_num_evol_data(path_folder::AbstractString, config::AbstractDict;
    var_specs::NamedTuple,
    validate_vars=(vars, folder) -> nothing,
    folder::AbstractString,
    label::AbstractString,
)
    names_var = keys(var_specs)
    values_var = map(values(var_specs)) do spec
        haskey(config, spec.config) || throw(ArgumentError("$label: missing configuration $(spec.config)"))
        raw_values = config[spec.config]
        raw_values isa AbstractVector && !isempty(raw_values) ||
            throw(ArgumentError("$label: $(spec.config) must be a nonempty vector"))
        spec.convert(raw_values)
    end
    vars = merge((rep=:auto,), NamedTuple{names_var}(values_var))
    validate_vars(vars, folder)

    config_to_name = Dict(spec.config => name for (name, spec) in pairs(var_specs))
    raw_order = get(config, "vars", nothing)
    raw_order isa AbstractVector || throw(ArgumentError("$label: vars must be a vector"))
    var_order = map(raw_order) do raw_name
        raw_name == "rep" ? :rep : get(config_to_name, string(raw_name)) do
            throw(ArgumentError("$label: unknown variable in vars: $raw_name"))
        end
    end |> Tuple
    expected_order = (:rep, names_var...)
    Set(var_order) == Set(expected_order) && length(var_order) == length(expected_order) ||
        throw(ArgumentError("$label: vars must order $(join(expected_order, ", ")) exactly once"))

    sources = get(config, "source", nothing)
    sources isa AbstractVector && !isempty(sources) && all(source -> source isa AbstractString, sources) ||
        throw(ArgumentError("$label: source must be a nonempty list of filenames"))
    sources = String.(sources)
    allunique(sources) || throw(ArgumentError("$label: source filenames must be unique"))
    all(source -> basename(source) == source, sources) ||
        throw(ArgumentError("$label: source entries must be filenames in the pair folder"))
    files = joinpath.(path_folder, sources)
    all(isfile, files) || begin
        missing_files = sources[.!isfile.(files)]
        throw(ArgumentError("$label: missing source files: $(join(missing_files, ", "))"))
    end
    date_runid = map(sources) do filename
        matched = match(r"^liferes (\d{4}) run(\d+)\.mat$", filename)
        isnothing(matched) &&
            throw(ArgumentError("$label: invalid liferes source filename: $filename"))
        (String(matched[1]), parse(Int, matched[2]))
    end
    allunique(date_runid) || throw(ArgumentError("$label: duplicate date/run IDs"))
    (; date_runid, vars, var_order, files, sources)
end

function read_num_evol_runinfos(path_root::AbstractString, folder::AbstractString;
    var_specs::NamedTuple,
    validate_vars=(vars, folder) -> nothing,
)
    path_folder = joinpath(path_root, folder)
    path_config = joinpath(path_folder, "config.yaml")
    config = YAML.load_file(path_config)
    config isa AbstractVector && !isempty(config) ||
        throw(ArgumentError("$folder: config.yaml must be a nonempty top-level sequence"))
    runinfos = map(enumerate(config)) do (idx_run, run_config)
        run_config isa AbstractDict ||
            throw(ArgumentError("$folder processing[$idx_run]: entry must be a mapping"))
        tag = get(run_config, "tag", nothing)
        tag isa AbstractString && !isempty(tag) ||
            throw(ArgumentError("$folder processing[$idx_run]: tag must be a nonempty string"))
        configs_data = get(run_config, "data", nothing)
        configs_data isa AbstractVector && !isempty(configs_data) ||
            throw(ArgumentError("$folder/$tag: data must be a nonempty sequence"))
        data = map(enumerate(configs_data)) do (idx_data, config_data)
            config_data isa AbstractDict ||
                throw(ArgumentError("$folder/$tag data[$idx_data]: entry must be a mapping"))
            read_num_evol_data(path_folder, config_data;
                var_specs, validate_vars, folder,
                label="$folder/$tag data[$idx_data]")
        end
        (; folder=String(folder), tag=String(tag), data)
    end
    tags = getproperty.(runinfos, :tag)
    allunique(tags) || throw(ArgumentError("$folder: processing tags must be unique"))
    runinfos
end

function calc_num_evol_block(runinfo_data::NamedTuple;
    label::AbstractString,
    bounds_sigmax_num::Tuple{<:Real,<:Real},
    bounds_sigmay_num::Tuple{<:Real,<:Real},
    num_max_num::Real,
    allow_partial_rep::Bool=false,
)
    for (name_bound, bounds) in
        (("sigmax", bounds_sigmax_num), ("sigmay", bounds_sigmay_num))
        lower, upper = bounds
        isfinite(lower) && !isnan(upper) && lower <= upper ||
            throw(ArgumentError("$label: $name_bound bounds must satisfy finite lower <= upper, got $bounds"))
    end
    isfinite(num_max_num) && num_max_num >= 0 ||
        throw(ArgumentError("$label: num_max_num must be finite and nonnegative"))

    shot_data = read_cres_fields(runinfo_data.files, ("atomnum", "sigmax", "sigmay"))
    num_data, sigmax_data, sigmay_data = shot_data
    all(isfinite, num_data) || throw(ArgumentError("$label: atomnum must contain only finite values"))
    len_data = length(num_data)
    n_variation = prod(length(getproperty(runinfo_data.vars, key))
        for key in keys(runinfo_data.vars) if key != :rep)
    len_data > 0 || throw(DimensionMismatch("$label: atom-number data must not be empty"))
    n_partial = rem(len_data, n_variation)
    iszero(n_partial) || allow_partial_rep ||
        throw(DimensionMismatch("$label: $len_data atom numbers must be a positive integer multiple of $n_variation variations"))
    n_rep = cld(len_data, n_variation)
    n_padding = n_rep * n_variation - len_data
    if n_padding > 0
        @warn "$label: padding incomplete final repetition with missing values" samples=len_data variations_per_rep=n_variation missing_slots=n_padding
    end
    vars = merge(runinfo_data.vars, (; rep=1:n_rep))
    name = keys(vars)

    name_acq = runinfo_data.var_order
    n_dims_acq = map(key -> length(getproperty(vars, key)), name_acq)
    mask_sigmax = isfinite.(sigmax_data) .&
        (sigmax_data .>= bounds_sigmax_num[1]) .& (sigmax_data .<= bounds_sigmax_num[2])
    mask_sigmay = isfinite.(sigmay_data) .&
        (sigmay_data .>= bounds_sigmay_num[1]) .& (sigmay_data .<= bounds_sigmay_num[2])
    mask_size = mask_sigmax .& mask_sigmay
    mask_num_low = num_data .>= 0
    mask_num_high = num_data .<= num_max_num
    mask_valid = mask_size .& mask_num_low .& mask_num_high
    num_data_masked = Vector{Union{Missing, Float64}}(num_data)
    num_data_masked[.!mask_valid] .= missing
    append!(num_data_masked, fill(missing, n_padding))

    num_acq = reshape(num_data_masked, reverse(n_dims_acq)) |>
        data -> permutedims(data, reverse(1:length(n_dims_acq)))
    num_fmt = permutedims(num_acq, indexin(collect(name), collect(name_acq)))
    idx_rep_axis = findfirst(==(:rep), name)
    num_stat = dropdims(mapslices(num_fmt; dims=idx_rep_axis) do values
        valid = collect(skipmissing(vec(values)))
        isempty(valid) ? NaN : mean(valid)
    end; dims=idx_rep_axis)
    std_num_stat = dropdims(mapslices(num_fmt; dims=idx_rep_axis) do values
        valid = collect(skipmissing(vec(values)))
        length(valid) < 2 ? NaN : std(valid)
    end; dims=idx_rep_axis)
    n_rep_stat = dropdims(sum(.!ismissing.(num_fmt); dims=idx_rep_axis); dims=idx_rep_axis)
    name_stat = Tuple(key for key in name if key != :rep)

    n_masked_size = count(!, mask_size)
    n_masked_sigmax = count(!, mask_sigmax)
    n_masked_sigmay = count(!, mask_sigmay)
    n_masked_num_low = count(!, mask_num_low)
    n_masked_num_high = count(!, mask_num_high)
    n_masked_total = count(!, mask_valid)
    println("$label: $len_data samples / $n_variation variations = $n_rep repetitions; " *
        "$n_masked_total rejected ($n_masked_size by size: $n_masked_sigmax sigmax, " *
        "$n_masked_sigmay sigmay; $n_masked_num_low below zero, " *
        "$n_masked_num_high above $num_max_num)")
    (; vars, name_stat, num_fmt, num_stat, std_num_stat, n_rep_stat, n_rep,
        mask_valid, n_masked_total)
end

function combine_num_evol_blocks(blocks::AbstractVector{<:NamedTuple};
    label::AbstractString,
)
    isempty(blocks) && throw(ArgumentError("$label: at least one data block is required"))
    name = keys(first(blocks).vars)
    :rep in name || throw(ArgumentError("$label: data blocks must include a rep axis"))

    for (idx_block, block) in enumerate(blocks)
        keys(block.vars) == name ||
            throw(ArgumentError("$label data[$idx_block]: variable axes $(keys(block.vars)) must match $name"))
        expected_size = Tuple(length(getproperty(block.vars, key)) for key in name)
        size(block.num_fmt) == expected_size ||
            throw(DimensionMismatch("$label data[$idx_block]: num_fmt size $(size(block.num_fmt)) must match variable dimensions $expected_size"))
        for key in name
            values_axis = getproperty(block.vars, key)
            allunique(values_axis) ||
                throw(ArgumentError("$label data[$idx_block]: duplicate values in $key: $(collect(values_axis))"))
        end
    end

    n_rep = sum(block.n_rep for block in blocks)
    values_combined = map(name) do key
        key == :rep && return 1:n_rep
        unique(vcat((collect(getproperty(block.vars, key)) for block in blocks)...))
    end
    vars = NamedTuple{name}(values_combined)
    size_combined = Tuple(length.(values_combined))
    num_fmt = Array{Union{Missing, Float64}}(undef, size_combined)
    fill!(num_fmt, missing)

    pos_rep = 0
    for (idx_block, block) in enumerate(blocks)
        indices = map(name) do key
            if key == :rep
                (pos_rep + 1):(pos_rep + block.n_rep)
            else
                values_axis = getproperty(block.vars, key)
                values_axis_combined = getproperty(vars, key)
                positions = indexin(values_axis, values_axis_combined)
                all(!isnothing, positions) ||
                    error("$label data[$idx_block]: failed to align $key values")
                Int.(positions)
            end
        end
        @views num_fmt[indices...] .= block.num_fmt
        pos_rep += block.n_rep
    end
    (; vars, num_fmt, n_rep)
end

function read_cres_fields(paths::AbstractVector{<:AbstractString}, fields::Tuple{Vararg{String}})
    isempty(fields) && throw(ArgumentError("at least one cres field is required"))
    names = Tuple(Symbol.(fields))
    values_by_file = map(paths) do path
        matopen(path) do file
            cres = read(file, "liferes")["cres"]
            entries = if cres isa MAT.MatlabStructArray
                [cres[idx] for idx in 1:length(cres[first(fields)])]
            elseif cres isa AbstractDict
                [cres]
            else
                vec(cres)
            end
            map(fields) do field
                map(entries) do entry
                    value = entry[field]
                    value = value isa AbstractArray ? only(value) : value
                    value isa Real ||
                        throw(ArgumentError("$path: $field must be a numeric scalar"))
                    Float64(value)
                end
            end |> Tuple
        end
    end
    values = map(eachindex(fields)) do idx
        vcat((values_file[idx] for values_file in values_by_file)...)
    end |> Tuple
    NamedTuple{names}(values)
end

function read_atomnums(paths::AbstractVector{<:AbstractString})
    nums = read_cres_fields(paths, ("atomnum",)).atomnum
    all(isfinite, nums) || throw(ArgumentError("atomnum must contain only finite values"))
    nums
end

function dualmot_curve_style(condition::NamedTuple)
    loadcfg, istp = condition.loadcfg, condition.istp
    hue = HUE_ISTP[istp]
    strokecolor = RGB(Oklch(LIGHTNESS_STROKE, CHROMA_STROKE, hue))
    markercolor = RGB(Oklch(LIGHTNESS_FACE, CHROMA_FACE, hue))
    color = strokecolor
    linewidth, strokewidth, markersize = 0.75, 0.75, 5
    marker = MARKER_LOADCFG[loadcfg]
    line_options = (; color, linewidth)
    marker_options = (; color=markercolor, markersize, strokecolor, strokewidth)
    errorbar_options = (; errorlinewidth=0.75)
    (; color, markercolor, linewidth, strokecolor, strokewidth, markersize, marker,
        line_options, marker_options, errorbar_options)
end

dualmot_ratio_style(istp::Symbol; numerator::Symbol=:DCS) =
    dualmot_curve_style((; loadcfg=numerator, istp))

function dualmot_axis_kwargs(; log_y::Bool=false, text_size::Real=8,
    compact_spacing::Bool=false)
    common = (
        xgridvisible=false,
        ygridvisible=false,
        xminorgridvisible=false,
        yminorgridvisible=false,
        xtickalign=1,
        ytickalign=1,
        xminortickalign=1,
        yminortickalign=1,
        xticksmirrored=true,
        yticksmirrored=true,
        xlabelsize=text_size,
        ylabelsize=text_size,
        xticklabelsize=text_size,
        yticklabelsize=text_size,
        titlesize=text_size,
        spinewidth=0.75,
        xtickwidth=0.75,
        ytickwidth=0.75,
        xminortickwidth=0.75,
        yminortickwidth=0.75,
    )
    compact_spacing && (common = merge(common, (
        xlabelpadding=2, ylabelpadding=3,
        xticklabelpad=1, yticklabelpad=2,
    )))
    log_y ? merge(common, (
        yticks=DUALMOT_LOG_MAJOR_TICKS,
        yminorticks=DUALMOT_LOG_MINOR_TICKS,
        yminorticksvisible=true,
    )) : merge(common, (
        yminorticks=IntervalsBetween(5),
        yminorticksvisible=true,
    ))
end

const DUALMOT_LEGEND_OPTIONS = (
    labelsize=22 / 3,
    patchsize=(6, 6),
    nbanks=2,
    rowgap=1,
    colgap=3,
    margin=(2, 2, 2, 2),
    padding=(2, 2, 2, 2),
    framevisible=false,
    backgroundcolor=:transparent,
)

function set_time_minor_ticks!(ax::Axis)
    reset_limits!(ax)
    n_major_ticks = length(ax.xaxis.tickvalues[])
    ax.xminorticks = IntervalsBetween(n_major_ticks < 4 ? 5 : 2)
    ax.xminorticksvisible = true
    nothing
end

balance_tag(bias::Real) = iszero(bias) ? "t" : "n"

function dualmot_num_evol_plot_spec(kind::AbstractString;
    key_x::Symbol,
    xlabel::Union{AbstractString,Function},
    file_head::AbstractString,
    scale_x::Real=1.0,
    transform_x=(values, condition, panel, idx_istp) -> values,
    axis_options=(;),
    allow_partial_rep::Bool=false,
    xautolimits=(condition, panel) -> true,
    yautolimits=(condition, panel) -> true,
    legend_position::Symbol=:rt,
    loadcfg_plot::Tuple=(:DDM, :DIS),
    ylabel::AbstractString="CMOT number",
    limits_linear=nothing,
)
    isfinite(scale_x) && scale_x > 0 ||
        throw(ArgumentError("scale_x must be finite and positive, got $scale_x"))
    xlabel_panel = xlabel isa AbstractString ? (_ -> xlabel) : xlabel
    (
        key_x,
        scale_x=Float64(scale_x),
        key_panel=:β_MOT,
        curves=(loadcfg=loadcfg_plot, istp=:all),
        scales=(:lin, :log),
        formats=("svg", "png"),
        size=(275, 205),
        frame_size=(227, 151),
        fontsize=8,
        fit_size=(570, 270),
        fit_fontsize=8,
        scale_num_linear=1e7,
        xlabel=xlabel_panel,
        transform_x,
        axis_options,
        allow_partial_rep,
        xautolimits,
        yautolimits,
        ylabel,
        limits_linear,
        title=(tag, bias, reps_used) ->
            "$tag · reps = $reps_used\nβ_MOT = $bias",
        filename=(scale, bias) ->
            "[$file_head].[$scale].[$(balance_tag(bias))-balanced]",
        curve_label=condition -> "$(condition.istp) $(condition.loadcfg)",
        curve_style=dualmot_curve_style,
        legend_position,
    )
end

dualmot_lifetime_plot_spec(kind::AbstractString; key_x::Symbol,
    xlabel::AbstractString, scale_x::Real=1.0, ylabel::AbstractString="CMOT number") =
    dualmot_num_evol_plot_spec(kind; key_x, xlabel, scale_x, ylabel,
        file_head="$kind.lifetime")
using LsqFit: curve_fit, stderror
using Printf: @sprintf

"""One- and two-body loss, with p = [N₀ (atoms), τ (s), κ (atom⁻¹ s⁻¹)]."""
function model_num_decay(t::AbstractVector, p::AbstractVector)
    n0, tau, kappa = p
    @. n0 * exp(-t / tau) / (1 + kappa * n0 * tau * (-expm1(-t / tau)))
end

"""Two-body-only loss (τ → ∞), with p = [N₀ (atoms), κ (atom⁻¹ s⁻¹)]."""
function model_num_decay_kappa(t::AbstractVector, p::AbstractVector)
    n0, kappa = p
    @. n0 / (1 + kappa * n0 * t)
end

"""One-body-only loss (κ = 0), with p = [N₀ (atoms), τ (s)]."""
function model_num_decay_tau(t::AbstractVector, p::AbstractVector)
    n0, tau = p
    @. n0 * exp(-t / tau)
end

"""Bounded, unweighted least squares; uncertainties are local residual-based 1σ errors."""
function fit_num_decay(t::AbstractVector{<:Real}, nums::AbstractVector{<:Real};
    mode::Symbol=:full, bounds_n0::Tuple{<:Real,<:Real}=(0.5, 2.0),
    bounds_tau::Tuple{<:Real,<:Real}=(0.1, 50.0),
    bounds_kappa::Tuple{<:Real,<:Real}=(eps(Float64), Inf),
)
    length(t) == length(nums) || throw(DimensionMismatch("time and number lengths differ"))
    mode in (:full, :kappa, :tau) ||
        throw(ArgumentError("decay fit mode must be :full, :kappa, or :tau, got $mode"))
    0 < bounds_n0[1] < bounds_n0[2] ||
        throw(ArgumentError("N₀ multiplier bounds must satisfy 0 < lower < upper, got $bounds_n0"))
    mode in (:full, :tau) && !(0 < bounds_tau[1] < bounds_tau[2]) &&
        throw(ArgumentError("τ bounds must satisfy 0 < lower < upper, got $bounds_tau"))
    mode in (:full, :kappa) && !(0 <= bounds_kappa[1] < bounds_kappa[2]) &&
        throw(ArgumentError("κ bounds must satisfy 0 ≤ lower < upper, got $bounds_kappa"))
    mask = isfinite.(t) .& isfinite.(nums) .& (nums .>= 0)
    ts, ns = Float64.(t[mask]), Float64.(nums[mask])
    length(unique(ts)) > 3 || throw(ArgumentError("decay fitting needs at least four distinct valid times"))
    minimum(ts) >= 0 || throw(ArgumentError("holding times must be nonnegative"))
    n_initial = ns[argmin(ts)]
    n_initial > 0 || throw(ArgumentError("initial number must be positive, got $n_initial"))
    # Scale numbers and κ to avoid ill-conditioned finite differences.
    model, scale, lower, upper, starts = if mode == :full
        (model_num_decay,
            [n_initial, 1.0, inv(n_initial)],
            [bounds_n0[1], bounds_tau[1], bounds_kappa[1] * n_initial],
            [bounds_n0[2], bounds_tau[2], bounds_kappa[2] * n_initial],
            ([1.0, tau, loss] for tau in (0.3, 3.0, 30.0) for loss in (0.01, 1.0)))
    elseif mode == :kappa
        (model_num_decay_kappa,
            [n_initial, inv(n_initial)],
            [bounds_n0[1], bounds_kappa[1] * n_initial],
            [bounds_n0[2], bounds_kappa[2] * n_initial],
            ([1.0, loss] for loss in (0.001, 0.01, 0.1, 1.0, 10.0)))
    else
        (model_num_decay_tau,
            [n_initial, 1.0],
            [bounds_n0[1], bounds_tau[1]],
            [bounds_n0[2], bounds_tau[2]],
            ([1.0, tau] for tau in (0.3, 3.0, 30.0)))
    end
    model_scaled = (x, p) -> model(x, p .* scale) ./ n_initial
    fits = []
    for start in starts
        p0 = clamp.(start, lower, upper)
        fit = curve_fit(model_scaled, ts, ns ./ n_initial, p0; lower, upper, maxIter=1000)
        fit.converged && push!(fits, fit)
    end
    isempty(fits) && error("decay fit did not converge for any initial guess")
    fit = fits[argmin([sum(abs2, f.resid) for f in fits])]
    params = fit.param .* scale
    errors = try
        stderror(fit) .* scale
    catch err
        @warn "Decay parameter covariance unavailable" exception=err
        fill(NaN, length(params))
    end
    at_bound = any(isapprox.(fit.param, lower; atol=1e-7, rtol=1e-4) .|
        isapprox.(fit.param, upper; atol=1e-7, rtol=1e-4))
    (; mode, model, params, errors, at_bound, fit, mask)
end

function format_fit_scientific(value::Real, error::Real)
    isfinite(value) && !iszero(value) || return "($(value) ± $(error))"
    exponent = floor(Int, log10(abs(value)))
    scale = 10.0^exponent
    if abs(round(value / scale; digits=2)) >= 10
        exponent += 1
        scale *= 10
    end
    value_text = @sprintf("%.2f", value / scale)
    error_scaled = abs(error) / scale
    error_text = if !isfinite(error_scaled)
        string(error_scaled)
    elseif round(error_scaled; digits=2) == 0
        "0"
    else
        @sprintf("%.2f", error_scaled)
    end
    "($value_text ± $error_text)e$exponent"
end

function label_num_decay(result)
    p, e = result.params, result.errors
    suffix = result.at_bound ? " *" : ""
    label_n0 = format_fit_scientific(p[1], e[1])
    result.mode == :full && return @sprintf(
        "N₀ = %s\nτ = (%.3g ± %.2g) s\nκ = %s atom⁻¹ s⁻¹%s",
        label_n0, p[2], e[2], format_fit_scientific(p[3], e[3]), suffix)
    result.mode == :kappa && return @sprintf(
        "N₀ = %s\nκ = %s atom⁻¹ s⁻¹%s",
        label_n0, format_fit_scientific(p[2], e[2]), suffix)
    @sprintf("N₀ = %s\nτ = (%.3g ± %.2g) s%s",
        label_n0, p[2], e[2], suffix)
end

function save_num_decay_results(path_root::AbstractString, dataset::AbstractString,
    fit_mode::Symbol, records::AbstractVector)
    isempty(records) && throw(ArgumentError("cannot save an empty $dataset decay-fit result"))
    path_results = joinpath(path_root, "Fit results")
    mkpath(path_results)
    created_at = Dates.now()
    stamp = Dates.format(created_at, dateformat"yyyymmdd-HHMMSS") * "-" *
        lpad(string(Dates.millisecond(created_at)), 3, '0')
    path = joinpath(path_results, "[$dataset.decay.fit].[$stamp].jld2")
    JLD2.jldsave(path;
        schema_version=1,
        dataset,
        fit_mode,
        created_at=string(created_at),
        records=collect(records),
    )
    println("Saved $dataset decay fits to $path")
    path
end

function load_latest_num_decay_results(path_root::AbstractString, dataset::AbstractString)
    path_results = joinpath(path_root, "Fit results")
    isdir(path_results) ||
        throw(ArgumentError("$dataset decay-fit folder does not exist: $path_results"))
    prefix = "[$dataset.decay.fit]."
    paths = filter(readdir(path_results; join=true)) do path
        isfile(path) && startswith(basename(path), prefix) && endswith(path, ".jld2")
    end
    isempty(paths) && throw(ArgumentError("no $dataset decay-fit JLD2 files in $path_results"))
    path = paths[argmax(mtime.(paths))]
    payload = JLD2.load(path)
    get(payload, "schema_version", nothing) == 1 ||
        throw(ArgumentError("unsupported decay-fit schema in $path"))
    get(payload, "dataset", nothing) == dataset ||
        throw(ArgumentError("expected $dataset decay fits in $path"))
    records = get(payload, "records", nothing)
    records isa AbstractVector || throw(ArgumentError("missing decay-fit records in $path"))
    (; path, fit_mode=Symbol(payload["fit_mode"]), records)
end
