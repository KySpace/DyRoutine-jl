using YAML
using XLSX
using MAT
using Statistics: mean, std
using CairoMakie
using Colors: Oklch, RGB
using Dates
using JLD2
using Printf: @sprintf

isdefined(@__MODULE__, :marker_errorbars!) ||
    include(joinpath(@__DIR__, "..", "snippets", "marker_errorbars.jl"))

function log_plot_limits(values::AbstractVector{<:Real}, errors::AbstractVector{<:Real};
    cap_factor::Real=10, margin::Real=1.04)
    length(values) == length(errors) || throw(DimensionMismatch(
        "log plot values and errors must have equal lengths"))
    cap_factor > 1 || throw(ArgumentError("log plot cap factor must exceed 1"))
    margin >= 1 || throw(ArgumentError("log plot margin must be at least 1"))
    positive = Float64[value for value in values if isfinite(value) && value > 0]
    isempty(positive) && throw(ArgumentError("log plot needs at least one positive value"))
    min_value, max_value = extrema(positive)
    floor_cap = min_value / cap_factor
    ceil_cap = max_value * cap_factor
    lows = Float64[]
    highs = Float64[]
    for (value_raw, error_raw) in zip(values, errors)
        value = Float64(value_raw)
        isfinite(value) && value > 0 || continue
        error = Float64(error_raw)
        if isfinite(error) && error >= 0
            push!(lows, max(value - error, floor_cap))
            push!(highs, min(value + error, ceil_cap))
        else
            push!(lows, value)
            push!(highs, value)
        end
    end
    low = max(minimum(lows) / margin, floor_cap)
    high = min(maximum(highs) * margin, ceil_cap)
    if !(low < high)
        low, high = floor_cap, ceil_cap
    end
    low, high
end

function marker_errorbars_log!(ax::Axis, x::AbstractVector{<:Real},
    y::AbstractVector{<:Real}, error::AbstractVector{<:Real};
    floor::Real, kwargs...)
    length(x) == length(y) == length(error) || throw(DimensionMismatch(
        "log error-bar coordinate and error lengths must match"))
    floor > 0 || throw(ArgumentError("log error-bar floor must be positive"))
    all(value -> isfinite(value) && value > floor, y) || throw(ArgumentError(
        "log error-bar values must be finite and above the axis floor"))
    error_low = [isfinite(err) && err >= 0 ? min(Float64(err), Float64(value) - floor) : 0.0
        for (value, err) in zip(y, error)]
    error_high = [isfinite(err) && err >= 0 ? Float64(err) : 0.0 for err in error]
    marker_errorbars!(ax, x, y, error_low, error_high; kwargs...)
end

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
const DUALMOT_CONFIG_VAR_NAMES = Dict{Symbol,String}(
    key => spec.config for specs in (DUALMOT_VAR_SPECS,
        DUALMOT_LOADING_VAR_SPECS, DUALMOT_LOADCFG_VAR_SPECS,
        DUALMOT_BALANCE_VAR_SPECS, DUALMOT_ODT_BFIELD_VAR_SPECS)
    for (key, spec) in pairs(specs))

const ODT_BFIELD_CALIBRATION = Dict(
    :x => (Q=1.48094, B0=0.276317),
    :z => (Q=1.71727, B0=-0.351795),
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
        subscript(rich(string(direction), font=:italic)), " (G)")
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
const DUALMOT_LOG_MINOR_TICKS = sort!([Float64(multiplier) * 10.0^exponent
    for exponent in -20:20 for multiplier in 2:9])
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

    config_to_name = Dict{String,Symbol}()
    for (name, spec) in pairs(var_specs)
        config_to_name[spec.config] = name
        config_to_name[string(name)] = name
    end
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
    selector = parse_num_evol_selector(get(config, "selector", nothing);
        var_specs, label)
    (; date_runid, vars, var_order, files, sources, selector)
end

function parse_num_evol_selector(raw_selector, ; var_specs::NamedTuple, label::AbstractString)
    isnothing(raw_selector) && return NamedTuple()
    raw_selector isa AbstractVector ||
        throw(ArgumentError("$label: selector must be a sequence of variable mappings"))
    config_to_name = Dict{String,Symbol}()
    for (name, spec) in pairs(var_specs)
        config_to_name[spec.config] = name
        config_to_name[string(name)] = name
    end
    selected = Pair{Symbol,NamedTuple}[]
    for (idx_selector, item) in enumerate(raw_selector)
        item isa AbstractDict && !isempty(item) ||
            throw(ArgumentError("$label selector[$idx_selector]: expected a variable mapping"))
        raw_key = first(keys(item))
        raw_predicates = item[raw_key]
        if isnothing(raw_predicates) && length(item) == 2
            # YAML also accepts the compact sibling form used in existing configs:
            # - rep:\n  value: "a -> ..."
            remaining = filter(pair -> first(pair) != raw_key, collect(item))
            raw_predicates = Dict(remaining)
        end
        key = get(config_to_name, string(raw_key)) do
            raw_key == "rep" ? :rep : throw(ArgumentError("$label selector[$idx_selector]: unknown variable $raw_key"))
        end
        raw_predicates isa AbstractDict && length(raw_predicates) == 1 ||
            throw(ArgumentError("$label selector[$idx_selector]: expected exactly one of value or index"))
        raw_kind, raw_expression = first(raw_predicates)
        kind = Symbol(raw_kind)
        kind in (:value, :index) ||
            throw(ArgumentError("$label selector[$idx_selector]: selector key must be value or index"))
        raw_expression isa AbstractString ||
            throw(ArgumentError("$label selector[$idx_selector]: $kind predicate must be a Julia expression string"))
        predicate = try
            Core.eval(@__MODULE__, Meta.parse(raw_expression))
        catch error
            throw(ArgumentError("$label selector[$idx_selector]: cannot parse $kind predicate: $(sprint(showerror, error))"))
        end
        predicate isa Function ||
            throw(ArgumentError("$label selector[$idx_selector]: $kind predicate must evaluate to a function"))
        push!(selected, key => NamedTuple{(kind,)}((predicate,)))
    end
    length(unique(first.(selected))) == length(selected) ||
        throw(ArgumentError("$label: selector may mention each variable only once"))
    NamedTuple{Tuple(first.(selected))}(last.(selected))
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
        (; folder=String(folder), tag=String(tag),
            for_load421=get(run_config, "for_load421", false), data)
    end
    tags = getproperty.(runinfos, :tag)
    allunique(tags) || throw(ArgumentError("$folder: processing tags must be unique"))
    runinfos
end

function latest_num_evol_runinfo(runinfos::AbstractVector{<:NamedTuple};
    label::AbstractString="processing entries",
)
    isempty(runinfos) && throw(ArgumentError("$label: no candidates"))
    source_key(runinfo) = maximum(data_key for data in runinfo.data
        for data_key in data.date_runid)
    keys_source = source_key.(runinfos)
    idx_latest = findmax(keys_source)[2]
    count(==(keys_source[idx_latest]), keys_source) == 1 ||
        throw(ArgumentError("$label: latest source date/run is not unique"))
    runinfos[idx_latest]
end

function select_for_load421_runinfo(runinfos::AbstractVector{<:NamedTuple},
    pair::AbstractString)
    length(runinfos) == 1 && return only(runinfos)
    selected = filter(runinfo -> runinfo.for_load421 === true, runinfos)
    length(selected) == 1 || throw(ArgumentError(
        "$pair MOT loading 421 comparison requires exactly one tagged group with for_load421: true when multiple groups exist; found $(length(selected)) among $(getproperty.(runinfos, :tag))"))
    only(selected)
end

filename_token(value::AbstractString) = replace(strip(value),
    r"[\[\]<>:\"/\\|?*\s]+" => "-")

function figure_xlsx_cell(value)
    ismissing(value) && return "NaN"
    value isa AbstractFloat && isnan(value) && return "NaN"
    value isa AbstractFloat && isinf(value) && return string(value)
    value
end

function figure_xlsx_source(runinfo::NamedTuple, condition::NamedTuple)
    sources = String[]
    for data in runinfo.data
        all(key -> !hasproperty(data.vars, key) ||
            getproperty(condition, key) in getproperty(data.vars, key), keys(condition)) || continue
        append!(sources, replace.(basename.(data.sources),
            r"^liferes\s+" => "", r"\.mat$" => ""))
    end
    join(unique(sources), ", ")
end

function figure_xlsx_sheet_name(value::AbstractString)
    name = filename_token(value)
    isempty(name) && (name = "data")
    first(name, min(length(name), 31))
end

function num_evol_xlsx_sheet_name(runinfo::NamedTuple, key_panel, panel;
    fit::Bool=false)
    isnothing(key_panel) && return figure_xlsx_sheet_name(runinfo.tag *
        (fit ? " fit" : ""))
    suffix = "__β$(@sprintf("%+.2f", Float64(panel)))" * (fit ? "_fit" : "")
    max_prefix = max(1, 31 - length(suffix))
    tag = figure_xlsx_sheet_name(runinfo.tag)
    tag = tag[1:min(lastindex(tag), max_prefix)]
    tag * suffix
end

function write_figure_workbook(path::AbstractString, sheets)
    isempty(sheets) && return nothing
    mkpath(dirname(path))
        names_written = Set{String}()
        XLSX.openxlsx(path; mode="w") do workbook
            for (idx_sheet, (raw_name, raw_matrix)) in enumerate(sheets)
                name_base = figure_xlsx_sheet_name(raw_name)
                name = name_base
                suffix = 2
                while name in names_written
                    tail = "-$suffix"
                    name = first(name_base, min(length(name_base), 31 - length(tail))) * tail
                    suffix += 1
                end
                push!(names_written, name)
                sheet = if idx_sheet == 1
                    first_sheet = workbook[1]
                    XLSX.renamesheet!(first_sheet, name)
                    first_sheet
                else
                    XLSX.addsheet!(workbook, name)
                end
                matrix = map(figure_xlsx_cell, raw_matrix)
                sheet["A1"] = matrix
            end
        end
    nothing
end

function write_pair_figure_workbooks(path_root::AbstractString,
    sheets_by_pair::AbstractDict, filename::AbstractString)
    for (pair, sheets) in sheets_by_pair
        isempty(sheets) && continue
        path_pair = joinpath(path_root, pair)
        mkpath(path_pair)
        write_figure_workbook(joinpath(path_pair, filename), sheets)
    end
    nothing
end

function num_evol_figure_matrix(runinfo::NamedTuple, key_x::Symbol, val_x,
    key_panel, panel, curves_plot::AbstractVector)
    headers = String["$(get(DUALMOT_CONFIG_VAR_NAMES, key_x, string(key_x))) (raw)"]
    append!(headers, reduce(vcat, map(curves_plot) do curve
        condition = join(string.(values(curve.condition)), " ")
        ["$(key_x) plot [$condition]", "mean [$condition]", "std [$condition]",
            "n_rep [$condition]", "sources [$condition]"]
    end; init=String[]))
    output = Matrix{Any}(undef, length(curves_plot[1].nums) + 1, length(headers))
    output[1, :] .= headers
    for idx_x in eachindex(curves_plot[1].nums)
        output[idx_x + 1, 1] = val_x[idx_x]
    end
    for (idx_curve, curve) in enumerate(curves_plot)
        col = 2 + 5 * (idx_curve - 1)
        condition_axes = merge(curve.condition,
            isnothing(key_panel) ? NamedTuple() : NamedTuple{(key_panel,)}((panel,)))
        for idx_x in eachindex(curve.nums)
            idx_row = idx_x + 1
            output[idx_row, col] = curve.val_x_curve[idx_x]
            output[idx_row, col + 1] = curve.nums[idx_x]
            output[idx_row, col + 2] = curve.stds[idx_x]
            output[idx_row, col + 3] = curve.n_reps[idx_x] == 0 ? NaN : curve.n_reps[idx_x]
            output[idx_row, col + 4] = figure_xlsx_source(runinfo,
                merge(condition_axes, NamedTuple{(key_x,)}((val_x[idx_x],))))
        end
    end
    output
end

function calc_num_evol_block(runinfo_data::NamedTuple;
    label::AbstractString,
    bounds_sigmax_num::Tuple{<:Real,<:Real},
    bounds_sigmay_num::Tuple{<:Real,<:Real},
    num_max_num::Real,
    field_num::AbstractString="atomnum",
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

    isempty(field_num) && throw(ArgumentError("$label: field_num must be nonempty"))
    shot_data = read_cres_fields(runinfo_data.files, (field_num, "sigmax", "sigmay"))
    num_data, sigmax_data, sigmay_data = shot_data
    all(value -> ismissing(value) || isfinite(value), num_data) ||
        throw(ArgumentError("$label: $field_num must contain only finite values or baddata entries"))
    len_data = length(num_data)
    n_variation = prod(length(getproperty(runinfo_data.vars, key))
        for key in keys(runinfo_data.vars) if key != :rep)
    len_data > 0 || throw(DimensionMismatch("$label: $field_num data must not be empty"))
    n_partial = rem(len_data, n_variation)
    iszero(n_partial) || allow_partial_rep ||
        throw(DimensionMismatch("$label: $len_data $field_num values must be a positive integer multiple of $n_variation variations"))
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
    mask_num_low = .!ismissing.(num_data) .& (coalesce.(num_data, -Inf) .>= - num_max_num / 100)
    mask_num_high = .!ismissing.(num_data) .& (coalesce.(num_data, Inf) .<= num_max_num)
    mask_valid = mask_size .& mask_num_low .& mask_num_high
    num_data_masked = Vector{Union{Missing, Float64}}(num_data)
    num_data_masked[.!mask_valid] .= missing
    append!(num_data_masked, fill(missing, n_padding))
    sigmax_masked = Vector{Union{Missing, Float64}}(sigmax_data)
    sigmax_masked[.!mask_valid] .= missing
    append!(sigmax_masked, fill(missing, n_padding))
    sigmay_masked = Vector{Union{Missing, Float64}}(sigmay_data)
    sigmay_masked[.!mask_valid] .= missing
    append!(sigmay_masked, fill(missing, n_padding))

    fmt_from_flat(data) = reshape(data, reverse(n_dims_acq)) |>
        values -> permutedims(values, reverse(1:length(n_dims_acq))) |>
        values -> permutedims(values, indexin(collect(name), collect(name_acq)))
    num_fmt = fmt_from_flat(num_data_masked)
    sigmax_fmt = fmt_from_flat(sigmax_masked)
    sigmay_fmt = fmt_from_flat(sigmay_masked)
    for key in keys(runinfo_data.selector)
        axis_values = getproperty(vars, key)
        predicate_spec = getproperty(runinfo_data.selector, key)
        kind = first(keys(predicate_spec))
        predicate = getproperty(predicate_spec, kind)
        selector_values = kind == :value ? axis_values : eachindex(axis_values)
        selected = predicate.(selector_values)
        selected isa AbstractVector{Bool} && length(selected) == length(axis_values) ||
            throw(ArgumentError("$label selector for $key.$kind must return one Boolean per axis value"))
        idx_axis = findfirst(==(key), name)
        for fmt in (num_fmt, sigmax_fmt, sigmay_fmt), idx in eachindex(selected)
            selected[idx] || (selectdim(fmt, idx_axis, idx) .= missing)
        end
    end
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
        sigmax_fmt, sigmay_fmt, mask_valid, n_masked_total)
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
        for field in (:num_fmt, :sigmax_fmt, :sigmay_fmt)
            size(getproperty(block, field)) == expected_size ||
                throw(DimensionMismatch("$label data[$idx_block]: $field size $(size(getproperty(block, field))) must match variable dimensions $expected_size"))
        end
        for key in name
            values_axis = getproperty(block.vars, key)
            allunique(values_axis) ||
                throw(ArgumentError("$label data[$idx_block]: duplicate values in $key: $(collect(values_axis))"))
        end
    end

    n_rep = sum(block.n_rep for block in blocks)
    values_combined = map(name) do key
        key == :rep && return 1:n_rep
        values = unique(vcat((collect(getproperty(block.vars, key)) for block in blocks)...))
        key in (:t_load, :t_hold, :ib) && sort!(values)
        values
    end
    vars = NamedTuple{name}(values_combined)
    size_combined = Tuple(length.(values_combined))
    combine_fmt(field) = begin
        output = Array{Union{Missing, Float64}}(undef, size_combined)
        fill!(output, missing)
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
            @views output[indices...] .= getproperty(block, field)
            pos_rep += block.n_rep
        end
        output
    end
    num_fmt = combine_fmt(:num_fmt)
    sigmax_fmt = combine_fmt(:sigmax_fmt)
    sigmay_fmt = combine_fmt(:sigmay_fmt)
    (; vars, num_fmt, sigmax_fmt, sigmay_fmt, n_rep)
end

function read_cres_fields(paths::AbstractVector{<:AbstractString}, fields::Tuple{Vararg{String}})
    isempty(fields) && throw(ArgumentError("at least one cres field is required"))
    names = Tuple(Symbol.(fields))
    values_by_file = map(paths) do path
        matopen(path) do file
            liferes = read(file, "liferes")
            cres = liferes["cres"]
            entries = if cres isa MAT.MatlabStructArray
                [cres[idx] for idx in 1:length(cres[first(fields)])]
            elseif cres isa AbstractDict
                [cres]
            else
                vec(cres)
            end
            field_values = map(fields) do field
                map(entries) do entry
                    value = entry[field]
                    value = value isa AbstractArray ? only(value) : value
                    value isa Real ||
                        throw(ArgumentError("$path: $field must be a numeric scalar"))
                    Float64(value)
                end
            end |> Tuple
            baddata = get(liferes, "baddata", Int[])
            bad_indices = baddata isa AbstractArray ? vec(baddata) : [baddata]
            # MATLAB liferes files use scalar zero to mean that no bad shots were recorded.
            bad_indices = length(bad_indices) == 1 && iszero(only(bad_indices)) ? eltype(bad_indices)[] : bad_indices
            all(index -> index isa Real && isinteger(index) && 1 <= index <= length(entries), bad_indices) ||
                throw(ArgumentError("$path: baddata must contain linear indices in 1:$(length(entries))"))
            for field_num in ("atomnum", "pixsum")
                field_num in fields || continue
                idx_num = findfirst(==(field_num), fields)
                nums = Union{Missing,Float64}[field_values[idx_num]...]
                nums[Int.(bad_indices)] .= missing
                field_values = ntuple(idx -> idx == idx_num ? nums : field_values[idx], length(fields))
            end
            field_values
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

function dualmot_title(tag, reps_used; bias=nothing)
    title_reps = rich(tag, " · ", rich("n", font=:italic),
        subscript("reps"), " = ", string(reps_used))
    isnothing(bias) && return title_reps
    rich(title_reps, " · ", rich("β", font=:italic), subscript("MOT"),
        " = ", string(bias))
end

function dualmot_num_evol_plot_spec(kind::AbstractString;
    key_x::Symbol,
    xlabel::Union{AbstractString,Function},
    file_head::AbstractString,
    scale_x::Real=1.0,
    transform_x=(values, condition, panel, idx_istp) -> values,
    axis_options=(;),
    draw_background=(ax, panel, scale) -> nothing,
    allow_partial_rep::Bool=false,
    field_num::AbstractString="atomnum",
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
        draw_background,
        allow_partial_rep,
        field_num=String(field_num),
        xautolimits,
        yautolimits,
        ylabel,
        limits_linear,
        title=(tag, bias, reps_used) -> dualmot_title(tag, reps_used; bias),
        filename=(scale, bias) ->
            "[$file_head].[$scale].[$(balance_tag(bias))-balanced]",
        curve_label=condition -> "$(condition.istp) $(condition.loadcfg)",
        curve_style=dualmot_curve_style,
        legend_position,
    )
end

dualmot_lifetime_plot_spec(kind::AbstractString; key_x::Symbol,
    xlabel::AbstractString, scale_x::Real=1.0, ylabel::AbstractString="CMOT number",
    field_num::AbstractString="atomnum") =
    dualmot_num_evol_plot_spec(kind; key_x, xlabel, scale_x, ylabel, field_num,
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

"""One-body-only loss (κ = 0), with p = [N₀, τ (s)]."""
function model_num_decay_tau(t::AbstractVector, p::AbstractVector)
    n0, tau = p
    @. n0 * exp(-t / tau)
end

"""Bounded, unweighted least squares; uncertainties are local residual-based 1σ errors."""
function fit_num_decay(t::AbstractVector{<:Real}, nums::AbstractVector{<:Real};
    mode::Symbol=:full, bounds_n0::Tuple{<:Real,<:Real}=(0.5, 2.0),
    bounds_tau::Tuple{<:Real,<:Real}=(0.1, 50.0),
    bounds_kappa::Tuple{<:Real,<:Real}=(eps(Float64), Inf),
    selector::Function=values -> trues(length(values)),
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
    mask_selected = selector(t)
    mask_selected isa AbstractVector{Bool} || throw(ArgumentError(
        "decay fit selector must return an AbstractVector{Bool}, got $(typeof(mask_selected))",
    ))
    length(mask_selected) == length(t) || throw(DimensionMismatch(
        "decay fit selector returned $(length(mask_selected)) values for $(length(t)) times",
    ))
    mask = isfinite.(t) .& isfinite.(nums) .& (nums .>= 0) .& mask_selected
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
    (; mode, model, params, errors, at_bound, fit, mask,
        mask_selected=collect(mask_selected))
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

function label_num_decay(result; number_unit::AbstractString="atom")
    p, e = result.params, result.errors
    suffix = result.at_bound ? " *" : ""
    label_n0 = format_fit_scientific(p[1], e[1])
    result.mode == :full && return @sprintf(
        "N₀ = %s\nτ = (%.3g ± %.2g) s\nκ = %s %s⁻¹ s⁻¹%s",
        label_n0, p[2], e[2], format_fit_scientific(p[3], e[3]), number_unit, suffix)
    result.mode == :kappa && return @sprintf(
        "N₀ = %s\nκ = %s %s⁻¹ s⁻¹%s",
        label_n0, format_fit_scientific(p[2], e[2]), number_unit, suffix)
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
    expected_mode = dataset == "CMOT" ? :kappa : dataset == "MOT" ? :tau : nothing
    candidates = sort(paths; by=mtime, rev=true)
    selected = findfirst(candidates) do candidate
        payload = JLD2.load(candidate)
        get(payload, "schema_version", nothing) == 1 ||
            throw(ArgumentError("unsupported decay-fit schema in $candidate"))
        get(payload, "dataset", nothing) == dataset ||
            throw(ArgumentError("expected $dataset decay fits in $candidate"))
        isnothing(expected_mode) || Symbol(payload["fit_mode"]) == expected_mode
    end
    isnothing(selected) && throw(ArgumentError(
        "no $dataset decay-fit result using expected mode $expected_mode in $path_results"))
    path = candidates[selected]
    payload = JLD2.load(path)
    records = get(payload, "records", nothing)
    records isa AbstractVector || throw(ArgumentError("missing decay-fit records in $path"))
    (; path, fit_mode=Symbol(payload["fit_mode"]), records)
end
