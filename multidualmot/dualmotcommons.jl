using YAML
using MAT
using Statistics: mean, std
using CairoMakie
using Colors: Oklch, RGB

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

const HUE_ISTP = Dict(Symbol(string(i)) => h for (i, h) in
    ((160, 195), (161, 306), (162, 21), (163, 90), (164, 259)))
const LIGHTNESS_STROKE, CHROMA_STROKE = 0.45, 0.10
const LIGHTNESS_FACE, CHROMA_FACE = 0.85, 0.06
const LIGHTNESS_DIS_LINE, CHROMA_DIS_LINE = 0.65, 0.08
const MARKER_LOADCFG = Dict(
    :DDM => :circle, :DIS => :utriangle, :DCS => :rect, :SCS => :circle)

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
    stroke = RGB(Oklch(LIGHTNESS_STROKE, CHROMA_STROKE, hue))
    face = RGB(Oklch(LIGHTNESS_FACE, CHROMA_FACE, hue))
    line_color = loadcfg in (:DIS, :SCS) ?
        RGB(Oklch(LIGHTNESS_DIS_LINE, CHROMA_DIS_LINE, hue)) : stroke
    (; color=line_color, markercolor=face, linewidth=2, strokecolor=stroke,
        strokewidth=1.5, markersize=11, marker=MARKER_LOADCFG[loadcfg])
end

balance_tag(bias::Real) = iszero(bias) ? "t" : "n"

function dualmot_num_evol_plot_spec(kind::AbstractString;
    key_x::Symbol,
    xlabel::AbstractString,
    file_head::AbstractString,
    scale_x::Real=1.0,
    legend_position::Symbol=:rt,
    loadcfg_plot::Tuple=(:DDM, :DIS),
    ylabel::AbstractString="CMOT number",
)
    isfinite(scale_x) && scale_x > 0 ||
        throw(ArgumentError("scale_x must be finite and positive, got $scale_x"))
    (
        key_x,
        scale_x=Float64(scale_x),
        key_panel=:β_MOT,
        curves=(loadcfg=loadcfg_plot, istp=:all),
        scales=(:lin, :log),
        formats=("svg", "png"),
        size=(500, 360),
        scale_num_linear=1e7,
        xlabel,
        ylabel,
        title=(tag, bias, reps_used) ->
            "$tag · β_MOT = $bias · $(balance_tag(bias))-balanced · reps = $reps_used",
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
