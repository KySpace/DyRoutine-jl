using TOML
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

const HUE_ISTP = Dict(Symbol(string(i)) => h for (i, h) in
    ((160, 195), (161, 306), (162, 21), (163, 90), (164, 259)))
const LIGHTNESS_STROKE, CHROMA_STROKE = 0.45, 0.10
const LIGHTNESS_FACE, CHROMA_FACE = 0.85, 0.06
const LIGHTNESS_DIS_LINE, CHROMA_DIS_LINE = 0.65, 0.08
const MARKER_LOADCFG = Dict(:DDM => :circle, :DIS => :utriangle)

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

function read_num_evol_runinfo(path_root::AbstractString, pair::AbstractString;
    var_specs::NamedTuple,
    excluded_date_runid=Tuple{String, Int}[],
    validate_vars=(vars, pair) -> nothing,
)
    path_pair = joinpath(path_root, pair)
    config = TOML.parsefile(joinpath(path_pair, "config.toml"))
    names_var = keys(var_specs)
    values_var = map(values(var_specs)) do spec
        haskey(config, spec.config) || throw(ArgumentError("$pair: missing configuration $(spec.config)"))
        raw_values = config[spec.config]
        raw_values isa AbstractVector && !isempty(raw_values) ||
            throw(ArgumentError("$pair: $(spec.config) must be a nonempty vector"))
        spec.convert(raw_values)
    end
    vars = merge((rep=:auto,), NamedTuple{names_var}(values_var))
    validate_vars(vars, pair)

    config_to_name = Dict(spec.config => name for (name, spec) in pairs(var_specs))
    raw_order = get(config, "vars", nothing)
    raw_order isa AbstractVector || throw(ArgumentError("$pair: vars must be a vector"))
    var_order = map(raw_order) do raw_name
        raw_name == "rep" ? :rep : get(config_to_name, string(raw_name)) do
            throw(ArgumentError("$pair: unknown variable in vars: $raw_name"))
        end
    end |> Tuple
    expected_order = (:rep, names_var...)
    Set(var_order) == Set(expected_order) && length(var_order) == length(expected_order) ||
        throw(ArgumentError("$pair: vars must order $(join(expected_order, ", ")) exactly once"))

    files = filter(readdir(path_pair)) do filename
        occursin(r"^liferes \d{4} run\d+\.mat$", filename)
    end
    isempty(files) && throw(ArgumentError("$pair: no renamed liferes MAT files found"))
    date_runid = map(files) do filename
        matched = match(r"^liferes (\d{4}) run(\d+)\.mat$", filename)
        (String(matched[1]), parse(Int, matched[2]))
    end
    allunique(date_runid) || throw(ArgumentError("$pair: duplicate date/run IDs"))
    keep = findall(id -> id ∉ excluded_date_runid, date_runid)
    files, date_runid = files[keep], date_runid[keep]
    isempty(files) && throw(ArgumentError("$pair: no liferes files remain after exclusions"))
    order = sortperm(date_runid)
    (; tag_head=String(pair), date_runid=date_runid[order], vars, var_order,
        files=joinpath.(path_pair, files[order]))
end

read_dualmot_runinfo(path_root::AbstractString, pair::AbstractString; kwargs...) =
    read_num_evol_runinfo(path_root, pair; var_specs=DUALMOT_VAR_SPECS,
        validate_vars=validate_dualmot_vars, kwargs...)

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
    line_color = loadcfg == :DIS ?
        RGB(Oklch(LIGHTNESS_DIS_LINE, CHROMA_DIS_LINE, hue)) : stroke
    (; color=line_color, markercolor=face, linewidth=2, strokecolor=stroke,
        strokewidth=1.5, markersize=11, marker=MARKER_LOADCFG[loadcfg])
end

balance_tag(bias::Real) = iszero(bias) ? "t" : "n"

function dualmot_lifetime_plot_spec(kind::AbstractString;
    key_x::Symbol,
    xlabel::AbstractString,
)
    file_head = "$kind.lifetime"
    (
        key_x,
        key_panel=:β_MOT,
        curves=(loadcfg=(:DDM, :DIS), istp=:all),
        scales=(:lin, :log),
        formats=("svg", "png"),
        size=(500, 360),
        scale_num_linear=1e7,
        xlabel,
        ylabel="$kind atom number",
        title=(tag, bias, reps_used) ->
            "$tag · β_MOT = $bias · $(balance_tag(bias))-balanced · reps = $reps_used",
        filename=(scale, bias) ->
            "[$file_head].[$scale].[$(balance_tag(bias))-balanced]",
        curve_label=condition -> "$(condition.istp) $(condition.loadcfg)",
        curve_style=dualmot_curve_style,
    )
end
