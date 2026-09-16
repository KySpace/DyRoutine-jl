using TOML
using MAT
using Statistics: mean, std
using CairoMakie
using Colors: Oklch, RGB

path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\CMOT lifetime"
val_pair = ["162-164", "162-160", "161-162", "161-164", "163-164", "162-163", "163-161"]
exclude_date_runid = Dict("161-162" => [("0805", 43)])

function read_cmot_runinfo(path_root::AbstractString, pair::AbstractString;
    excluded_date_runid=Tuple{String, Int}[])
    path_pair = joinpath(path_root, pair)
    config = TOML.parsefile(joinpath(path_pair, "config.toml"))
    for key in ("tbiasmot", "t_hold", "loadcfg", "istp")
        haskey(config, key) || throw(ArgumentError("$pair: missing configuration $key"))
        config[key] isa AbstractVector && !isempty(config[key]) ||
            throw(ArgumentError("$pair: $key must be a nonempty vector"))
    end
    for key in ("tbiasmot", "t_hold")
        all(x -> x isa Real && !(x isa Bool) && isfinite(x), config[key]) ||
            throw(ArgumentError("$pair: $key must contain finite numbers"))
    end
    vars = (rep=:auto, β_MOT=Float64.(config["tbiasmot"]),
        t_hold=Float64.(config["t_hold"]),
        loadcfg=Symbol.(string.(config["loadcfg"])),
        istp=Symbol.(string.(config["istp"])))
    var_order = Symbol.(replace.(string.(config["vars"]), "tbiasmot" => "β_MOT"))
    Set(var_order) == Set((:rep, :β_MOT, :t_hold, :loadcfg, :istp)) && length(var_order) == 5 ||
        throw(ArgumentError("$pair: vars must order rep, tbiasmot, t_hold, loadcfg, and istp exactly once"))
    join(string.(vars.istp), "-") == pair ||
        throw(ArgumentError("$pair: istp values must match the pair in order"))
    length(vars.istp) == 2 || throw(ArgumentError("$pair: expected two isotopes"))
    all(in(vars.loadcfg), (:DDM, :DIS)) ||
        throw(ArgumentError("$pair: loadcfg must include DDM and DIS"))
    for key in (:β_MOT, :loadcfg, :istp)
        values = getproperty(vars, key)
        allunique(values) || throw(ArgumentError("$pair: duplicate values in $key"))
    end
    balances = [iszero(b) ? "t" : "n" for b in vars.β_MOT]
    allunique(balances) || throw(ArgumentError("$pair: multiple biases share an output filename"))
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

runinfos_grouped = [read_cmot_runinfo(path_root, pair) for pair in val_pair]
runinfos = runinfos_grouped
ids_runinfo = eachindex(runinfos)

hue_istp = Dict(Symbol(string(i)) => h for (i, h) in
    ((160, 195), (161, 306), (162, 21), (163, 90), (164, 259)))
lightness_stroke, chroma_stroke = 0.45, 0.10
lightness_face, chroma_face = 0.85, 0.06
lightness_dis_line, chroma_dis_line = 0.65, 0.08
marker_loadcfg = Dict(:DDM => :circle, :DIS => :utriangle)
size_figure = (500, 360)

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = runinfo.tag_head
    global path_output = joinpath(path_root, tag_head)
    global vars = runinfo.vars
    global name = keys(vars)
    global val_vars = values(vars)
    global val_tbiasmot = vars.β_MOT
    global val_t_hold = vars.t_hold
    global val_loadcfg = vars.loadcfg
    global val_istp = vars.istp
    println("Processing: $tag_head")
    include(joinpath(@__DIR__, "anlz_cmot_lifetime.jl"))
end
