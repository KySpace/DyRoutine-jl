include(joinpath(@__DIR__, "dualmotcommons.jl"))
using CSV

path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading balance"
val_pair = ["162-164", "160-162", "161-162", "161-164", "163-164",
    "162-163", "161-163"]
var_specs_num = DUALMOT_BALANCE_VAR_SPECS
key_x_num = :β_MOT
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8

path_balance = joinpath(path_root, "balance.csv")
rows_balance = collect(CSV.File(path_balance;
    header=["Pair", "tbiasmot"], skipto=2, stripwhitespace=true))
canonical_pair_name(pair) = begin
    isotopes = parse.(Int, split(strip(string(pair)), "-"))
    length(isotopes) == 2 ||
        throw(ArgumentError("invalid isotope pair in $path_balance: $pair"))
    join(sort(isotopes), "-")
end
balance_by_pair = Dict{String, Float64}()
for row in rows_balance
    pair = canonical_pair_name(row.Pair)
    haskey(balance_by_pair, pair) &&
        throw(ArgumentError("duplicate isotope pair in $path_balance: $pair"))
    bias = Float64(row.tbiasmot)
    isfinite(bias) || throw(ArgumentError("nonfinite tbiasmot for $pair in $path_balance"))
    balance_by_pair[pair] = bias
end
Set(keys(balance_by_pair)) == Set(val_pair) ||
    throw(ArgumentError("$path_balance must contain exactly: $(join(val_pair, ", "))"))

function balance_axis_options(values::AbstractVector{<:Real})
    step_major = 0.2
    idx_min = floor(Int, minimum(values) / step_major)
    idx_max = ceil(Int, maximum(values) / step_major)
    (
        xticks=collect(idx_min:idx_max) .* step_major,
        xminorticks=IntervalsBetween(4),
        xminorticksvisible=true,
    )
end

runinfos_grouped = [read_num_evol_runinfos(path_root, pair;
    var_specs=var_specs_num,
    validate_vars=validate_balance_vars,
) for pair in val_pair]
runinfos = vcat(runinfos_grouped...)
ids_runinfo = eachindex(runinfos)
plot_num_evol = merge(dualmot_num_evol_plot_spec("MOT";
    key_x=key_x_num,
    xlabel=_ -> rich("β", subscript("MOT")),
    ylabel="CMOT number",
    file_head="MOT.loading.balance",
    loadcfg_plot=(:DDM, :DIS),
    allow_partial_rep=true,
), (
    key_panel=nothing,
    scales=(:lin,),
    title=(tag, _, reps_used) -> dualmot_title(tag, reps_used),
    filename=(_, _) -> "[MOT.loading.balance].[lin].[balance]",
))
plot_num_evol_target = nothing

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_root, runinfo.folder)
    values_beta = reduce(vcat, [data.vars.β_MOT for data in runinfo.data])
    balance_bias = balance_by_pair[runinfo.folder]
    global plot_num_evol = merge(plot_num_evol, (
        axis_options=balance_axis_options(values_beta),
        draw_background=(ax, _, _) -> vlines!(ax, [balance_bias];
            color=RGBAf(0.75, 0.75, 0.75, 1), linewidth=0.75),
    ))
    println("Processing: $tag_head")
    GC.gc()
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
    include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))
end
