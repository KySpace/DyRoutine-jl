include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_dataset = joinpath(path_root, "MOT loading 421 monofreq")
val_pair = ["162-163", "161-163"]
var_specs_num = DUALMOT_LOADING_VAR_SPECS
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8
γ_active = 0.43

function validate_monofreq_vars(vars::NamedTuple, folder::AbstractString)
    pair = Symbol.(split(folder, "-"))
    length(pair) == 2 || throw(ArgumentError("$folder: expected an isotope-pair folder"))
    selected_pair = filter(istp -> istp in vars.istp, pair)
    vars.istp == selected_pair || throw(ArgumentError(
        "$folder: configured isotopes $(collect(vars.istp)) must follow pair order $(collect(selected_pair))"))
    Set(vars.loadcfg) == Set((:SDS, :SIS)) || throw(ArgumentError(
        "$folder: loadcfg must contain SDS and SIS"))
    for key in (:β_MOT, :t_load, :loadcfg, :istp)
        values = getproperty(vars, key)
        allunique(values) || throw(ArgumentError("$folder: duplicate values in $key"))
    end
    nothing
end

runinfos_grouped = [read_num_evol_runinfos(path_dataset, pair;
    var_specs=var_specs_num,
    validate_vars=validate_monofreq_vars,
) for pair in val_pair]
runinfos = vcat(runinfos_grouped...)
ids_runinfo = eachindex(runinfos)
figure_data_sheets = Dict{String,Vector{Pair{String,Matrix{Any}}}}()
plot_num_evol_base = dualmot_num_evol_plot_spec("MOT";
    key_x=:t_load,
    scale_x=1,
    xlabel=bias -> iszero(bias) ?
        "Effective loading time (s)" : "Equivalent loading time (s)",
    transform_x=(values, condition, bias, idx_istp) ->
        iszero(bias) ? values .* γ_active : values,
    ylabel="CMOT number",
    file_head="MOT.loading.421.monofreq",
    legend_position=:rb,
    loadcfg_plot=(:SDS, :SIS),
)
plot_num_evol_target = nothing

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_dataset, runinfo.folder)
    local filename_run = (scale, bias) ->
        "[MOT.loading.421.monofreq].[$scale].[$(balance_tag(bias))-balanced.$(filename_token(runinfo.tag))]"
    global plot_num_evol = merge(plot_num_evol_base, (; filename=filename_run))
    println("Processing: $tag_head")
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
    include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))
end
write_pair_figure_workbooks(path_dataset, figure_data_sheets, "MOT loading 421 monofreq.xlsx")
