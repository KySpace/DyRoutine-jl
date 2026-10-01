include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_dataset = joinpath(path_root, "CMOT lifetime")
val_pair = ["162-164", "160-162", "161-162", "161-164", "163-164", "162-163", "161-163"]
var_specs_num = DUALMOT_VAR_SPECS
key_x_num = :t_hold
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8
fit_num_decay_config = (
    mode=:kappa,
    bounds_n0=(0.5, 2.0),
    bounds_kappa=(eps(Float64), Inf),
    selector=times -> times .<= 2.0,
)

runinfos_grouped = [read_num_evol_runinfos(path_dataset, pair;
    var_specs=var_specs_num,
    validate_vars=validate_dualmot_vars,
) for pair in val_pair]
runinfos = vcat(runinfos_grouped...)
ids_runinfo = eachindex(runinfos)
figure_data_sheets = Dict{String,Vector{Pair{String,Matrix{Any}}}}()
plot_num_evol = dualmot_lifetime_plot_spec("CMOT";
    field_num="atomnum",
    key_x=key_x_num,
    scale_x=1000,
    xlabel="CMOT holding time (s)",
    ylabel="CMOT number",
)
plot_num_evol_target = nothing
fit_records_num_decay = NamedTuple[]

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_dataset, runinfo.folder)
    println("Processing: $tag_head")
    GC.gc() # Release figures from preceding datasets before rendering more.
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
    include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))
    include(joinpath(@__DIR__, "anlz_num_evol_decay.jl"))
    include(joinpath(@__DIR__, "anlz_num_evol_size_output.jl"))
end
path_fit_results = save_num_decay_results(
    path_dataset, "CMOT", fit_num_decay_config.mode, fit_records_num_decay)
write_pair_figure_workbooks(path_dataset, figure_data_sheets, "CMOT lifetime.xlsx")
