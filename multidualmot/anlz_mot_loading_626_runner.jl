include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading 626"
val_pair = ["162-164", "160-162", "161-162", "161-164", "163-164", "162-163", "161-163"]
var_specs_num = DUALMOT_LOADCFG_VAR_SPECS
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8
γ_active = 0.43

runinfos_grouped = [read_num_evol_runinfos(path_root, pair;
    var_specs=var_specs_num,
    validate_vars=validate_loadcfg_comparison_vars,
) for pair in val_pair]
all(length(runinfos) == 1 for runinfos in runinfos_grouped) ||
    throw(ArgumentError("MOT loading 626 requires exactly one processing tag per pair folder"))
runinfos = only.(runinfos_grouped)
ids_runinfo = eachindex(runinfos)
plot_cmpr_loadcfg = (
    file_head="MOT.loading.626",
    xlabel="Effective loading time (s)",
    ylabel_num="CMOT number",
    γ_active,
    scale_num=1e7,
    formats=("svg", "png"),
    size=(275, 205),
    frame_size=(227, 151),
    fontsize=8,
)
plot_cmpr_loadcfg_target = nothing

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_root, runinfo.folder)
    println("Processing: $tag_head")
    include(joinpath(@__DIR__, "anlz_cmpr_loadcfg.jl"))
end
