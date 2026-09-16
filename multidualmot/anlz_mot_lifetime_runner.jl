include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT lifetime"
val_pair = ["162-164", "160-162", "161-162", "161-164", "163-164", "162-163", "161-163"]
var_specs_num = DUALMOT_VAR_SPECS
key_x_num = :t_hold
size_min_num = 4e-4
num_max_num = 1e8

runinfos_grouped = [read_num_evol_runinfos(path_root, pair;
    var_specs=var_specs_num,
    validate_vars=validate_dualmot_vars,
) for pair in val_pair]
runinfos = vcat(runinfos_grouped...)
ids_runinfo = eachindex(runinfos)
plot_num_evol = dualmot_lifetime_plot_spec("MOT";
    key_x=key_x_num,
    scale_x=1000,
    xlabel="MOT holding time (s)",
    ylabel="CMOT number",
)

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_root, runinfo.folder)
    println("Processing: $tag_head")
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
end
