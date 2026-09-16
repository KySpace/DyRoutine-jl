include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT lifetime"
val_pair = ["162-164", "162-160", "161-162", "161-164", "163-164", "162-163 n-balanced", "162-163 t-balanced", "163-161"]
var_specs_num = DUALMOT_VAR_SPECS
key_x_num = :t_hold
size_min_num = 4e-4

runinfos_grouped = [read_num_evol_runinfo(path_root, pair;
    var_specs=var_specs_num,
    validate_vars=validate_dualmot_vars,
) for pair in val_pair]
runinfos = runinfos_grouped
ids_runinfo = eachindex(runinfos)
plot_num_evol = dualmot_lifetime_plot_spec("MOT";
    key_x=key_x_num,
    xlabel="t_hold (ms)",
)

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = runinfo.tag_head
    global path_output = joinpath(path_root, tag_head)
    println("Processing: $tag_head")
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
end
