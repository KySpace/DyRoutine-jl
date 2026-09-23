include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\ODT BField"
val_pair = ["162-164", "160-162", "161-162", "161-164", "163-164",
    "162-163", "161-163"]
var_specs_num = DUALMOT_ODT_BFIELD_VAR_SPECS
key_x_num = :ib
bounds_sigmax_num = (2e-4, 10e-4)
bounds_sigmay_num = (1e-4, 5e-4)
num_max_num = 2e8

runinfos_grouped = [read_num_evol_runinfos(path_root, pair;
    var_specs=var_specs_num,
    validate_vars=validate_odt_bfield_vars,
) for pair in val_pair]
runinfos = vcat(runinfos_grouped...)
ids_runinfo = eachindex(runinfos)
plot_num_evol_base = dualmot_num_evol_plot_spec("ODT";
    key_x=key_x_num,
    scale_x=1,
    xlabel="",
    ylabel="ODT number",
    file_head="ODT.BField",
    legend_position=:rt,
    loadcfg_plot=(:DDM, :DIS),
)
plot_num_evol_target = nothing

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    direction = Symbol(runinfo.tag)
    direction in (:x, :z) ||
        throw(ArgumentError("$(runinfo.folder): ODT BField tag must be x or z, got $(runinfo.tag)"))
    current_values = reduce(vcat, [data.vars.ib for data in runinfo.data])
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_root, runinfo.folder)
    global plot_num_evol = merge(plot_num_evol_base, (
        key_panel=nothing,
        xlabel=_ -> odt_bfield_xlabel(direction),
        transform_x=(values, condition, panel, idx_istp) ->
            odt_bfield_values(values, direction),
        axis_options=odt_bfield_axis_options(current_values, direction),
        title=(tag, _, reps_used) -> dualmot_title(tag, reps_used),
        filename=(scale, _) -> "[ODT.BField].[$scale].[$direction]",
    ))
    println("Processing: $tag_head")
    GC.gc()
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
    include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))
end
