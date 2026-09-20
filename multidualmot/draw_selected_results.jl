include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_data = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data"
path_result = joinpath(path_data, "Result")
pair = "162-164"
γ_active = 0.43
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8

fig_result = Figure(size=(325, 300), fontsize=8, figure_padding=1)
colgap!(fig_result.layout, 0)
rowgap!(fig_result.layout, 0)
axis_options_ticks = (xticksize=10 / 3, yticksize=10 / 3,
    xminorticksize=2, yminorticksize=2)
axis_options_numbers = merge(axis_options_ticks,
    (yticks=0:5:10, yminorticks=IntervalsBetween(5)))
axis_options_loading = merge(axis_options_numbers,
    (xticks=0:5:20, xminorticks=IntervalsBetween(2),
        xminorticksvisible=true))
limits_loading = (x=(-0.5, 18.0), y=(0.0, 11.5))

# 421 loading, zero-bias panel.
path_root = joinpath(path_data, "MOT loading 421")
runinfos = read_num_evol_runinfos(path_root, pair;
    var_specs=DUALMOT_LOADING_VAR_SPECS,
    validate_vars=validate_dualmot_vars)
runinfo = only(filter(info -> length(info.data) == 1 &&
    0.0 in only(info.data).vars.β_MOT, runinfos))
tag_head = "$(runinfo.folder) $(runinfo.tag)"
plot_num_evol = dualmot_num_evol_plot_spec("MOT";
    key_x=:t_load,
    xlabel="Effective loading time (s)",
    file_head="MOT.loading.421",
    transform_x=(values, condition, panel, idx_istp) -> values .* γ_active,
    legend_position=:rb,
    loadcfg_plot=(:DDM, :DIS, :DCS))
plot_num_evol_target = (fig=fig_result, slot=fig_result[1, 1],
    panel=0.0, scale=:lin, axis_options=axis_options_loading,
    limits=limits_loading, show_legend=false)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

# 626 loading, number curves only.
path_root = joinpath(path_data, "MOT loading 626")
runinfo = only(read_num_evol_runinfos(path_root, pair;
    var_specs=DUALMOT_LOADCFG_VAR_SPECS,
    validate_vars=validate_loadcfg_comparison_vars))
tag_head = "$(runinfo.folder) $(runinfo.tag)"
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
plot_cmpr_loadcfg_target = (fig=fig_result, slot=fig_result[1, 2],
    axis_options=axis_options_loading, limits=limits_loading)
include(joinpath(@__DIR__, "anlz_cmpr_loadcfg.jl"))

# CMOT lifetime, zero-bias panel.
path_root = joinpath(path_data, "CMOT lifetime")
runinfo = only(filter(info -> length(info.data) == 1 &&
    0.0 in only(info.data).vars.β_MOT,
    read_num_evol_runinfos(path_root, pair;
        var_specs=DUALMOT_VAR_SPECS,
        validate_vars=validate_dualmot_vars)))
tag_head = "$(runinfo.folder) $(runinfo.tag)"
plot_num_evol = dualmot_lifetime_plot_spec("CMOT";
    key_x=:t_hold,
    scale_x=1000,
    xlabel="CMOT holding time (s)",
    ylabel="CMOT number")
plot_num_evol_target = (fig=fig_result, slot=fig_result[2, 1],
    panel=0.0, scale=:lin,
    axis_options=merge(axis_options_numbers,
        (xminorticks=IntervalsBetween(2), xminorticksvisible=true)),
    limits=(x=nothing, y=(0.0, 11.5)), show_legend=false)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

# ODT B-field z component.
path_root = joinpath(path_data, "ODT BField")
runinfo = only(filter(info -> info.tag == "z",
    read_num_evol_runinfos(path_root, pair;
        var_specs=DUALMOT_ODT_BFIELD_VAR_SPECS,
        validate_vars=validate_odt_bfield_vars)))
tag_head = "$(runinfo.folder) $(runinfo.tag)"
bounds_sigmax_num = (2e-4, 10e-4)
bounds_sigmay_num = (1e-4, 5e-4)
plot_num_evol = merge(dualmot_num_evol_plot_spec("ODT";
    key_x=:ib,
    xlabel="B field current in z (A)",
    ylabel="ODT number",
    file_head="ODT.BField",
    loadcfg_plot=(:DDM, :DIS)),
    (key_panel=nothing,))
plot_num_evol_target = (fig=fig_result, slot=fig_result[2, 2],
    panel=nothing, scale=:lin, axis_options=axis_options_ticks,
    limits=nothing, show_legend=false)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

mkpath(path_result)
for format in ("svg", "png")
    save(joinpath(path_result, "selected_162-164.$format"), fig_result)
end
println("Saved selected 162-164 panels to $path_result")
