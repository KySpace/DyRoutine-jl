include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_result = joinpath(path_root, "Result")
pair = "162-164"
γ_active = 0.43
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8

# At 96 dpi, this gives a 178 mm wide exported figure.
fig_width = round(Int, 178 / 25.4 * 96)
fig_height = 240
fig_result_loading = Figure(size=(fig_width, fig_height), fontsize=8,
    figure_padding=4)
fig_result_decay = Figure(size=(fig_width, fig_height), fontsize=8,
    figure_padding=4)
axis_options_ticks = (xticksize=10 / 3, yticksize=10 / 3,
    xminorticksize=2, yminorticksize=2,
    xlabelsize=7 / 0.75, ylabelsize=7 / 0.75)
axis_options_numbers = merge(axis_options_ticks,
    (yticks=0:5:10, yminorticks=IntervalsBetween(5)))
axis_options_loading = merge(axis_options_numbers,
    (xticks=0:5:20,
        xminorticks=setdiff(collect(0:1:20), collect(0:5:20)),
        xminorticksvisible=true))
limits_loading = (x=nothing, y=(0.0, 11.5))
limits_decay = (x=nothing, y=(2e6, 1e8))
axis_options_decay = merge(axis_options_ticks,
    (yticks=DUALMOT_LOG_MAJOR_TICKS,
        yminorticks=DUALMOT_LOG_MINOR_TICKS, yminorticksvisible=true))
mask_loading_421 = values -> values .!= 40.0
mask_loading_626 = values -> values .<= 30.0
panel_target(fig, col; panel=0.0, scale=:lin, axis_options,
    limits=limits_loading, mask_x=values -> trues(length(values))) =
    (fig=fig, slot=fig[1, col], panel, scale, axis_options, limits,
        show_legend=false, mask_x, compact_spacing=false)

# MOT loading 421, zero-bias panel. Use the alternative 0712 dataset.
path_dataset = joinpath(path_root, "MOT loading 421")
runinfos = read_num_evol_runinfos(path_dataset, pair;
    var_specs=DUALMOT_LOADING_VAR_SPECS,
    validate_vars=validate_dualmot_vars)
runinfos_421_0712 = filter(info -> info.tag == "0712" &&
    any(data -> 0.0 in data.vars.β_MOT, info.data), runinfos)
length(runinfos_421_0712) == 1 || throw(ArgumentError(
    "$pair selected MOT loading 421 requires one zero-bias 0712 dataset; found $(length(runinfos_421_0712))"))
runinfo = only(runinfos_421_0712)
tag_head = "$(runinfo.folder) $(runinfo.tag)"
plot_num_evol = dualmot_num_evol_plot_spec("MOT";
    key_x=:t_load,
    xlabel="Effective loading time (s)",
    file_head="MOT.loading.421",
    transform_x=(values, condition, panel, idx_istp) -> values .* γ_active,
    legend_position=:rb,
    loadcfg_plot=(:DDM, :DIS, :DCS),
    ylabel=rich("CMOT ", rich("N", subscript("i"); font=:italic)))
plot_num_evol_target = panel_target(fig_result_loading, 1;
    axis_options=axis_options_loading, limits=limits_loading,
    mask_x=mask_loading_421)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

# MOT loading 626, number curves only; omit t_load > 30 s.
path_dataset = joinpath(path_root, "MOT loading 626")
runinfo = latest_num_evol_runinfo(read_num_evol_runinfos(path_dataset, pair;
    var_specs=DUALMOT_LOADCFG_VAR_SPECS,
    validate_vars=validate_loadcfg_comparison_vars);
    label="$pair selected MOT loading 626")
tag_head = "$(runinfo.folder) $(runinfo.tag)"
plot_cmpr_loadcfg = (
    file_head="MOT.loading.626",
    xlabel="Effective loading time (s)",
    ylabel_num=rich("CMOT ", rich("N", subscript("i"); font=:italic)),
    γ_active,
    scale_num=1e7,
    formats=("svg", "png"),
    size=(275, 205),
    frame_size=(227, 151),
    fontsize=8,
)
plot_cmpr_loadcfg_target = (fig=fig_result_loading,
    slot=fig_result_loading[1, 2], axis_options=axis_options_loading,
    limits=limits_loading, mask_x=mask_loading_626,
    show_legend=false, compact_spacing=false)
include(joinpath(@__DIR__, "anlz_cmpr_loadcfg.jl"))

# MOT loading balance.
using CSV
path_dataset = joinpath(path_root, "MOT loading balance")
runinfo = latest_num_evol_runinfo(read_num_evol_runinfos(path_dataset, pair;
    var_specs=DUALMOT_BALANCE_VAR_SPECS,
    validate_vars=validate_balance_vars);
    label="$pair selected MOT loading balance")
tag_head = "$(runinfo.folder) $(runinfo.tag)"
rows_balance = collect(CSV.File(joinpath(path_dataset, "balance.csv");
    header=["Pair", "tbiasmot"], skipto=2, stripwhitespace=true))
balance_bias = only(Float64(row.tbiasmot) for row in rows_balance
    if strip(string(row.Pair)) == pair)
values_beta = reduce(vcat, [data.vars.β_MOT for data in runinfo.data])
step_major = 0.5
idx_min = floor(Int, minimum(values_beta) / step_major)
idx_max = ceil(Int, maximum(values_beta) / step_major)
axis_options_balance = merge(axis_options_ticks,
    (xticks=collect(idx_min:idx_max) .* step_major,
        xminorticks=IntervalsBetween(5), xminorticksvisible=true))
plot_num_evol = merge(dualmot_num_evol_plot_spec("MOT";
    key_x=:β_MOT,
    xlabel=_ -> rich("β", subscript("MOT")),
    ylabel=rich("CMOT ", rich("N", subscript("i"); font=:italic)),
    file_head="MOT.loading.balance",
    loadcfg_plot=(:DDM, :DIS)),
    (key_panel=nothing,
        draw_background=(ax, _, _) -> vlines!(ax, [balance_bias];
            color=RGBAf(0.75, 0.75, 0.75, 1), linewidth=0.75)))
plot_num_evol_target = panel_target(fig_result_loading, 3;
    panel=nothing, axis_options=axis_options_balance, limits=nothing)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
# Add zero-valued DDM endpoints to the statistics so their curve segments,
# markers, and zero-sized error bars use the ordinary plotting path.
idx_beta_axis = findfirst(==(:β_MOT), name_stat)
idx_loadcfg_axis = findfirst(==(:loadcfg), name_stat)
idx_istp_axis = findfirst(==(:istp), name_stat)
beta_values_old = vars.β_MOT
beta_values = sort!(unique(vcat(Float64.(beta_values_old), [-1.0, 1.0])))
idx_beta_new = [findfirst(==(value), beta_values) for value in beta_values_old]
function insert_beta_axis(values::AbstractArray, fill_value)
    dims = ntuple(idx -> idx == idx_beta_axis ? length(beta_values) : size(values, idx),
        ndims(values))
    result = fill(convert(eltype(values), fill_value), dims)
    for (idx_old, idx_new) in enumerate(idx_beta_new)
        src = ntuple(idx -> idx == idx_beta_axis ? idx_old : Colon(), ndims(values))
        dst = ntuple(idx -> idx == idx_beta_axis ? idx_new : Colon(), ndims(values))
        result[dst...] .= values[src...]
    end
    result
end
num_stat = insert_beta_axis(num_stat, NaN)
std_num_stat = insert_beta_axis(std_num_stat, NaN)
n_rep_stat = insert_beta_axis(n_rep_stat, 0)
vars = merge(vars, (; β_MOT=beta_values))
val_vars = values(vars)
val_x = beta_values
val_x_plot = val_x ./ plot_num_evol.scale_x
for (istp, β_MOT) in ((Symbol("164"), -1.0), (Symbol("162"), 1.0))
    indices = Any[Colon() for _ in name_stat]
    indices[idx_beta_axis] = findfirst(==(β_MOT), beta_values)
    indices[idx_loadcfg_axis] = findfirst(==(:DDM), vars.loadcfg)
    indices[idx_istp_axis] = findfirst(==(istp), vars.istp)
    num_stat[indices...] = 0.0
    std_num_stat[indices...] = 0.0
    n_rep_stat[indices...] = 0
end
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

# CMOT decay, zero-bias panel.
path_dataset = joinpath(path_root, "CMOT lifetime")
runinfo = only(filter(info -> length(info.data) == 1 &&
    0.0 in only(info.data).vars.β_MOT,
    read_num_evol_runinfos(path_dataset, pair;
        var_specs=DUALMOT_VAR_SPECS,
        validate_vars=validate_dualmot_vars)))
tag_head = "$(runinfo.folder) $(runinfo.tag)"
plot_num_evol = dualmot_lifetime_plot_spec("CMOT";
    key_x=:t_hold,
    scale_x=1000,
    xlabel="CMOT holding time (s)",
    ylabel=rich("CMOT ", rich("N", subscript("i"); font=:italic)))
axis_options_cmot_decay = merge(axis_options_decay,
    (xticks=0:0.5:2, xminorticks=IntervalsBetween(4),
        xminorticksvisible=true))
plot_num_evol_target = panel_target(fig_result_decay, 2;
    scale=:log, axis_options=axis_options_cmot_decay, limits=limits_decay)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

# MOT decay, zero-bias panel; retain times through 50 s (including 50 s).
path_dataset = joinpath(path_root, "MOT lifetime")
runinfo = only(filter(info -> length(info.data) == 1 &&
    0.0 in only(info.data).vars.β_MOT,
    read_num_evol_runinfos(path_dataset, pair;
        var_specs=DUALMOT_VAR_SPECS,
        validate_vars=validate_dualmot_vars)))
tag_head = "$(runinfo.folder) $(runinfo.tag)"
plot_num_evol = dualmot_lifetime_plot_spec("MOT";
    key_x=:t_hold,
    scale_x=1000,
    xlabel="MOT holding time (s)",
    ylabel=rich("CMOT ", rich("N", subscript("i"); font=:italic)))
axis_options_mot_decay = merge(axis_options_decay,
    (xticks=0:20:60, xminorticks=IntervalsBetween(4),
        xminorticksvisible=true))
plot_num_evol_target = panel_target(fig_result_decay, 1;
    scale=:log, axis_options=axis_options_mot_decay, limits=limits_decay,
    mask_x=values -> values .<= 50.0)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

# ODT B-field z component, displayed in gauss.
path_dataset = joinpath(path_root, "ODT BField")
runinfo = only(filter(info -> info.tag == "z",
    read_num_evol_runinfos(path_dataset, pair;
        var_specs=DUALMOT_ODT_BFIELD_VAR_SPECS,
        validate_vars=validate_odt_bfield_vars)))
tag_head = "$(runinfo.folder) $(runinfo.tag)"
bounds_sigmax_num = (2e-4, 10e-4)
bounds_sigmay_num = (1e-4, 5e-4)
plot_num_evol = merge(dualmot_num_evol_plot_spec("ODT";
    key_x=:ib,
    xlabel=_ -> odt_bfield_xlabel(:z),
    transform_x=(values, condition, panel, idx_istp) ->
        odt_bfield_values(values, :z),
    axis_options=odt_bfield_axis_options(
        reduce(vcat, [data.vars.ib for data in runinfo.data]), :z),
    ylabel=rich("ODT ", rich("N", subscript("i"); font=:italic)),
    file_head="ODT.BField",
    loadcfg_plot=(:DDM, :DIS)),
    (key_panel=nothing,))
plot_num_evol_target = panel_target(fig_result_decay, 3;
    panel=nothing, axis_options=axis_options_ticks, limits=nothing)
include(joinpath(@__DIR__, "anlz_num_evol.jl"))
include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))

mkpath(path_result)
for (name, fig) in (("selected_162-164_loading", fig_result_loading),
    ("selected_162-164_decay", fig_result_decay)), format in ("svg", "png", "pdf")
    save(joinpath(path_result, "$name.$format"), fig)
end
println("Saved selected 162-164 figures to $path_result")
