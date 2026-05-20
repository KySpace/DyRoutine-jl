# using HDF5
# using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
# using GLMakie
# using JLD2
# using Printf
# GLMakie.activate!()
# include(joinpath(@__DIR__, "..", "src", "helper.jl"))
# include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
# include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
# include(joinpath(@__DIR__, "..", "src", "percond.jl"))
# include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
# include(joinpath(@__DIR__, "..", "src", "corr.jl"))

# year_test = 2026
# path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
# rep = 1:3
# t_hold = 6:2:200
# istp = ["162", "164"]
# vars = (; rep, t_hold, istp)
# runinfos = [
#     (date="0325", runid=95, IB=5.311, tag_head="CFNM"),
#     (date="0325", runid=82, IB=5.313, tag_head="CFNM"),
#     (date="0325", runid=52, IB=5.316, tag_head="CFNM"),
#     (date="0325", runid=80, IB=5.318, tag_head="CFNM"),
#     (date="0325", runid=96, IB=5.321, tag_head="CFNM"),
#     (date="0325", runid=67, IB=5.322, tag_head="CFNM"),
#     (date="0325", runid=68, IB=5.328, tag_head="CFNM"),
#     (date="0325", runid=50, IB=5.328, tag_head="CFNM"),
#     (date="0325", runid=81, IB=5.332, tag_head="CFNM"),
#     (date="0325", runid=51, IB=5.333, tag_head="CFNM"),
#     (date="0325", runid=79, IB=5.336, tag_head="CFNM"),
#     (date="0325", runid=53, IB=5.338, tag_head="CFNM"),
# ]


# title_anlz = "[05.12].37.Correlations"
# runinfo = merge(runinfos[10], (; vars))
date, runid, tag = runinfo.date, runinfo.runid, gen_run_tag(runinfo)
println("Processing: $tag")
dir_test = gen_date_path(date, year_test)
file_data = gen_h5name(date, runid)
path_input = joinpath(path_root, dir_test, @sprintf("run%02d", runid), file_data)
path_output = joinpath(path_root, dir_test, "AnlzRoutine", title_anlz);
if !isdir(path_output)
    mkpath(path_output)
end

wh_corner = (10, 10)
smwh_peak = (30, 60)
wh_peak = smwh_peak .* 2 .+ 1
smw_peak, smh_peak = smwh_peak
smw_ft = 5
px_in_um = 6.5 / 22.06

name = propertynames(runinfo.vars)
r = format_dens_runinfo(runinfo; path_root, year_test, wh_corner, smwh_roi=smwh_peak)
val = r.val
n_dim_vars = r.n_dim_vars
n_variation = prod(n_dim_vars)
n_rep, n_main, n_istp = n_dim_vars
wh_dens = r.wh_dens
xy_peak_px = r.xy_peak_px
dens_full_fmt = r.dens_full_fmt

# A lite version for tests
# rng_lite = 1:50;
# val = (
#     rep=collect(1:3),
#     t_hold=collect(6:2:200)[rng_lite],
#     istp=["162", "164"],
# )
# n_dim_vars = map(length, val);
# n_variation = prod(n_dim_vars)
# n_rep, n_main, n_istp = n_dim_vars
# dens_full_fmt = dens_full_fmt[:, rng_lite, :]

# Statistics on number sum
# num_fmt = sum.(dens_full_fmt)
# stat_n_fmt = num_fmt |> a -> mapslices(calc_mean_std, a; dims=(1))

# fig_num, ax_num = set_axis!("number vs t hold")
# for (i, istp) in enumerate(val.istp)
#     plot_num_stat_evo!(ax_num, val.t_hold, stat_n_fmt[1, :, i], val.istp[i])
# end
# display(fig_num)

step_posi = px_in_um
step_modl = 1 / (2 * smwh_peak[2] * px_in_um)
x_vec, y_vec = smwh_peak |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi
x_modl, y_modl = (x_vec, y_vec) .* step_modl
essn_2d_fmt = map(d -> calc_solo_essn_2d(d, smwh_peak .+ 1, smwh_peak, smw_ft, px_in_um), dens_full_fmt)
info_fmt = [Dict("istp" => val.istp[i], "t_hold" => val.t_hold[t], "repeat" => val.rep[r]) for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]]
essn_stacked_over_rep = [
    calc_stacked_essn(@view essn_2d_fmt[:, t, i])
    for t in axes(essn_2d_fmt, 2), i in axes(essn_2d_fmt, 3)
]
essn_stacked_over_rep_t = [
    calc_stacked_essn((@view essn_2d_fmt[:, :, i]) |> vec)
    for i in axes(essn_2d_fmt, 3)
]
fit_prfl_modl_over_rep_t_1d = [
    essn_stacked_over_rep_t[istp] |>
    e -> fit_prfl_modl_twinpeak_decay_1d(y_modl, e.prfl_modl_norm_px, (y_modl .> 0.02))
    for istp in axes(essn_stacked_over_rep_t, 1)
]
extr_fmt = [
    essn_2d_fmt[r, t, i] |> e -> calc_solo_extr(e, fit_prfl_modl_over_rep_t_1d[i])
    for r in axes(essn_2d_fmt, 1), t in axes(essn_2d_fmt, 2), i in axes(essn_2d_fmt, 3)
]
extr_stacked_over_rep = [
    essn_stacked_over_rep[t, i] |> e -> calc_solo_extr(e, fit_prfl_modl_over_rep_t_1d[i])
    for t in axes(essn_stacked_over_rep, 1), i in axes(essn_stacked_over_rep, 2)
]
modl2d_side = essn_2d_fmt |> f -> map(a -> a.modl2d, f) |>
                                  m ->
    map(a -> a[smwh_peak[2]+1+8:smwh_peak[2]+1+15, smwh_peak[1]+1-smw_ft:smwh_peak[1]+1+smw_ft], m) |>
    # d -> [permutedims(stack(@view d[i, j, :]), (3, 1, 2))
    #     for i in axes(d, 1), j in axes(d, 2)];
    d -> d
modes_pca_modl2d = [modl2d_side[:, :, i] |> m -> fit_pca_modes(8, m) for i in 1:n_istp]

selector_t_sidepeak = t -> 0 .< t .< 20
selector_t_envelope = t -> 0 .< t .< 20
trend_sidepeak_nvlp = [
    extr_fmt[r, :, i] |> e -> anlz_trend_from_extr(val.t_hold, e, 1:1:100; selector_t_sidepeak, selector_t_envelope)
    for r in axes(extr_fmt, 1), i in axes(extr_fmt, 3)
]

trend_stacked_over_rep = [
    extr_stacked_over_rep[:, i] |> e -> anlz_trend_from_extr(val.t_hold, e, 1:1:100; selector_t_sidepeak, selector_t_envelope)
    for i in axes(extr_fmt, 3)
]
##  saving data, still problematic

# @save joinpath(path_output, @sprintf("%s_data.jld2", tag))
# val
# essn_2d_fmt
# info_fmt
# essn_stacked_over_rep
# essn_stacked_over_rep_t
# fit_prfl_modl_over_rep_t_1d
# extr_fmt[1,1,1].fit_tailess
# extr_stacked_over_rep
# modl2d_side
# modes_pca_modl2d

## Overall plots
fig_trend, axs_trend = set_axis_sidepeak_nvlp!(n_dim_vars, set_panel_trend_sidepeak_nvlp!, runinfo)
fig_nvlp, axs_nvlp = set_axis_stack_all!(n_dim_vars, set_panel_trend_nvlp!, runinfo)
for i in 1:n_istp
    trend = trend_sidepeak_nvlp[:, i]
    trend_stacked = trend_stacked_over_rep[i]
    val_istp = val.istp[i]
    plot_trend_all!(axs_trend, trend, trend_stacked, val_istp)
    plot_trend_nvlp!(axs_nvlp, trend, trend_stacked, val_istp)
    resize_to_layout!(fig_trend)
    resize_to_layout!(fig_nvlp)
    for format in ["pdf", "png"]
        fig_trend |> f -> save(joinpath(path_output, @sprintf("%s_%s_trend.%s", tag, val_istp, format)), f; backend=CairoMakie)
        fig_nvlp |> f -> save(joinpath(path_output, @sprintf("%s_%s_trend_nvlp.%s", tag, val_istp, format)), f; backend=CairoMakie)
    end
end
# fig_trend |> display
##

fig_pca, axs_pca = set_axis_pca_dual_4x2!()
for idx_mode in 1:8, idx_istp in 1:n_istp
    plot_mode_evol_freq_solo!(axs_pca[idx_istp, idx_mode], modes_pca_modl2d[idx_istp][idx_mode], val.t_hold)
end
resize_to_layout!(fig_pca)
fig_pca |> f -> save(joinpath(path_output, @sprintf("%s_pca.pdf", tag)), f; backend=CairoMakie)
# display(fig_pca)


## Large file generation for all shots

# fig_full, axs_solo, axs_stacked = set_axis_full(n_dim_vars, set_panel_solo_essn_2d!)
# fig_full, axs_solo, axs_stacked = set_axis_full(n_dim_vars, set_panel_solo_modl!)
# for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
#     info = info_fmt[r, t, i]
#     print("\r\033[2Kplotting for rep $r, $(info["t_hold"]) ms, $(info["istp"])")
#     draw_solo_modl!(axs_solo[r, t, i], extr_fmt[r, t, i], info)
#     # draw_solo_modl!(axs_live, extr_fmt[r, t, i], info)
# end
# println("Full axes ready: dimensions $(n_dim_vars)")
# for t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
#     info = info_fmt[1, t, i] |> d -> merge(d, Dict("repeat" => "stacked"))
#     print("\r\033[2Kplotting for stacked $(info["t_hold"]) ms, $(info["istp"])")
#     draw_solo_modl!(axs_stacked[t, i], extr_stacked_over_rep[t, i], info)
#     # draw_solo_modl!(axs_live, extr_stacked_over_rep[t, i], info)
# end
# println("Full modulation table drawn.")
# resize_to_layout!(fig_full)

# fig_full |> f -> save(joinpath(path_output, @sprintf("%s_essn_table.pdf", tag)), f; backend=CairoMakie)
# println("Full modulation plot saved.")
