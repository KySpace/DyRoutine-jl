using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))

year_test = 2026
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
title_anlz = "DevTest"

date, runid = "0325", 80
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
smw_ft = 10
px_in_um = 6.5 / 22.06

name = ["repeat", "t_hold", "istp"]
val = (
    collect(1:3),
    collect(6:2:200),
    ["162", "164"],
)
n_variation = length(val[1]) * length(val[2]) * length(val[3])
n_dim_vars = map(length, val);
n_rep, n_main, n_istp = n_dim_vars
h5open(path_input, "r") do f
    global dens = f["/od"] |>
                  read |>
                  x -> permutedims(x, (3, 2, 1)) |>
                       x -> stack(
                      map(d -> subtract_corner_mean(d, wh_corner), eachslice(x; dims=1));
                      dims=1,
                  )
    ndims(dens) == 3 || error("Expected /od to have 3 dimensions, got $(ndims(dens)).")
end

_, h_dens, w_dens = size(dens)
wh_dens = (w_dens, h_dens)
dens_mean = dropdims(mean(dens; dims=1); dims=1)
xy_peak_px = find_positive_cluster_center(dens_mean; smwh=smwh_peak) |> cent -> round.(Int, cent)
dens_full_fmt = dens |>
                ds -> mapslices(d -> crop_center(d, xy_peak_px, smwh_peak), ds; dims=(2, 3)) |>
                      ds -> reshape(ds, (reverse(n_dim_vars)..., reverse(wh_peak)...)) |>
                            ds -> permutedims(ds, (3, 2, 1, 4, 5))

# A lite version for tests
rng_lite = 1:50;
val = (
    collect(1:3),
    collect(6:2:200)[rng_lite],
    ["162", "164"],
)
n_variation = length(val[1]) * length(val[2]) * length(val[3])
n_dim_vars = map(length, val);
n_rep, n_main, n_istp = n_dim_vars
dens_full_fmt = dens_full_fmt[:,rng_lite,:,:,:]

# Statistics on number sum
# num_fmt = dens_full_fmt |> ds -> mapslices(calc_dens_sum, ds; dims=(4, 5)) |> n -> dropdims(n; dims=(4, 5));
# stat_n_fmt = num_fmt |> a -> mapslices(calc_mean_std, a; dims=(1))

# fig_num, ax_num = set_axis!("number vs t hold")
# for (i, istp) in enumerate(val[3])
#     plot_num_stat_evo!(ax_num, val[2], stat_n_fmt[1, :, i], val[3][i])
# end
# display(fig_num)

fig_live = Figure()
gl = GridLayout()
fig_live[1, 1] = gl
axs_live = set_panel_solo_modl!(gl);
fig_live |> display

essn_2d_fmt = dens_full_fmt |> ds -> mapslices(d -> calc_solo_essn_2d(d, smwh_peak .+ 1, smwh_peak, smw_ft, px_in_um), ds; dims=(4, 5)) |> e -> dropdims(e; dims=(4, 5));
info_fmt = [Dict("istp" => val[3][i], "t_hold" => val[2][t], "repeat" => val[1][r]) for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]]
essn_stacked_over_rep = [
    calc_stacked_essn(@view essn_2d_fmt[:, t, i])
    for t in axes(essn_2d_fmt, 2), i in axes(essn_2d_fmt, 3)
]
essn_stacked_over_rep_t = [
    calc_stacked_essn((@view essn_2d_fmt[:, :, i]) |> vec)
    for i in axes(essn_2d_fmt, 3)
]

modl2d_side = essn_2d_fmt |> f -> map(a -> a.modl2d, f) |>
                                  m ->
    map(a -> a[smwh_peak[2]+1+8:smwh_peak[2]+1+15, smwh_peak[1]+1-smw_ft:smwh_peak[1]+1+smw_ft], m) |>
    # d -> [permutedims(stack(@view d[i, j, :]), (3, 1, 2))
    #     for i in axes(d, 1), j in axes(d, 2)];
    d -> d
modes_pca_modl2d = [modl2d_side[:, :, i] |> m -> fit_pca_modes(8, m) for i in 1:n_istp]

fig_pca, axs_pca = set_axis_pca_dual_4x2!()
for idx_mode in 1:8, istp in 1:n_istp
    plot_mode_evol_freq_solo!(axs_pca[istp, idx_mode], modes_pca_modl2d[istp][idx_mode], val[2])
end
resize_to_layout!(fig_pca)
display(fig_pca)
# fig_full, axs_solo, axs_stacked = set_axis_full(n_dim_vars, set_panel_solo_essn_2d!)
fig_full, axs_solo, axs_stacked = set_axis_full(n_dim_vars, set_panel_solo_modl!)
for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
    info = info_fmt[r, t, i]
    print("\rplotting for rep $r, $(info["t_hold"]) ms, $(info["istp"])")
    draw_solo_modl!(axs_solo[r, t, i], essn_2d_fmt[r, t, i], info)
    draw_solo_modl!(axs_live, essn_2d_fmt[r, t, i], info)
end
for t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
    info = info_fmt[1, t, i] |> d -> merge(d, Dict("repeat" => "stacked"))
    print("\rplotting for stacked $(info["t_hold"]) ms, $(info["istp"])")
    draw_solo_modl!(axs_stacked[t, i], essn_stacked_over_rep[t, i], info)
    draw_solo_modl!(axs_live, essn_stacked_over_rep[t, i], info)
end
resize_to_layout!(fig_full)

fig_full |> f -> save(joinpath(path_output, "full_essn_CFNM_5.318_rastr.pdf"), f; backend=CairoMakie)
