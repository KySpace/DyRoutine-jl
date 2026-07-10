using CairoMakie
using HDF5
using Printf
using Statistics
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

# after running partially from anlz_ssntfr_src.jl
fig_live = Figure()
axs_dens = Axis(fig_live[1, 1]; aspect=DataAspect())
axs_dens_mean = Axis(fig_live[1, 2]; aspect=DataAspect())
axs_strip = Axis(fig_live[2, 1:2]; aspect=DataAspect())
axs_prfl = Axis(fig_live[3, 1:2]; width=200)
rowsize!(fig_live.layout, 1, 200)
clrmap = [gen_clrmap_solo(hue_theme_istp[istp]; alpha_base=0.2, thres_alpha=0.1) for istp in ("162", "164")]

cfg = get_prfl_modl_1d_config(smwh_src)

ib, istp, rep = (5, 1, 95)
dens = dens_src_raw_fmt[ib, istp, rep]
dens_mean = dens_src_raw_fmt[ib, istp, :] |> mean
dens_core = crop_center(dens, xy_fixed_src, smwh_src) |> copy
dens_core_mean = crop_center(dens_mean, xy_fixed_src, smwh_src) |> copy

hm = heatmap!(axs_dens, x_dens, x_dens, dens_core'; colormap=clrmap[istp], colorrange=(0.0, maximum(vec(dens_core))))
heatmap!(axs_dens_mean, x_dens, x_dens, dens_core_mean'; colormap=clrmap[istp], colorrange=(0.0, maximum(vec(dens_core))))

smw, smh = smwh_src
tucky_prfl = tucky1d(smh; alpha=0.2)
idx_strip = (cfg.smh_dens_strip |> s -> (-s:1:s) .+ smh .+ 1)
idx_modl = (cfg.smw_modl |> s -> (-s:1:s) .+ smw .+ 1)
dens_strip = @view dens_core[idx_strip, :]

hm = heatmap!(axs_strip, x_dens, x_dens[idx_strip], dens_strip'; colormap=clrmap[istp], colorrange=(0.0, maximum(vec(dens_core))))
prfl = @pipe dens_core |> calc_prfl_modl_1d(_, smwh_src; step_modl)
lines!(axs_prfl, x_modl, prfl)
fig_live |> resize_to_layout!
fig_live |> display
