using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "pershot.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\2026-03\0325\run50\d0325r50.h5"
path_plot = joinpath(@__DIR__, "probe_temp_number_vs_t_hold.svg")
path_plot_peak = joinpath(@__DIR__, "probe_temp_avg_density_peak.svg")
path_plot_duet = joinpath(@__DIR__, "probe_temp_duet.svg")
path_plot_sheet = joinpath(@__DIR__, "probe_temp_contact_sheet.pdf")
wh_corner = (10, 10)
smwh_peak = (30, 60)
wh_peak = smwh_peak .* 2 .+ 1
smw_peak, smh_peak = smwh_peak
duet_color_max = 40.0
sheet_color_max = 40.0
sheet_gap = 6

name = ["repeat", "t_hold", "istp"]
val = (
    collect(1:3),
    collect(6:2:200),
    ["162", "164"],
)
n_variation = length(val[1]) * length(val[2]) * length(val[3])
n_dim_vars = map(length, val);
n_rep, n_main, n_istp = n_dim_vars
h5open(path, "r") do f
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

# Statistics on number sum
num_fmt = dens_full_fmt |> ds -> mapslices(calc_dens_sum, ds; dims=(4, 5)) |> n -> dropdims(n; dims=(4, 5));
stat_n_fmt = num_fmt |> a -> mapslices(calc_mean_std, a; dims=(1))

fig_num, ax_num = set_axis!("number vs t hold")
for (i, istp) in enumerate(val[3])
plot_num_stat_evo!(ax_num, val[2], stat_n_fmt[1,:,i], val[3][i])
end
display(fig_num)
