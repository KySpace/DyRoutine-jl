using HDF5
using Statistics: mean
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save
include(joinpath(@__DIR__, "..", "src", "pershot.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\2026-03\0324\run64\d0324r64.h5"
path_plot = joinpath(@__DIR__, "probe_temp_number_vs_t_hold.svg")
path_plot_peak = joinpath(@__DIR__, "probe_temp_avg_density_peak.svg")
corner_height = 10
corner_width = 10
smwh_peak = (40, 120)

name = ["repeat", "t_hold", "istp"]
val = (
    collect(1:3),
    collect(6:2:200),
    [5, 0],
)
variation = length(val[1]) * length(val[2]) * length(val[3])

h5open(path, "r") do f
    global dens = f["/od"] |>
                  read |>
                  x -> permutedims(x, (3, 2, 1)) |>
                  x -> stack(
                      map(d -> subtract_corner_mean(d, corner_height, corner_width), eachslice(x; dims=1));
                      dims=1,
                  )
end

ndims(dens) == 3 || error("Expected /od to have 3 dimensions, got $(ndims(dens)).")
size(dens) == (variation, 401, 201) || error(
    "Expected permuted /od size to be ($variation, 401, 201), got $(size(dens)).",
)

_, height, width = size(dens)

stat_number = dens |> x -> summarize_repeat_number(x, val)
dens_by_variation = stat_number.dens_by_variation
val_number = stat_number.val_number
err_number = stat_number.err_number

write_number_plot(path_plot, val[2], val_number, err_number, val[3])

dens_mean = dropdims(mean(dens; dims=1); dims=1)
cx_peak, cy_peak = find_positive_cluster_center(dens_mean; smwh=smwh_peak)

smw_peak, smh_peak = smwh_peak
left_peak = cx_peak - smw_peak
right_peak = cx_peak + smw_peak
top_peak = cy_peak - smh_peak
bottom_peak = cy_peak + smh_peak

fig_peak = Figure(size=(820, 980))
ax_peak = Axis(
    fig_peak[1, 1];
    title="Average Density with Peak-Finding Window",
    xlabel="x",
    ylabel="y",
    yreversed=true,
    aspect=DataAspect(),
)
hm = heatmap!(ax_peak, 1:width, 1:height, dens_mean'; colormap=:viridis)
lines!(
    ax_peak,
    [left_peak, right_peak, right_peak, left_peak, left_peak],
    [top_peak, top_peak, bottom_peak, bottom_peak, top_peak];
    color=:white,
    linewidth=2.5,
)
scatter!(ax_peak, [cx_peak], [cy_peak]; color=:tomato, markersize=16)
Colorbar(fig_peak[1, 2], hm, label="mean density")
save(path_plot_peak, fig_peak)

println("name = ", name)
println("val = ", val)
println("variation = ", variation)
println("permuted /dens size = ", size(dens))
println("corner subtraction = ", (corner_height, corner_width))
println("dens_by_variation size = ", size(dens_by_variation))
println("Indexing order: dens_by_variation[repeat, t_hold, istp, height, width]")
println("val_number size = ", size(val_number))
println("err_number size = ", size(err_number))
println("plot path = ", path_plot)
println("mean density size = ", size(dens_mean))
println("peak center (x, y) = ", (cx_peak, cy_peak))
println("peak crop smwh = ", smwh_peak)
println("peak plot path = ", path_plot_peak)
