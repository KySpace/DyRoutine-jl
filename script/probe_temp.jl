using HDF5
using Statistics: mean
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
include(joinpath(@__DIR__, "..", "src", "pershot.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\2026-03\0324\run64\d0324r64.h5"
path_plot = joinpath(@__DIR__, "probe_temp_number_vs_t_hold.svg")
path_plot_peak = joinpath(@__DIR__, "probe_temp_avg_density_peak.svg")
path_plot_duet = joinpath(@__DIR__, "probe_temp_duet.svg")
path_plot_sheet = joinpath(@__DIR__, "probe_temp_contact_sheet.pdf")
corner_height = 10
corner_width = 10
smwh_peak = (30, 60)
duet_color_max = 40.0
sheet_color_max = 60.0
sheet_gap = 6

name = ["repeat", "t_hold", "istp"]
val = (
    collect(1:3),
    collect(6:2:200),
    ["162", "164"],
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
xy_peak_px = round.(Int, (cx_peak, cy_peak))
cropy_mean = crop_center(dens_mean, xy_peak_px, smwh_peak)

smw_peak, smh_peak = smwh_peak
left_peak = xy_peak_px[1] - smw_peak
right_peak = xy_peak_px[1] + smw_peak
top_peak = xy_peak_px[2] - smh_peak
bottom_peak = xy_peak_px[2] + smh_peak

fig_peak = Figure(size=(1320, 760))
ax_peak_full = Axis(
    fig_peak[1, 1];
    title="Average Density with Peak-Finding Window",
    xlabel="x",
    ylabel="y",
    yreversed=true,
    aspect=DataAspect(),
)
hm_full = heatmap!(ax_peak_full, 1:width, 1:height, dens_mean'; colormap=:viridis)
lines!(
    ax_peak_full,
    [left_peak, right_peak, right_peak, left_peak, left_peak],
    [top_peak, top_peak, bottom_peak, bottom_peak, top_peak];
    color=:white,
    linewidth=2.5,
)
scatter!(ax_peak_full, [cx_peak], [cy_peak]; color=:tomato, markersize=16)

ax_peak_crop = Axis(
    fig_peak[1, 2];
    title="Cropped Mean Density",
    xlabel="x",
    ylabel="y",
    yreversed=true,
    aspect=DataAspect(),
)
hm_crop = heatmap!(
    ax_peak_crop,
    left_peak:right_peak,
    top_peak:bottom_peak,
    Matrix(cropy_mean)';
    colormap=:viridis,
)
scatter!(ax_peak_crop, [cx_peak], [cy_peak]; color=:tomato, markersize=16)

Colorbar(fig_peak[1, 3], hm_full, label="mean density")
save(path_plot_peak, fig_peak)

duet_a = crop_center(@view(dens_by_variation[1, 1, 1, :, :]), xy_peak_px, smwh_peak)
duet_b = crop_center(@view(dens_by_variation[1, 1, 2, :, :]), xy_peak_px, smwh_peak)

fig_duet = Figure(size=(1180, 620))
ax_duet_a = Axis(
    fig_duet[1, 1];
    title="Duet Shot 1",
    xlabel="x",
    ylabel="y",
    yreversed=true,
    aspect=DataAspect(),
)
hm_duet_a = heatmap!(
    ax_duet_a,
    left_peak:right_peak,
    top_peak:bottom_peak,
    Matrix(duet_a)';
    colormap=:RdPu_9,
    colorrange=(0, duet_color_max),
)

ax_duet_b = Axis(
    fig_duet[1, 2];
    title="Duet Shot 2",
    xlabel="x",
    ylabel="y",
    yreversed=true,
    aspect=DataAspect(),
)
hm_duet_b = heatmap!(
    ax_duet_b,
    left_peak:right_peak,
    top_peak:bottom_peak,
    Matrix(duet_b)';
    colormap=:PuBu_9,
    colorrange=(0, duet_color_max),
)

Colorbar(fig_duet[1, 3], hm_duet_a, label="shot 1 density")
Colorbar(fig_duet[1, 4], hm_duet_b, label="shot 2 density")
save(path_plot_duet, fig_duet)

sheet_colormaps = Dict(1 => :RdPu_9, 2 => :PuBu_9)
sheet_labels = Dict(1 => string(val[3][1]), 2 => string(val[3][2]))
sheet_crop_size = 2 .* collect(reverse(smwh_peak)) .+ 1
sheet_ncols = length(val[1]) * length(val[3])
fig_sheet = Figure(
    size=(
        round(Int, sheet_ncols * sheet_crop_size[2] * 1.35 + (sheet_ncols - 1) * sheet_gap + 220),
        round(Int, length(val[2]) * sheet_crop_size[1] * 1.25 + (length(val[2]) - 1) * sheet_gap + 220),
    ),
)

for (row_idx, t_hold) in enumerate(val[2])
    for (col_idx, repeat_val) in enumerate(val[1])
        for istp_idx in eachindex(val[3])
            sheet_col_idx = (col_idx - 1) * length(val[3]) + istp_idx
            cropped_panel = crop_center(@view(dens_by_variation[col_idx, row_idx, istp_idx, :, :]), xy_peak_px, smwh_peak)
            ax_panel = Axis(
                fig_sheet[row_idx, sheet_col_idx];
                xticksvisible=false,
                yticksvisible=false,
                xticklabelsvisible=false,
                yticklabelsvisible=false,
                xgridvisible=false,
                ygridvisible=false,
                topspinevisible=false,
                rightspinevisible=false,
                bottomspinevisible=false,
                leftspinevisible=false,
                yreversed=true,
                aspect=DataAspect(),
            )
            heatmap!(
                ax_panel,
                left_peak:right_peak,
                top_peak:bottom_peak,
                Matrix(cropped_panel)';
                colormap=sheet_colormaps[istp_idx],
                colorrange=(0, sheet_color_max),
            )
            text!(
                ax_panel,
                (left_peak + right_peak) / 2,
                top_peak + 4;
                text="$(repeat_val) rep | $(t_hold) ms",
                align=(:center, :top),
                color=:black,
                fontsize=10,
            )
        end
    end
end

rowgap!(fig_sheet.layout, sheet_gap)
colgap!(fig_sheet.layout, sheet_gap)
save(path_plot_sheet, fig_sheet)

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
println("peak center pixel (x, y) = ", xy_peak_px)
println("peak crop smwh = ", smwh_peak)
println("cropped mean density size = ", size(cropy_mean))
println("peak plot path = ", path_plot_peak)
println("duet crop sizes = ", (size(duet_a), size(duet_b)))
println("duet colorrange = ", (0, duet_color_max))
println("duet plot path = ", path_plot_duet)
println("sheet colorrange = ", (0, sheet_color_max))
println("sheet gap = ", sheet_gap)
println("sheet path = ", path_plot_sheet)
