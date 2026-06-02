using HDF5: DatasetOrAttribute
using MAT
using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))


path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"01.[2025.06.01]"
step_t = 0.1798
step_in_μm = 0.2613

path_test = joinpath(path_root, dir_test)

pattern_filename_data = Regex(raw"nxy_t=(?<idx_time>\d+).mat")
filename_xy = "XY.mat"

filenames_data = readdir(path_test) |> fs -> filter(f -> occursin(pattern_filename_data, f), fs)
ids_t = Vector{Int}(undef, length(filenames_data))
dens_raw = Array{Float64}(undef, (length(filenames_data), 2, w, h))

for (i, fn) in enumerate(filenames_data)
    file = matopen(joinpath(path_test, fn))
    ids_t[i] = parse(Int, match(pattern_filename_data, fn)["idx_time"])
    dens_raw[i, 1, :, :] = read(file, "n1")
    dens_raw[i, 2, :, :] = read(file, "n2")
    close(file)
end
(step_x, step_y) =
    begin
        file_xy = matopen(joinpath(path_test, filename_xy))
        x_mat = read(file_xy, "X")
        y_mat = read(file_xy, "Y")
        close(file_xy)
        (
            x_mat[1, :] |> diff |> unique |> x -> x[1],
            y_mat[2, :] |> diff |> unique |> y -> y[1]
        )
    end
(w, h) = wh = (256, 256)
px_in_um = (0.7812, 0.2344) .* step_in_μm
smwh_core = smwh_roi = (50, 100)
step_posi = px_in_um
step_modl = 1 ./ (2 .* smwh_roi .* px_in_um)
smw_ft = 5
x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi
x_modl, y_modl = (x_vec, y_vec) .* step_modl


proc_sidepeak = true
proc_envelope = true
selector_moment = y -> (y .> 0.10) .& (y .< 0.50)
selector_sidepeak = y -> (y .> 0.1) .& (y .< 0.5)
selector_t_sidepeak = t -> 0 .< t .< 80
selector_t_envelope = t -> 0 .< t .< 80
selector_tail_stack = y -> y .> 0.02

perm_t = sortperm(ids_t)
ids_t = ids_t[perm_t]
dens_raw = dens_raw[perm_t, :, :, :]
istp = ["162", "164"]

var_vals = (;
    t_hold=ids_t .* step_t,
    istp,
)
xy_peak_roi = (w, h) |> s -> map(s -> round((s + 1) / 2), s)
dens_full_fmt = [copy(@view dens_raw[t, i, :, :]) for t in axes(dens_raw, 1), i in axes(dens_raw, 2)] |>
                ds -> map(d -> crop_center(xy_peak_roi, smwh_roi), ds)
info_fmt = [
    Dict(
        "repeat" => 1,
        "istp" => istp[i],
        "t_hold" => ids_t[t] * step_t,
    )
    for t in 1:n_dim_vars[3], i in 1:n_dim_vars[4]
]

t_stage = log_step("calculating solo essentials for $(length(dens_full_fmt)) shots")
essn_2d_fmt = map(
    d -> calc_solo_essn_2d(d, smwh_roi .+ 1, smwh_roi, smw_ft, px_in_um, smwh_roi .+ 1, smwh_core),
    dens_full_fmt,
)
log_done("calculated solo essentials", t_stage)
essn_stacked_over_rep_t = [
    begin
        essns_t = [essn_2d_fmt[t, i] for t in axes(essn_2d_fmt, 1)] |> vec
        print("\r  [$tag] stacking over t istp_idx=$i n=$(length(essns_rt))")
        flush(stdout)
        calc_stacked_essn(essns_t)
    end
    for i in axes(essn_2d_fmt, 2)
]
println()
log_done("stacked essentials", t_stage)
t_stage = log_step("fitting stacked modulation tails")
fit_prfl_modl_over_rep_t_1d = [
    essn_stacked_over_rep_t[istp] |>
    e -> fit_prfl_modl_twinpeak_decay_1d(y_modl, e.prfl_modl_norm_px, selector_tail_stack(y_modl))
    for istp in axes(essn_stacked_over_rep_t, 2)
]
log_done("fit stacked modulation tails", t_stage)

t_stage = log_step("extracting per-shot sidepeak/envelope values")
extr_fmt = [
    begin
        print("\r  [$tag] extracting shots IB_idx=$c rep-$r t_idx=$t istp_idx=$i")
        flush(stdout)
        essn_2d_fmt[c, r, t, i] |> e -> calc_solo_extr(
            e,
            fit_prfl_modl_over_rep_t_1d[i];
            proc_sidepeak,
            proc_envelope,
            selector_moment,
            selector_sidepeak,
        )
    end
    for i in axes(essn_2d_fmt, 2), t in axes(essn_2d_fmt, 1)
] |> e -> permutedims(e, reverse(1:ndims(e)))
println()
log_done("extracted per-shot sidepeak/envelope values", t_stage)

fig_live = Figure()
axs_162 = Axis(fig_live[1, 1], title="162")
axs_164 = Axis(fig_live[1, 2], title="164")
clrmap = [gen_clrmap_solo(hue_theme) for hue_theme in [hue_theme_istp["162"], hue_theme_istp["164"]]]
for (t, timestamp) in enumerate(ids_t), i in axes(dens_raw, 2)
    ax = i == 1 ? axs_162 : axs_164
    [ax] |> clear_axes!
    heatmap!(ax, imgs[t, i]; colorrange=(0, 36), colormap=clrmap[i])
    ax.aspect = DataAspect()
    fig_live |> resize_to_layout!
    if i == 2
        fig_live |> display
        sleep(0.2)
    end
end
