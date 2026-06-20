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

tag = "SIMU-NTRC"
log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"01.[2026.06.01]"
step_t = 0.1798
step_in_μm = 0.2613
freq_query = 1:1:140

# commit
title = "Anlz.08.Simu-01.[2025.06.02].[30-100ms]"
path_test = joinpath(path_root, dir_test)
path_this = @__FILE__
path_output = joinpath(path_root, title)
isdir(path_output) || mkpath(path_output)
cp(path_this, joinpath(path_output, basename(path_this)); force=true)

pattern_filename_data = Regex(raw"nxy_t=(?<idx_time>\d+).mat")
filename_xy = "XY.mat"

filenames_data = readdir(path_test) |> fs -> filter(f -> occursin(pattern_filename_data, f), fs)
ids_t = Vector{Int}(undef, length(filenames_data))
(w, h) = wh = (256, 256)
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
        x_mat = read(file_xy, "Y")
        y_mat = read(file_xy, "X")
        close(file_xy)
        (
            x_mat[1, :] |> diff |> unique |> x -> x[1],
            y_mat[2, :] |> diff |> unique |> y -> y[1]
        )
    end
px_in_um = (0.2344, 0.7812) .* step_in_μm
smwh_core = smwh_roi = (80, 125)
step_posi = px_in_um
step_modl = 1 ./ (2 .* smwh_roi .* px_in_um)
x_modl, y_modl = smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_modl
x_posi, y_posi = smwh_roi |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_posi


proc_sidepeak = true
proc_envelope = true
selector_moment = y -> (y .> 0.10) .& (y .< 0.60)
selector_sidepeak = y -> (y .> 0.1) .& (y .< 0.6)
selector_t_sidepeak = t -> 30 .< t .< 100
selector_t_envelope = t -> 30 .< t .< 100
selector_tail_stack = y -> y .> 0.02

perm_t = sortperm(ids_t)
ids_t = ids_t[perm_t]
dens_raw = dens_raw[perm_t, :, :, :]
istp = ["162", "164"]

sel_t = 1:2:201 # ids_t |> ids -> findall(i -> (i * step_t > 100) & (mod(i, 12) == 0), ids)
val_vars = (;
    t_hold=ids_t[sel_t] .* step_t,
    istp,
)
n_dim_vars = val_vars |> vs -> map(length, vs) |> Tuple
xy_peak_roi = (w, h) |> s -> map(s -> round((s + 1) / 2) |> Int, s)
dens_full_fmt = [copy(@view dens_raw[t, i, :, :]) |> transpose for t in axes(dens_raw, 1)[sel_t], i in axes(dens_raw, 2)] |>
                ds -> map(d -> crop_center(d, xy_peak_roi, smwh_roi), ds)
info_fmt = [
    Dict(
        "repeat" => 1,
        "istp" => istp[i],
        "t_hold" => ids_t[t] * step_t,
    )
    for t in 1:n_dim_vars[1], i in 1:n_dim_vars[2]
]

t_stage = log_step("calculating solo essentials for $(length(dens_full_fmt)) shots")
essn_2d_fmt = map(
    d -> calc_solo_essn_2d(d, smwh_roi .+ 1, smwh_roi, px_in_um, smwh_roi .+ 1, smwh_core),
    dens_full_fmt,
)
log_done("calculated solo essentials", t_stage)
essn_stacked_over_rep_t = [
    begin
        essns_t = [essn_2d_fmt[t, i] for t in axes(essn_2d_fmt, 1)] |> vec
        print("\r  [$tag] stacking over t istp_idx=$i n=$(length(essns_t))")
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
    e -> fit_prfl_modl_twinpeak_decay_1d(y_modl, e.prfl_modl.side.normed_px, selector_tail_stack(y_modl))
    for istp in axes(essn_stacked_over_rep_t, 1)
]
log_done("fit stacked modulation tails", t_stage)

t_stage = log_step("extracting per-shot sidepeak/envelope values")
extr_fmt = [
    begin
        print("\r  [$tag] extracting shots t_idx=$t istp_idx=$i")
        flush(stdout)
        essn_2d_fmt[t, i] |> e -> calc_solo_extr(
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

t_stage = log_step("analyzing per-shot trends")
trend_sidepeak_nvlp = [
    extr_fmt[:, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_sidepeak, selector_t_envelope)
    for r in axes(extr_fmt, 2), i in axes(extr_fmt, 4)
]
log_done("analyzed per-shot trends", t_stage)

t_stage = log_step("building trend figures for $tag")
fig_trend, axs_trend = begin
    fig = Figure()
    axs = Array{Dict}(undef, length(val_vars.istp))
    for i = 1:length(val_vars.istp)
        gl = GridLayout()
        fig[1, i] = gl
        axs[i] = set_panel_trend_sidepeak_nvlp!(gl, i)
    end
    (fig, axs)
end
log_done("built trend figures for $tag", t_stage)
t_plot_stage = log_step("plotting and saving trends for $tag ")
for i in 1:length(val_vars.istp)
    trend = trend_sidepeak_nvlp[i]
    val_istp = val_vars.istp[i]
    plot_trends_sidepeak!(axs_trend[i], trend, val_istp; to_clean=true, alpha=1.0, to_legend=true)
    plot_trends_nvlp!(axs_trend[i], trend, val_istp; to_clean=true, alpha=1.0, to_legend=true)
end
resize_to_layout!(fig_trend)
for format in ["pdf", "png"]
    fig_trend |> f -> save(joinpath(path_output, @sprintf("%s_trend.%s", tag, format)), f; backend=CairoMakie)
end
log_done("saved trends for $tag", t_plot_stage)

# fig_full, axs_solo = set_axis_full((1, n_dim_vars...), set_panel_solo_modl!; to_plot_stacked=false)
# println("Full axes ready: dimensions $(n_dim_vars)")
# for t in 1:n_dim_vars[1], i in 1:n_dim_vars[2]
#     info = info_fmt[t, i]
#     print("\r\033[2Kplotting for $(info["t_hold"]) ms, $(info["istp"])")
#     draw_solo_modl!(axs_solo[1, t, i], extr_fmt[t, i], info; dens_max=64.0, peak_height_max=3.0)
# end
# println("Full modulation table drawn.")
# resize_to_layout!(fig_full)
# fig_full |> f -> save(joinpath(path_output, @sprintf("%s_essn_table.pdf", tag)), f; backend=CairoMakie)
# println("Full modulation plot saved.")
