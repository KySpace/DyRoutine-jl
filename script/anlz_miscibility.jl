using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
using Statistics
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))
include(joinpath(@__DIR__, "..", "src", "visduet.jl"))

year_test = 2026
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\SingDrplMisc"
runinfos = [
    (date="0513", runids=71:75, IB=5.378, tag_head="ImbaEvol", rep_each=1, bias=0.1:0.05:0.6, t_hold=6:5:56),
    (date="0513", runids=76:80, IB=5.376, tag_head="ImbaEvol", rep_each=1, bias=0.1:0.05:0.6, t_hold=6:5:56),
    (date="0513", runids=81, IB=[5.392, 5.386], tag_head="ImbaEvol", rep_each=1, bias=0.1:0.1:0.6, t_hold=6:10:116),
    (date="0513", runids=82, IB=[5.378, 5.372], tag_head="ImbaEvol", rep_each=1, bias=0.1:0.1:0.6, t_hold=6:10:116),
][1:2]
n_istp = 2
title_anlz = "[05.15].01.DevTest"
path_output = joinpath(path_root, "AnlzRoutine", title_anlz);
name_dims = ["IB" "repeat" "bias" "t_hold" "istp"]
val_istp = ["162", "164"]
if !isdir(path_output)
    mkpath(path_output)
end

wh_corner = (10, 10)
smwh_peak = (30, 30)
wh_peak = smwh_peak .* 2 .+ 1
smw_peak, smh_peak = smwh_peak
smw_ft = 5
px_in_um = 6.5 / 22.06
step_posi = px_in_um
step_modl = 1 / (2 * smwh_peak[2] * px_in_um)
x_vec, y_vec = smwh_peak |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi
x_modl, y_modl = (x_vec, y_vec) .* step_modl

as_vector(x) = x isa AbstractArray ? collect(x) : [x]

function gen_run_tag(runinfo)
    str_runids = runinfo.runids |> a -> "$(a)" |> s -> replace(s, ":" => "-")
    ib_tag = as_vector(runinfo.IB) |> ib -> join(map(v -> @sprintf("%.3f", v), ib), "_")
    return @sprintf("%s_IB%s_r%s", runinfo.tag_head, ib_tag, str_runids)
end

function load_dens_run(date::AbstractString, runid::Integer)
    dir_data = gen_date_path(date, year_test)
    file_data = gen_h5name(date, runid)
    path_input_routine = joinpath(path_root, dir_data, @sprintf("run%02d", runid), file_data)
    path_input_misc = joinpath(path_root, @sprintf("run%02d", runid), file_data)
    path_input = isfile(path_input_routine) ? path_input_routine : path_input_misc

    h5open(path_input, "r") do f
        dens_run = f["/od"] |>
                   read |>
                   x_vec -> permutedims(x_vec, (3, 2, 1)) |>
                            x_vec -> stack(
                       map(d -> subtract_corner_mean(d, wh_corner), eachslice(x_vec; dims=1));
                       dims=1,
                   )
        ndims(dens_run) == 3 || error("Expected /od in $path_input to have 3 dimensions after formatting, got $(ndims(dens_run)).")
        return dens_run
    end
end

function format_dens_runinfo(runinfo)
    runids = as_vector(runinfo.runids)
    val_ib = as_vector(runinfo.IB)
    val_bias = as_vector(runinfo.bias)
    val_thold = as_vector(runinfo.t_hold)

    n_ib = length(val_ib)
    n_rep = runinfo.rep_each * length(runids)
    n_bias = length(val_bias)
    n_thold = length(val_thold)
    n_istp = length(val_istp)

    dens = map(runid -> load_dens_run(runinfo.date, runid), runids) |>
           ds -> cat(ds...; dims=1)
    n_shot, h_dens, w_dens = size(dens)
    val = (
        ib=val_ib,
        rep=collect(1:n_rep),
        bias=val_bias,
        t_hold=val_thold,
        istp=val_istp,
    )
    n_dim_vars = map(length, val)
    n_variation = prod(n_dim_vars)
    n_shot == n_variation || throw(DimensionMismatch("Loaded $n_shot shots for $(gen_run_tag(runinfo)), but expected $n_variation from dimensions $n_dim_vars."))

    dens_mean = dropdims(mean(dens; dims=1); dims=1)
    xy_peak_px = find_positive_cluster_center(dens_mean; smwh=smwh_peak) |> cent -> round.(Int, cent)
    dens_full_fmt = dens |>
                    ds -> mapslices(d -> crop_center(d, xy_peak_px, smwh_peak), ds; dims=(2, 3)) |>
                          ds -> reshape(ds, (reverse(n_dim_vars)..., reverse(wh_peak)...)) |>
                                ds -> permutedims(ds, (5, 4, 3, 2, 1, 6, 7))

    # A lite version for tests
    # runinfo_lite = (date=runinfo.date, runids=runinfo.runids[1:3], IB=runinfo.IB, tag_head="ImbaEvol", rep_each=1, bias=runinfo.bias[5:7], t_hold=runinfo.t_hold[1:4])
    # val_lite = (
    #     ib=val_ib,
    #     rep=collect(1:3),
    #     bias=val_bias[5:7],
    #     t_hold=val_thold[1:4],
    #     istp=val_istp,
    # )
    # n_dim_vars = map(length, val_lite)
    # n_variation = prod(n_dim_vars)
    # dens_full_fmt = dens_full_fmt[:, 1:3, 5:7, 1:4, :, :, :]

    # size(dens_full_fmt)[1:5] == (n_dim_vars |> Tuple) || throw(DimensionMismatch("Formatted dimensions $(size(dens_full_fmt)[1:5]) do not match expected $n_dim_vars for $(gen_run_tag(runinfo))."))
    # return (; runinfo=runinfo_lite, val=val_lite, dens_full_fmt, wh_dens=(w_dens, h_dens), xy_peak_px, n_dim_vars)
    return (; runinfo, val, dens_full_fmt, wh_dens=(w_dens, h_dens), xy_peak_px, n_dim_vars)
end

for (idx_runinfo, runinfo) in enumerate(runinfos)
    println("Processing set $idx_runinfo: $(gen_run_tag(runinfo))")
    global r = format_dens_runinfo(runinfo)
    println("  val lengths ($(join(vec(name_dims), ", "))): $(map(length, r.val))")
    println("  dens_full_fmt size: $(size(r.dens_full_fmt))")
    println("  xy_peak_px: $(r.xy_peak_px), wh_dens: $(r.wh_dens)")
    global essn_2d_fmt = r.dens_full_fmt |> ds -> mapslices(d -> calc_solo_essn_2d(d, smwh_peak .+ 1, smwh_peak, smw_ft, px_in_um), ds; dims=(6, 7)) |> e -> dropdims(e; dims=(6, 7))
    global info_fmt = [Dict("istp" => r.val.istp[i], "t_hold" => r.val.t_hold[t], "repeat" => r.val.rep[rep], "ib" => r.val.ib[c], "bias" => r.val.bias[b])
                for c in 1:r.n_dim_vars[1], rep in 1:r.n_dim_vars[2], b in 1:r.n_dim_vars[3], t in 1:r.n_dim_vars[4], i in 1:r.n_dim_vars[5]]
    for c in 1:r.n_dim_vars[1]
        global fig_full_duets, axs_full_duets = set_axes_v_t_rep!(Tuple(r.n_dim_vars)[2:end], set_panel_misc_duet_2d!, r.runinfo, info_fmt[c, :, :, :, :])
        for rep in 1:r.n_dim_vars[2], b in 1:r.n_dim_vars[3], t in 1:r.n_dim_vars[4]
            draw_misc_duet_2d!(axs_full_duets[rep, b, t], essn_2d_fmt[c, rep, b, t, :])
            print("\r\033[2K\rdrawing duet at $rep, $b, $t.")
        end
        println("\r\033[2K\rdrawing complete for $c.")
        fig_full_duets |> resize_to_layout!
        fig_full_duets |> f -> save(joinpath(path_output, @sprintf("%s_essn_table.pdf", gen_run_tag(runinfo))), f; backend=CairoMakie)
    end
end
