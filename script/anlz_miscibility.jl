using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
using Statistics
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
include(joinpath(@__DIR__, "..", "src", "visduet.jl"))

year_test = 2026
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\SingDrplMisc"
istp = ["162", "164"]
runinfos = [
    (date="0513", runids=71:75, tag_head="ImbaEvol", vars=(IB=5.378, rep=1:5, bias=0.1:0.05:0.6, t_hold=6:5:56, istp)),
    (date="0513", runids=76:80, tag_head="ImbaEvol", vars=(IB=5.376, rep=1:5, bias=0.1:0.05:0.6, t_hold=6:5:56, istp)),
    (date="0513", runids=81:82, tag_head="ImbaEvol", vars=(IB=[5.392, 5.386, 5.378, 5.372], rep=1:1, bias=0.1:0.1:0.6, t_hold=6:10:116, istp)),
]
title_anlz = "[05.21].05.StackedDuet"
path_output = joinpath(path_root, "AnlzRoutine", title_anlz);
if !isdir(path_output)
    mkpath(path_output)
end

wh_corner = (10, 10)
smwh_roi = (30, 30)
smwh_core = (20, 20)
smwh_strip = (2, 20)
wh_peak = smwh_roi .* 2 .+ 1
smw_peak, smh_peak = smwh_roi
smw_ft = 5
px_in_um = 6.5 / 22.06
step_posi = px_in_um
step_modl = 1 / (2 * smwh_roi[2] * px_in_um)
x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi
x_modl, y_modl = (x_vec, y_vec) .* step_modl

for (idx_runinfo, runinfo) in enumerate(runinfos)
    println("Processing set $idx_runinfo: $(gen_run_tag(runinfo))")
    global r = format_dens_runinfo(runinfo; path_root, year_test, wh_corner, smwh_roi)
    println("  val lengths ($(join(string.(propertynames(r.val)), ", "))): $(map(length, r.val))")
    println("  dens_full_fmt size: $(size(r.dens_full_fmt))")
    println("  image size: $(size(first(r.dens_full_fmt)))")
    println("  xy_peak_px: $(r.xy_peak_px), wh_dens: $(r.wh_dens)")
    global xy_peak_duet = r.dens_full_fmt |>
                          ds -> mapslices(
        imgs -> mean(imgs) |>
                d -> fit_dens2d_gaussian_round_disk(1:wh_peak[1], 1:wh_peak[2], d, :)["params"] |>
                     p -> (round(Int, p[2]), round(Int, p[3])),
        ds;
        dims=ndims(ds),
    ) |> p -> repeat(p, inner=ntuple(i -> i == 5 ? 2 : 1, length(r.n_dim_vars)))
    global essn_2d_fmt = map((d, xy) -> calc_solo_essn_2d(d, smwh_roi .+ 1, smwh_roi, smw_ft, px_in_um, xy, smwh_core; smwh_strip), r.dens_full_fmt, xy_peak_duet)
    global extr_2d_fmt = essn_2d_fmt |> es -> map(e -> calc_solo_extr(e, nothing; proc_sidepeak=false, proc_envelope=true), es)
    global info_fmt = [Dict("istp" => r.val.istp[i], "t_hold" => r.val.t_hold[t], "repeat" => r.val.rep[rep], "ib" => r.val.IB[c], "bias" => r.val.bias[b])
                       for c in 1:r.n_dim_vars[1], rep in 1:r.n_dim_vars[2], b in 1:r.n_dim_vars[3], t in 1:r.n_dim_vars[4], i in 1:r.n_dim_vars[5]]
    # Statistics on number sum
    global num_fmt = sum.(r.dens_full_fmt)
    global stat_n_fmt = num_fmt |> a -> mapslices(calc_mean_std, a; dims=(2))
    global essn_2d_stacked_over_rep = essn_2d_fmt |> es -> mapslices(calc_stacked_essn, es; dims=2)
    global info_stacked_over_rep = info_fmt[:, 1:1, :, :, :]
    fig_sizes, axs_sizes = set_axes_2axes!(r.runinfo.vars |> NamedTuple{(:IB, :bias)}, set_panel_single_axis, r.runinfo)
    for (c, ib) in enumerate(r.val.IB), (b, bias) in enumerate(r.val.bias)
        ax = axs_sizes[c, b]["ax"]
        [ax] |> clear_axes!
        for (i, istp) in enumerate(r.val.istp), rep in r.val.rep
            clr_theme = Oklch(0.52, 0.14, hue_theme_istp[istp])
            sizes = extr_2d_fmt[c, rep, b, :, i] |> es -> map(e -> e.envelope.params_round["size"], es)
            lines!(ax, r.val.t_hold, sizes; color=(clr_theme, 0.65))
        end
        ylims!(ax, 0, 6.0)
    end
    fig_sizes |> resize_to_layout!
    fig_sizes |> f -> save(joinpath(path_output, @sprintf("%s_sizes_t.png", gen_run_tag(r.runinfo))), f; backend=CairoMakie)
    fig_sizes, axs_sizes = set_axes_2axes!(r.runinfo.vars |> NamedTuple{(:IB, :t_hold)}, set_panel_single_axis, r.runinfo)
    for (c, ib) in enumerate(r.val.IB), (t, t_hold) in enumerate(r.val.t_hold)
        ax = axs_sizes[c, t]["ax"]
        [ax] |> clear_axes!
        for (i, istp) in enumerate(r.val.istp), rep in r.val.rep
            clr_theme = Oklch(0.52, 0.14, hue_theme_istp[istp])
            sizes = extr_2d_fmt[c, rep, :, t, i] |> es -> map(e -> e.envelope.params_round["size"], es)
            lines!(ax, r.val.bias, sizes; color=(clr_theme, 0.65))
        end
        ylims!(ax, 0, 6.0)
    end
    fig_sizes |> resize_to_layout!
    fig_sizes |> f -> save(joinpath(path_output, @sprintf("%s_sizes_bias.png", gen_run_tag(r.runinfo))), f; backend=CairoMakie)
    # for c in 1:r.n_dim_vars[1], b in 1:r.n_dim_vars[3]
    #     tag = @sprintf("Top View Number Stat [IB = %.3fA | bias = %.2f]", r.val.IB[c], r.val.bias[b])
    #     fig_num, axs_num = set_axis!(tag)
    #     [axs_num] |> clear_axes!
    #     for istp in 1:r.n_dim_vars[5]
    #         plot_num_stat_evo!(axs_num, r.val.t_hold, stat_n_fmt[c, 1, b, :, istp], r.val.istp[istp])
    #     end
    #     ylims!(axs_num, 0, 8000.0)
    #     fig_num |> f -> save(joinpath(path_output, @sprintf("%s_num_stat_[IB=%.3fA'bias=%.2f].png", gen_run_tag(r.runinfo), r.val.IB[c], r.val.bias[b])), f; backend=CairoMakie)
    # end
    # println("\r\033[2K\rNow drawing table for stacked over rep.")
    # for c in 1:r.n_dim_vars[1]
    #     global fig_stacked_duets, axs_stacked_duets = set_axes_v_t_rep!(Base.setindex(r.n_dim_vars, 1, 2)[2:end], set_panel_misc_duet_2d!, r.runinfo, info_stacked_over_rep[c, :, :, :, :]; partidx=c)
    #     for b in 1:r.n_dim_vars[3], t in 1:r.n_dim_vars[4]
    #         draw_misc_duet_core_2d!(axs_stacked_duets[1, b, t], essn_2d_stacked_over_rep[c, 1, b, t, :])
    #         print("\r\033[2K\rdrawing duet at, $b, $t.")
    #     end
    #     println("\r\033[2K\rdrawing complete for $c.")
    #     fig_stacked_duets |> resize_to_layout!
    #     fig_stacked_duets |> f -> save(joinpath(path_output, @sprintf("%s_[IB=%.3fA]_essn_table_stacked.pdf", gen_run_tag(r.runinfo), r.val.IB[c])), f; backend=CairoMakie)
    # end
    # println("\r\033[2K\rNow drawing table for full run.")
    # for c in 1:r.n_dim_vars[1]
    #     global fig_full_duets, axs_full_duets = set_axes_v_t_rep!(Tuple(r.n_dim_vars)[2:end], set_panel_misc_duet_2d!, r.runinfo, info_fmt[c, :, :, :, :]; partidx=c)
    #     for rep in 1:r.n_dim_vars[2], b in 1:r.n_dim_vars[3], t in 1:r.n_dim_vars[4]
    #         draw_misc_duet_2d!(axs_full_duets[rep, b, t], essn_2d_fmt[c, rep, b, t, :])
    #         print("\r\033[2K\rdrawing duet at $rep, $b, $t.")
    #     end
    #     println("\r\033[2K\rdrawing complete for $c.")
    #     fig_full_duets |> resize_to_layout!
    #     fig_full_duets |> f -> save(joinpath(path_output, @sprintf("%s_[IB=%.3fA]_essn_table.pdf", gen_run_tag(r.runinfo), r.val.IB[c])), f; backend=CairoMakie)
    # end
end
