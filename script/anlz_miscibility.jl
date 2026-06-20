println("Processing set $idx_runinfo: $(gen_run_tag(runinfo))")
fmt_dens = format_dens_runinfo(runinfo; path_root, year_test, wh_corner, smwh_roi, len_avg_peak, sel_vars)
runinfo = fmt_dens.runinfo
val_vars = fmt_dens.val_vars
(; dens_full_fmt, wh_dens, xy_peak_px, n_dim_vars) = fmt_dens
println("  val_vars lengths ($(join(string.(propertynames(val_vars)), ", "))): $(map(length, val_vars))")
println("  dens_full_fmt size: $(size(dens_full_fmt))")
println("  image size: $(size(first(dens_full_fmt)))")
println("  xy_peak_px: $(xy_peak_px), wh_dens: $(wh_dens)")

xy_peak_duet = dens_full_fmt |>
               ds -> mapslices(
                    imgs -> mean(imgs) |>
                            d -> fit_dens2d_gaussian_round_disk(1:wh_peak[1], 1:wh_peak[2], d, :; fit_round_kwargs...).params |>
                                p -> (round(Int, p[2]), round(Int, p[3])),
                    ds;
                    dims=ndims(ds),
                ) |> p -> repeat(p, inner=ntuple(i -> i == idx_istp_axis ? n_istp_per_condition : 1, length(n_dim_vars)))

essn_2d_fmt = map(
    (d, xy) -> calc_solo_essn_2d(d, smwh_roi .+ 1, smwh_roi, px_in_um, xy, smwh_core; smwh_strip),
    dens_full_fmt,
    xy_peak_duet,
)
extr_2d_fmt = essn_2d_fmt |>
              es -> map(
    e -> calc_solo_extr(
        e,
        nothing;
        proc_sidepeak,
        proc_envelope,
        selector_moment,
        selector_sidepeak,
        fit_asymm_kwargs,
        fit_round_kwargs,
    ),
    es,
)
info_fmt = [
    Dict("istp" => val_vars.istp[i], "t_hold" => val_vars.t_hold[t], "repeat" => val_vars.rep[rep], "ib" => val_vars.IB[c], "bias" => val_vars.bias[b])
    for c in 1:n_dim_vars[1], rep in 1:n_dim_vars[2], b in 1:n_dim_vars[3], t in 1:n_dim_vars[4], i in 1:n_dim_vars[5]
]

# Statistics on number sum
num_fmt = sum.(dens_full_fmt)
stat_n_fmt = num_fmt |> a -> mapslices(calc_mean_std, a; dims=(idx_rep_axis))
essn_2d_stacked_over_rep = essn_2d_fmt |> es -> mapslices(calc_stacked_essn, es; dims=idx_rep_axis)
info_stacked_over_rep = info_fmt[:, 1:1, :, :, :]

fig_sizes, axs_sizes = set_axes_2axes!(runinfo.vars |> NamedTuple{(:IB, :bias)}, set_panel_single_axis, runinfo)
for (c, ib) in enumerate(val_vars.IB), (b, bias) in enumerate(val_vars.bias)
    ax = axs_sizes[c, b]["ax"]
    [ax] |> clear_axes!
    for (i, istp_iter) in enumerate(val_vars.istp), rep in val_vars.rep
        clr_theme = Oklch(0.52, 0.14, hue_theme_istp[istp_iter])
        sizes = extr_2d_fmt[c, rep, b, :, i] |> es -> map(e -> e.envelope.params_round.size, es)
        lines!(ax, val_vars.t_hold, sizes; color=(clr_theme, 0.65))
    end
    ylims!(ax, 0, 6.0)
end
fig_sizes |> resize_to_layout!
fig_sizes |> f -> save(joinpath(path_output, @sprintf("%s_sizes_t.png", gen_run_tag(runinfo))), f; backend=CairoMakie)

fig_sizes, axs_sizes = set_axes_2axes!(runinfo.vars |> NamedTuple{(:IB, :t_hold)}, set_panel_single_axis, runinfo)
for (c, ib) in enumerate(val_vars.IB), (t, t_hold_iter) in enumerate(val_vars.t_hold)
    ax = axs_sizes[c, t]["ax"]
    [ax] |> clear_axes!
    for (i, istp_iter) in enumerate(val_vars.istp), rep in val_vars.rep
        clr_theme = Oklch(0.52, 0.14, hue_theme_istp[istp_iter])
        sizes = extr_2d_fmt[c, rep, :, t, i] |> es -> map(e -> e.envelope.params_round.size, es)
        lines!(ax, val_vars.bias, sizes; color=(clr_theme, 0.65))
    end
    ylims!(ax, 0, 6.0)
end
fig_sizes |> resize_to_layout!
fig_sizes |> f -> save(joinpath(path_output, @sprintf("%s_sizes_bias.png", gen_run_tag(runinfo))), f; backend=CairoMakie)

# for c in 1:n_dim_vars[1], b in 1:n_dim_vars[3]
#     tag = @sprintf("Top View Number Stat [IB = %.3fA | bias = %.2f]", val_vars.IB[c], val_vars.bias[b])
#     fig_num, axs_num = set_axis!(tag)
#     [axs_num] |> clear_axes!
#     for istp in 1:n_dim_vars[5]
#         plot_num_stat_evo!(axs_num, val_vars.t_hold, stat_n_fmt[c, 1, b, :, istp], val_vars.istp[istp])
#     end
#     ylims!(axs_num, 0, 8000.0)
#     fig_num |> f -> save(joinpath(path_output, @sprintf("%s_num_stat_[IB=%.3fA'bias=%.2f].png", gen_run_tag(runinfo), val_vars.IB[c], val_vars.bias[b])), f; backend=CairoMakie)
# end
# println("\r\033[2K\rNow drawing table for stacked over rep.")
# for c in 1:n_dim_vars[1]
#     fig_stacked_duets, axs_stacked_duets = set_axes_v_t_rep!(Base.setindex(n_dim_vars, 1, 2)[2:end], set_panel_misc_duet_2d!, runinfo, info_stacked_over_rep[c, :, :, :, :]; partidx=c)
#     for b in 1:n_dim_vars[3], t in 1:n_dim_vars[4]
#         draw_misc_duet_core_2d!(axs_stacked_duets[1, b, t], essn_2d_stacked_over_rep[c, 1, b, t, :])
#         print("\r\033[2K\rdrawing duet at, $b, $t.")
#     end
#     println("\r\033[2K\rdrawing complete for $c.")
#     fig_stacked_duets |> resize_to_layout!
#     fig_stacked_duets |> f -> save(joinpath(path_output, @sprintf("%s_[IB=%.3fA]_essn_table_stacked.pdf", gen_run_tag(runinfo), val_vars.IB[c])), f; backend=CairoMakie)
# end
# println("\r\033[2K\rNow drawing table for full run.")
# for c in 1:n_dim_vars[1]
#     fig_full_duets, axs_full_duets = set_axes_v_t_rep!(Tuple(n_dim_vars)[2:end], set_panel_misc_duet_2d!, runinfo, info_fmt[c, :, :, :, :]; partidx=c)
#     for rep in 1:n_dim_vars[2], b in 1:n_dim_vars[3], t in 1:n_dim_vars[4]
#         draw_misc_duet_2d!(axs_full_duets[rep, b, t], essn_2d_fmt[c, rep, b, t, :])
#         print("\r\033[2K\rdrawing duet at $rep, $b, $t.")
#     end
#     println("\r\033[2K\rdrawing complete for $c.")
#     fig_full_duets |> resize_to_layout!
#     fig_full_duets |> f -> save(joinpath(path_output, @sprintf("%s_[IB=%.3fA]_essn_table.pdf", gen_run_tag(runinfo), val_vars.IB[c])), f; backend=CairoMakie)
# end
