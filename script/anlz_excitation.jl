log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))

get_bind_date(runinfo, idx_bind) = hasproperty(runinfo, :date_runid) ? first(runinfo.date_runid[idx_bind]) : runinfo.date
get_bind_runid(runinfo, idx_bind) = if hasproperty(runinfo, :date_runid)
    last(runinfo.date_runid[idx_bind])
elseif hasproperty(runinfo, :runids)
    as_vector(runinfo.runids)[idx_bind]
else
    as_vector(runinfo.runid)[idx_bind]
end
get_bind_runinfo(runinfo, val_vars, idx_bind) = merge(
    runinfo,
    (;
        date=get_bind_date(runinfo, idx_bind),
        runid=get_bind_runid(runinfo, idx_bind),
        IB=val_vars.IB[idx_bind],
    ),
)

name = propertynames(runinfo.vars)
t_stage = log_step("formatting density data")
(; val_vars, dens_full_fmt, wh_dens, xy_peak_px, n_dim_vars, name_dims) = format_dens_runinfo(runinfo; path_root, year_test, wh_corner, smwh_roi, len_avg_peak)
n_variation = prod(n_dim_vars)
n_dim_vars_per_IB = Tuple(n_dim_vars[2:end])
n_IB, n_rep, n_main, n_istp = n_dim_vars
log_done("formatted density data: axes $(name_dims) dims $(n_dim_vars), per-IB dims $(n_dim_vars_per_IB), image $(size(first(dens_full_fmt)))", t_stage)

# Statistics on number sum
# num_fmt = sum.(dens_full_fmt)
# stat_n_fmt = num_fmt |> a -> mapslices(calc_mean_std, a; dims=(idx_rep_axis))

# fig_num, ax_num = set_axis!("number vs t hold")
# for (i, istp_iter) in enumerate(val_vars.istp)
#     plot_num_stat_evo!(ax_num, val_vars.t_hold, stat_n_fmt[1, :, i], istp_iter)
# end
# display(fig_num)

xy_peak_core = smwh_roi .+ 1
t_stage = log_step("calculating solo essentials for $(length(dens_full_fmt)) shots")
essn_2d_fmt = map(
    d -> calc_solo_essn_2d(d, smwh_roi .+ 1, smwh_roi, smw_ft, px_in_um, xy_peak_core, smwh_core),
    dens_full_fmt,
)
log_done("calculated solo essentials", t_stage)
info_fmt = [
    Dict(
        "istp" => val_vars.istp[i],
        "t_hold" => val_vars.t_hold[t],
        "repeat" => val_vars.rep[r],
        "runid" => get_bind_runid(runinfo, c),
        "IB" => val_vars.IB[c],
    )
    for c in 1:n_dim_vars[1], r in 1:n_dim_vars[2], t in 1:n_dim_vars[3], i in 1:n_dim_vars[4]
]

t_stage = log_step("stacking essentials over repeats and full time traces")
essn_stacked_over_rep = [
    begin
        print("\r  [$tag] stacking over rep IB_idx=$c t_idx=$t istp_idx=$i")
        flush(stdout)
        calc_stacked_essn(@view essn_2d_fmt[c, :, t, i])
    end
    for c in axes(essn_2d_fmt, 1), t in axes(essn_2d_fmt, 3), i in axes(essn_2d_fmt, 4)
]
println()
essn_stacked_over_rep_t = [
    begin
        essns_rt = [essn_2d_fmt[c, r, t, i] for r in axes(essn_2d_fmt, 2), t in axes(essn_2d_fmt, 3)] |> vec
        print("\r  [$tag] stacking over rep+t IB_idx=$c istp_idx=$i n=$(length(essns_rt))")
        flush(stdout)
        calc_stacked_essn(essns_rt)
    end
    for c in axes(essn_2d_fmt, 1), i in axes(essn_2d_fmt, 4)
]
println()
log_done("stacked essentials", t_stage)

t_stage = log_step("fitting stacked modulation tails")
fit_prfl_modl_over_rep_t_1d = [
    essn_stacked_over_rep_t[c, istp] |>
    e -> fit_prfl_modl_twinpeak_decay_1d(y_modl, e.prfl_modl_norm_px, selector_tail_stack(y_modl); fit_stack_kwargs...)
    for c in axes(essn_stacked_over_rep_t, 1), istp in axes(essn_stacked_over_rep_t, 2)
]
log_done("fit stacked modulation tails", t_stage)

t_stage = log_step("extracting per-shot sidepeak/envelope values")
extr_fmt = [
    begin
        # if r == first(axes(essn_2d_fmt, 2)) && (t == first(axes(essn_2d_fmt, 3)) || t % 25 == 0 || t == last(axes(essn_2d_fmt, 3)))
            print("\r  [$tag] extracting shots IB_idx=$c t_idx=$t istp_idx=$i")
            flush(stdout)
        # end
        essn_2d_fmt[c, r, t, i] |> e -> calc_solo_extr(
            e,
            fit_prfl_modl_over_rep_t_1d[c, i];
            proc_sidepeak,
            proc_envelope,
            selector_moment,
            selector_sidepeak,
            fit_tailess_kwargs,
            fit_asymm_kwargs,
            fit_round_kwargs,
        )
    end
    for i in axes(essn_2d_fmt, 4), t in axes(essn_2d_fmt, 3), r in axes(essn_2d_fmt, 2), c in axes(essn_2d_fmt, 1)
] |> e -> permutedims(e, reverse(1:ndims(e)))
println()
log_done("extracted per-shot sidepeak/envelope values", t_stage)

t_stage = log_step("extracting stacked-over-repeat values")
extr_stacked_over_rep = [
    begin
        if t == first(axes(essn_stacked_over_rep, 2)) || t % 25 == 0 || t == last(axes(essn_stacked_over_rep, 2))
            print("\r  [$tag] extracting stacked IB_idx=$c t_idx=$t istp_idx=$i")
            flush(stdout)
        end
        essn_stacked_over_rep[c, t, i] |> e -> calc_solo_extr(
            e,
            fit_prfl_modl_over_rep_t_1d[c, i];
            proc_sidepeak,
            proc_envelope,
            selector_moment,
            selector_sidepeak,
            fit_tailess_kwargs,
            fit_asymm_kwargs,
            fit_round_kwargs,
        )
    end
    for c in axes(essn_stacked_over_rep, 1), t in axes(essn_stacked_over_rep, 2), i in axes(essn_stacked_over_rep, 3)
]
println()
log_done("extracted stacked-over-repeat values", t_stage)

t_stage = log_step("preparing PCA samples")
modl2d_side = essn_2d_fmt |> f -> map(a -> a.modl2d, f) |>
                                  m ->
    map(a -> a[smwh_roi[2]+1+8:smwh_roi[2]+1+15, smwh_roi[1]+1-smw_ft:smwh_roi[1]+1+smw_ft], m) |>
    # d -> [permutedims(stack(@view d[i, j, :]), (3, 1, 2))
    #     for i in axes(d, 1), j in axes(d, 2)];
    d -> d
log_done("prepared PCA samples", t_stage)

t_stage = log_step("fitting PCA modes")
modes_pca_modl2d = [
    begin
        println("  [$tag] fitting PCA IB_idx=$c istp_idx=$i")
        flush(stdout)
        modl2d_side[c, :, :, i] |> m -> fit_pca_modes(n_pca_modes, m)
    end
    for c in axes(modl2d_side, 1), i in axes(modl2d_side, 4)
]
log_done("fit PCA modes", t_stage)

t_stage = log_step("analyzing per-shot trends")
trend_sidepeak_nvlp = [
    extr_fmt[c, r, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_sidepeak, selector_t_envelope, query_weight_kwargs)
    for c in axes(extr_fmt, 1), r in axes(extr_fmt, 2), i in axes(extr_fmt, 4)
]
log_done("analyzed per-shot trends", t_stage)

t_stage = log_step("analyzing stacked trends")
trend_stacked_over_rep = [
    extr_stacked_over_rep[c, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_sidepeak, selector_t_envelope, query_weight_kwargs)
    for c in axes(extr_stacked_over_rep, 1), i in axes(extr_stacked_over_rep, 3)
]
log_done("analyzed stacked trends", t_stage)
##  saving data, still problematic

# @save joinpath(path_output, @sprintf("%s_data.jld2", tag))
# val_vars
# essn_2d_fmt
# info_fmt
# essn_stacked_over_rep
# essn_stacked_over_rep_t
# fit_prfl_modl_over_rep_t_1d
# extr_fmt[1,1,1,1].sidepeak.fit_tailess
# extr_stacked_over_rep
# modl2d_side
# modes_pca_modl2d

## Overall plots
# for (c, IB) in enumerate(val_vars.IB)
#     tag_IB = gen_run_tag(get_bind_runinfo(runinfo, val_vars, c))
#     runinfo_plot = get_bind_runinfo(runinfo, val_vars, c)

#     t_stage = log_step("building trend figures for $tag_IB")
#     fig_trend, axs_trend = set_axis_sidepeak_nvlp!(n_dim_vars_per_IB, set_panel_trend_sidepeak_nvlp!, runinfo_plot)
#     fig_nvlp, axs_nvlp = set_axis_stack_all!(n_dim_vars_per_IB, set_panel_trend_nvlp!, runinfo_plot)
#     log_done("built trend figures for $tag_IB", t_stage)
#     for i in 1:n_istp
#         t_plot_stage = log_step("plotting and saving trends for $tag_IB istp=$(val_vars.istp[i])")
#         trend = trend_sidepeak_nvlp[c, :, i]
#         trend_stacked = trend_stacked_over_rep[c, i]
#         val_istp = val_vars.istp[i]
#         plot_trend_all!(axs_trend, trend, trend_stacked, val_istp)
#         plot_trend_nvlp!(axs_nvlp, trend, trend_stacked, val_istp)
#         resize_to_layout!(fig_trend)
#         resize_to_layout!(fig_nvlp)
#         for format in ["pdf", "png"]
#             fig_trend |> f -> save(joinpath(path_output, @sprintf("%s_%s_trend.%s", tag_IB, val_istp, format)), f; backend=CairoMakie)
#             fig_nvlp |> f -> save(joinpath(path_output, @sprintf("%s_%s_trend_nvlp.%s", tag_IB, val_istp, format)), f; backend=CairoMakie)
#         end
#         log_done("saved trends for $tag_IB istp=$(val_vars.istp[i])", t_plot_stage)
#     end

#     t_stage = log_step("building and saving PCA figure for $tag_IB")
#     fig_pca, axs_pca = set_axis_pca_dual_4x2!()
#     for idx_mode in 1:n_pca_modes, idx_istp in 1:n_istp
#         plot_mode_evol_freq_solo!(axs_pca[idx_istp, idx_mode], modes_pca_modl2d[c, idx_istp][idx_mode], val_vars.t_hold)
#     end
#     resize_to_layout!(fig_pca)
#     fig_pca |> f -> save(joinpath(path_output, @sprintf("%s_pca.pdf", tag_IB)), f; backend=CairoMakie)
#     log_done("saved PCA figure for $tag_IB", t_stage)
# end
# fig_trend |> display
##


## Large file generation for all shots

# fig_full, axs_solo, axs_stacked = set_axis_full(n_dim_vars_per_IB, set_panel_solo_essn_2d!)
# fig_full, axs_solo, axs_stacked = set_axis_full(n_dim_vars_per_IB, set_panel_solo_modl!)
# for c in 1:n_dim_vars[1], r in 1:n_dim_vars[2], t in 1:n_dim_vars[3], i in 1:n_dim_vars[4]
#     info = info_fmt[c, r, t, i]
#     print("\r\033[2Kplotting for runid $(info["runid"]), rep $r, $(info["t_hold"]) ms, $(info["istp"])")
#     draw_solo_modl!(axs_solo[r, t, i], extr_fmt[c, r, t, i], info)
#     # draw_solo_modl!(axs_live, extr_fmt[c, r, t, i], info)
# end
# println("Full axes ready: dimensions $(n_dim_vars_per_IB)")
# for c in 1:n_dim_vars[1], t in 1:n_dim_vars[3], i in 1:n_dim_vars[4]
#     info = info_fmt[c, 1, t, i] |> d -> merge(d, Dict("repeat" => "stacked"))
#     print("\r\033[2Kplotting for stacked runid $(info["runid"]), $(info["t_hold"]) ms, $(info["istp"])")
#     draw_solo_modl!(axs_stacked[t, i], extr_stacked_over_rep[c, t, i], info)
#     # draw_solo_modl!(axs_live, extr_stacked_over_rep[c, t, i], info)
# end
# println("Full modulation table drawn.")
# resize_to_layout!(fig_full)

# fig_full |> f -> save(joinpath(path_output, @sprintf("%s_essn_table.pdf", tag)), f; backend=CairoMakie)
# println("Full modulation plot saved.")
