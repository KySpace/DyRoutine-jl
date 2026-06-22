log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))

get_bind_date(runinfo, idx_bind) = hasproperty(runinfo, :date_runid) ? first(runinfo.date_runid[idx_bind]) : runinfo.date
get_bind_runid(runinfo, idx_bind) =
    if hasproperty(runinfo, :date_runid)
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
fmt_dens = if @isdefined(format_dens_runinfo_kwargs)
    local formatter = @isdefined(format_dens_runinfo_fn) ? format_dens_runinfo_fn : format_dens_runinfo
    formatter(runinfo; format_dens_runinfo_kwargs...)
else
    format_dens_runinfo(runinfo; path_root, year_test, wh_corner, smwh_roi, len_avg_peak, sel_vars)
end
runinfo = fmt_dens.runinfo
(; val_vars, dens_full_fmt, wh_dens, xy_peak_px, n_dim_vars, name_dims) = fmt_dens
n_variation = prod(n_dim_vars)
n_dim_vars_per_IB = Tuple(n_dim_vars[2:end])
n_IB, n_rep, n_main, n_istp = n_dim_vars
log_done("formatted density data: axes $(name_dims) dims $(n_dim_vars), per-IB dims $(n_dim_vars_per_IB), image $(size(first(dens_full_fmt)))", t_stage)

xy_peak_core_per_IB_rep = if @isdefined(xy_peak_core_fixed)
    fill(xy_peak_core_fixed, n_IB, n_rep)
else
    [
        dens_full_fmt[c, r, :, :] |> mean |> ds ->
            find_positive_cluster_center(ds, smwh_core; len_avg=len_avg_peak) |> cent -> round.(Int, cent)
        for c in 1:n_IB, r in 1:n_rep
    ]
end
xy_peak_core = xy_peak_core_per_IB_rep |> p -> repeat(p, inner=ntuple(i -> i in (3, 4) ? n_dim_vars[i] : 1, length(n_dim_vars)))

t_stage = log_step("calculating solo essentials for $(length(dens_full_fmt)) shots")
essn_2d_fmt = map(
    (d, xy) -> calc_solo_essn_2d(d, smwh_roi .+ 1, smwh_roi, px_in_um, xy, smwh_core; smwh_strip=smwh_core, mask_modl),
    dens_full_fmt,
    xy_peak_core,
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
        essns_r = [essn_2d_fmt[c, r, t, i] for r in axes(essn_2d_fmt, 2)] |> vec
        print("\r  [$tag] stacking over rep IB_idx=$c t_hold=$t istp_idx=$i n=$(length(essns_r))")
        flush(stdout)
        calc_stacked_essn(essns_r)
    end
    for c in axes(essn_2d_fmt, 1), t in axes(essn_2d_fmt, 3), i in axes(essn_2d_fmt, 4)
]
println()
log_done("stacked essentials over rep only", t_stage)

t_stage = log_step("fitting stacked modulation tails")
fit_prfl_modl_over_rep_1d = [
    essn_stacked_over_rep[c, t, i] |>
    e -> fit_prfl_modl_sidepeak_decay_1d(y_modl, e.prfl_modl.side.normed_px, selector_tail_sidepeak(y_modl); fit_stack_kwargs...)
    for c in axes(essn_2d_fmt, 1), t in axes(essn_2d_fmt, 3), i in axes(essn_2d_fmt, 4)
]
log_done("fit stacked modulation tails", t_stage)

t_stage = log_step("extracting per-shot sidepeak/envelope values")
extr_fmt = [
    begin
        print("\r  [$tag] extracting shots IB_idx=$c rep-$r t_idx=$t istp_idx=$i")
        flush(stdout)
        essn_2d_fmt[c, r, t, i] |> e -> calc_solo_extr(
            e,
            fit_prfl_modl_over_rep_1d[c, t, i];
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
        print("\r  [$tag] extracting stacked IB_idx=$c t_idx=$t istp_idx=$i")
        flush(stdout)
        essn_stacked_over_rep[c, t, i] |> e -> calc_solo_extr(
            e,
            fit_prfl_modl_over_rep_1d[c, t, i];
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

tag_IBs = [gen_run_tag(get_bind_runinfo(runinfo, val_vars, c)) for c in axes(extr_fmt, 1)]
runinfo_plots = [get_bind_runinfo(runinfo, val_vars, c) for c in axes(extr_fmt, 1)]
selected_values(vals, selector) = vals[selector(vals)]
selected_range(vals, selector) = selected_values(vals, selector) |> extrema
info_fitting = (;
    stack_tail=(;
        fit_function="fit_prfl_modl_sidepeak_decay_1d",
        model_function="fit_prfl_modl_sidepeak_decay_1d_model",
        tail_function="fit_prfl_modl_sidepeak_decay_1d_tail",
        variable="y_modl",
        val_fit=selected_values(y_modl, selector_tail_sidepeak),
        range_fit=selected_range(y_modl, selector_tail_sidepeak),
        range_full=extrema(y_modl),
    ),
    tailess=(;
        fit_function="fit_prfl_modl_sidepeak_1d",
        model_function="fit_prfl_modl_sidepeak_1d_model",
        variable="y_modl",
        val_fit=selected_values(y_modl, selector_sidepeak),
        range_fit=selected_range(y_modl, selector_sidepeak),
        range_full=extrema(y_modl),
    ),
    moment=(;
        variable="y_modl",
        val_fit=selected_values(y_modl, selector_moment),
        range_fit=selected_range(y_modl, selector_moment),
        range_full=extrema(y_modl),
    ),
    envelope_asymm=(;
        fit_function="fit_dens2d_gaussian_elliptic_disk",
        model_function="fit_dens2d_gaussian_elliptic_disk_model",
        variables=("x_posi", "y_posi"),
        range_x=extrema(x_posi),
        range_y=extrema(y_posi),
    ),
    envelope_round=(;
        fit_function="fit_dens2d_gaussian_round_disk",
        model_function="fit_dens2d_gaussian_round_disk_model",
        variables=("x_posi", "y_posi"),
        range_x=extrema(x_posi),
        range_y=extrema(y_posi),
    ),
)

meta_extr = (;
    kind="excitation_extr",
    tag,
    path_output,
    runinfo,
    val_vars,
    name_dims,
    n_dim_vars,
    n_dim_vars_per_IB,
    n_IB,
    n_rep,
    n_main,
    n_istp,
    n_pca_modes,
    fmt_dens,
    px_in_um,
    smwh_roi,
    smwh_core,
    y_modl,
    x_posi,
    y_posi,
    freq_query,
    freq_query_pca,
    tag_IBs,
    runinfo_plots,
    info_fmt,
    info_fitting,
    trend_property_specs,
    trend_panel_per_IB_kwargs,
    trend_panel_per_prop_kwargs,
    trend_all_IB_groups,
    trend_spectrum_IB_groups,
    trend_spectrum_IB_kwargs,
    trend_spectrum_IB_plot_kwargs,
)

t_stage = log_step("saving excitation extraction cache")
path_cache_extr = joinpath(path_output, @sprintf("%s_essn_extr.jld2", tag))
JLD2.jldsave(
    path_cache_extr;
    meta_extr,
    essn_2d_fmt,
    essn_stacked_over_rep,
    fit_prfl_modl_over_rep_1d,
    extr_fmt,
    extr_stacked_over_rep,
)
log_done("saved excitation extraction cache", t_stage)
