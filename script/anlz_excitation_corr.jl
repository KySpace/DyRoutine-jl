t_stage = log_step("fitting PCA modes")
selector_t_pca_dens = @isdefined(selector_t_pca_dens) ? selector_t_pca_dens : (@isdefined(selector_t_pca) ? selector_t_pca : (t -> trues(length(t))))
selector_t_pca_modl = @isdefined(selector_t_pca_modl) ? selector_t_pca_modl : selector_t_pca_dens
n_pca_modes_prfl_modl = @isdefined(n_pca_modes_prfl_modl) ? n_pca_modes_prfl_modl : 12
modes_pca_dens2d = [
    begin
        println("  [$tag] fitting PCA IB_idx=$c")
        flush(stdout)
        essn_2d_fmt[c, :, :, :] |> es -> map(a -> a.dens2d_core |> filter_core_pca, es) |> es -> eachslice(es; dims=(1, 2)) |> m -> fit_pca_modes(n_pca_modes, m)
    end
    for c in axes(essn_2d_fmt, 1)
]
pca_spectra = [[
    begin
        mode = modes_pca_dens2d[c][m]
        spectral_weight = calc_spct_rep_evol(eachslice(mode.weight; dims=1), val_vars.t_hold, freq_query_pca; sel_evol=selector_t_pca_dens)
        peaks_prominent = spectral_weight.spct_mean_mask |> spct -> get_spectrum_peaks(freq_query_pca, spct; min_prom=0.2)
        (; spectral_weight, peaks_prominent)
    end
    for m in 1:n_pca_modes
] for c in axes(modes_pca_dens2d, 1)
]
log_done("fit density PCA modes", t_stage)

t_stage = log_step("analyzing per-shot trends")
trend_sidepeak_nvlp = [
    extr_fmt[c, r, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_spectrum, query_weight_kwargs)
    for c in axes(extr_fmt, 1), r in axes(extr_fmt, 2), i in axes(extr_fmt, 4)
]
log_done("analyzed per-shot trends", t_stage)

t_stage = log_step("analyzing stacked trends")
trend_extr_stacked_over_rep = [
    extr_stacked_over_rep[c, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_spectrum, query_weight_kwargs)
    for c in axes(extr_stacked_over_rep, 1), i in axes(extr_stacked_over_rep, 3)
]
trend_stacked_over_rep = [
    trend_sidepeak_nvlp[c, :, i] |> mean_dict
    for c in axes(trend_sidepeak_nvlp, 1), i in axes(trend_sidepeak_nvlp, 3)
]
log_done("analyzed stacked trends", t_stage)

t_stage = log_step("composing FT sidepeak profile evolution")
prfl_evol = [
    [
        extr_fmt[c, r, t, i].sidepeak.prfl_norm_tailess_px
        for t in axes(extr_fmt, 3)
    ] |> prfls -> reduce(hcat, prfls)
    for c in axes(extr_fmt, 1), r in axes(extr_fmt, 2), i in axes(extr_fmt, 4)
]
prfl_evol_stacked = [
    [
        extr_fmt[c, r, t, i].sidepeak.prfl_norm_tailess_px
        for r in axes(extr_fmt, 2), t in axes(extr_fmt, 3)
    ] |> prfls -> mean(prfls; dims=1) |> vec |> prfls -> reduce(hcat, prfls)
    for c in axes(extr_fmt, 1), i in axes(extr_fmt, 4)
]
log_done("finished composing FT sidepeak profile evolution", t_stage)

t_stage = log_step("fitting modulation profile PCA modes")
mask_y_modl_pca = (0.06 .<= y_modl .<= 0.6)
any(mask_y_modl_pca) || throw(ArgumentError("modulation profile PCA wavenum selector 0.06-0.6 selected no y_modl values."))
y_modl_pca = y_modl[mask_y_modl_pca]
idx_y_modl_pca = findall(mask_y_modl_pca)
samples_pca_prfl_modl = [
    [
        prfl_evol[c, r, i][idx_y, t]
        for i in axes(prfl_evol, 3), idx_y in idx_y_modl_pca
    ]
    for c in axes(prfl_evol, 1), t in axes(first(prfl_evol), 2), r in axes(prfl_evol, 2)
]
modes_pca_prfl_modl = fit_pca_modes(n_pca_modes_prfl_modl, samples_pca_prfl_modl)
pca_spectra_prfl_modl = [
    begin
        mode = modes_pca_prfl_modl[m]
        spectral_weight = calc_spct_rep_evol(
            [vec(mode.weight[c, :, r]) for r in axes(mode.weight, 3)],
            val_vars.t_hold,
            freq_query_pca_modl;
            sel_evol=selector_t_pca_modl,
        )
        peaks_prominent = spectral_weight.spct_mean_mask |> spct -> get_spectrum_peaks(freq_query_pca_modl, spct; min_prom=0.2)
        (; spectral_weight, peaks_prominent)
    end
    for m in 1:n_pca_modes_prfl_modl, c in axes(samples_pca_prfl_modl, 1)
]
log_done("fit modulation profile PCA modes", t_stage)

meta_corr = merge(
    meta_extr,
    (;
        kind="excitation_corr",
        path_output,
        n_pca_modes_prfl_modl,
        y_modl_pca,
        freq_query_pca_modl,
        selector_t_pca_dens_val=val_vars.t_hold[selector_t_pca_dens(val_vars.t_hold)],
        selector_t_pca_modl_val=val_vars.t_hold[selector_t_pca_modl(val_vars.t_hold)],
        selector_t_spectrum_val=NamedTuple{propertynames(selector_t_spectrum)}(
            Tuple(selector_t_spectrum[key](val_vars.t_hold) |> mask -> val_vars.t_hold[mask] for key in propertynames(selector_t_spectrum))
        ),
    ),
)

t_stage = log_step("saving excitation correlation cache")
path_cache_corr = joinpath(path_output, @sprintf("%s_corr.jld2", tag))
JLD2.jldsave(
    path_cache_corr;
    meta_corr,
    trend_sidepeak_nvlp,
    trend_extr_stacked_over_rep,
    trend_stacked_over_rep,
    prfl_evol,
    prfl_evol_stacked,
    modes_pca_dens2d,
    pca_spectra,
    modes_pca_prfl_modl,
    pca_spectra_prfl_modl,
)
log_done("saved excitation correlation cache", t_stage)
