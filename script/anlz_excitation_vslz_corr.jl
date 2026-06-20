## correlation visualization
t_stage = log_step("building and saving spectrum-vs-IB figures")
trend_spectrum_groups = (;
    stacked=trend_extr_stacked_over_rep,
    all=trend_stacked_over_rep,
)
for spec in trend_property_specs
    fig_spectrum, axs_spectrum = set_axis_spectrum_property_IB_istp!(
        val_vars.IB,
        val_vars.istp,
        spec,
        "$tag | spectrum vs IB | $(spec.name)";
        groups=trend_spectrum_IB_groups,
        trend_spectrum_IB_kwargs...,
    )
    plot_spectrum_property_IB_istp!(
        axs_spectrum,
        trend_spectrum_groups,
        val_vars.IB,
        val_vars.istp,
        spec;
        trend_spectrum_IB_plot_kwargs...,
    )
    resize_to_layout!(fig_spectrum)
    for format in ["pdf", "png"]
        fig_spectrum |> f -> save(
            joinpath(path_output, @sprintf("bands_IB_%s_[%s].%s", spec.name, tag, format)),
            f;
            backend=CairoMakie,
        )
    end
end
log_done("saved spectrum-vs-IB figures", t_stage)

for (c, tag_IB) in enumerate(tag_IBs)
    runinfo_plot = runinfo_plots[c]

    local t_stage = log_step("sidepeak distribution evolution for $tag_IB")
    fig_prfl_evol, axs_prfl_evol = set_axis_prfl_modl_evol!(
        val_vars.rep,
        val_vars.istp,
        "$tag_IB modulation sidepeak profile";
        width=30 * size(prfl_evol[c, 1, 1], 2),
        height=400,
    )
    plot_prfl_modl_evol!(
        axs_prfl_evol,
        prfl_evol[c, :, :],
        prfl_evol_stacked[c, :],
        val_vars.t_hold,
        y_modl,
        val_vars.istp,
    )
    fig_prfl_evol |> resize_to_layout!
    for format in ["svg", "png"]
        fig_prfl_evol |> f -> save(joinpath(path_output, @sprintf("prfl_modl_evol_[%s].%s", tag_IB, format)), f; backend=CairoMakie)
    end
    log_done("finished sidepeak distribution for $tag_IB", t_stage)

    local t_stage = log_step("building trend figures for $tag_IB")
    panel_setter = (gl, col; extra=false) -> set_panel_trend_properties!(
        gl,
        trend_property_specs;
        col,
        extra,
        trend_panel_per_IB_kwargs...,
    )
    fig_trend, axs_trend = set_axis_sidepeak_nvlp!(n_dim_vars_per_IB, panel_setter, runinfo_plot)
    log_done("built trend figures for $tag_IB", t_stage)
    for i in 1:n_istp
        val_istp = val_vars.istp[i]
        t_plot_stage = log_step("plotting and saving trends for $tag_IB istp=$val_istp")
        trend_reps = trend_sidepeak_nvlp[c, :, i]
        trend_stacked = trend_extr_stacked_over_rep[c, i]
        plot_trend_all!(axs_trend, trend_reps, trend_stacked, val_istp; property_specs=trend_property_specs)
        resize_to_layout!(fig_trend)
        for format in ["pdf", "png"]
            fig_trend |> f -> save(joinpath(path_output, @sprintf("property_trends_[%s.%s].%s", tag_IB, val_istp, format)), f; backend=CairoMakie)
        end
        log_done("saved trends for $tag_IB istp=$val_istp", t_plot_stage)
    end

    t_stage = log_step("building and saving PCA figure for $tag_IB")
    path_pca = joinpath(path_output, "PCA modes", tag_IB)
    isdir(path_pca) || mkpath(path_pca)
    fig_pca_mode = Figure()
    for idx_mode in 1:n_pca_modes
        mode = modes_pca_dens2d[c][idx_mode]
        spct_weight, peaks = pca_spectra[c][idx_mode]
        fig_pca_mode.layout |> clean_gridlayout!
        gl_pca_mode = GridLayout()
        fig_pca_mode[1, 1] = gl_pca_mode
        axs_pca_mode = set_panel_pca_duet!(gl_pca_mode)
        gl_pca_mode[0, 1] = Label(fig_pca_mode, "$tag_IB | #$idx_mode"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
        plot_mode_evol_spct_duet!(axs_pca_mode, mode, spct_weight, peaks, val_vars.istp; step_posi=px_in_um, smwh=smwh_core)
        gl_pca_mode |> l -> rowgap!(l, 0)
        resize_to_layout!(fig_pca_mode)
        fig_pca_mode |> f -> save(joinpath(path_pca, @sprintf("%s_%d.png", tag_IB, idx_mode)), f; backend=CairoMakie)
    end
    log_done("saved PCA figure for $tag_IB", t_stage)
end

for spec in trend_property_specs
    title_property = "$tag | $(spec.name)"
    fig_property, axs_IB_istp = set_axis_trend_property_IB_istp!(
        val_vars.IB,
        val_vars.istp,
        n_rep,
        spec,
        title_property;
        groups=trend_all_IB_groups,
        trend_panel_per_prop_kwargs...,
    )
    for c in axes(trend_sidepeak_nvlp, 1), i in axes(trend_sidepeak_nvlp, 3)
        plot_trend_all!(
            axs_IB_istp[c, i],
            trend_sidepeak_nvlp[c, :, i],
            trend_extr_stacked_over_rep[c, i],
            val_vars.istp[i];
            property_specs=[spec],
        )
    end
    resize_to_layout!(fig_property)
    for format in ["pdf", "png"]
        fig_property |> f -> save(
            joinpath(path_output, @sprintf("trend_%s_[%s].%s", spec.name, tag, format)),
            f;
            backend=CairoMakie,
        )
    end
end
