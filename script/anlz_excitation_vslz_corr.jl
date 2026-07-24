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

function save_prfl_comparison!(rows, title, name, tag_IB, config, val_t, val_istp)
    width = config.width_to_time * (maximum(val_t) - minimum(val_t))
    fig, axs, colorbars = set_axis_prfl_comparison!(
        rows,
        val_istp,
        title;
        width,
        height=config.height,
    )
    plot_prfl_comparison!(
        fig,
        axs,
        colorbars,
        rows,
        val_t,
        val_istp;
        ylims=config.ylims,
        colorrange=config.colorrange,
    )
    for format in ("svg", "png")
        save(joinpath(path_output, @sprintf("%s_[%s].%s", name, tag_IB, format)), fig; backend=CairoMakie)
    end
end

for (c, tag_IB) in enumerate(tag_IBs)
    local t_stage = log_step("profile comparison figures for $tag_IB")
    runinfo_plot = runinfo_plots[c]
    val_t = val_vars.t_hold
    val_istp = val_vars.istp
    modl_stack = reshape([prfl_evol_stacked[c, i] for i in axes(prfl_evol_stacked, 2)], 1, :)
    modl_rows_stack = [(label="stack", evol=modl_stack, pos=y_modl)]
    modl_rows_reps = [
        (
            label="rep $(val_vars.rep[r])",
            evol=reshape([prfl_evol[c, r, i] for i in axes(prfl_evol, 3)], 1, :),
            pos=y_modl,
        )
        for r in axes(prfl_evol, 2)
    ]
    push!(modl_rows_reps, modl_rows_stack[1])
    save_prfl_comparison!(modl_rows_stack, "$tag_IB modulation profile stack", "prfl_modl_stack", tag_IB, vis_evol_prfl_modl, val_t, val_istp)
    save_prfl_comparison!(modl_rows_reps, "$tag_IB modulation profile reps and stack", "prfl_modl_reps_stack", tag_IB, vis_evol_prfl_modl, val_t, val_istp)

    essn_ref = essn_2d_fmt[c, 1, 1, 1]
    pos_axial = ((-essn_ref.smwh_core[2]):essn_ref.smwh_core[2]) .* essn_ref.step_posi[2]
    pos_radial = ((-essn_ref.smwh_core[1]):essn_ref.smwh_core[1]) .* essn_ref.step_posi[1]
    height_radial = vis_evol_prfl_core.height * (maximum(pos_radial) - minimum(pos_radial)) / (maximum(pos_axial) - minimum(pos_axial))
    core_stack = [
        (label="axial stack", evol=reshape([prfl_axial_evol_stacked[c, i] for i in axes(prfl_axial_evol_stacked, 2)], 1, :), pos=pos_axial, height=vis_evol_prfl_core.height),
        (label="radial stack", evol=reshape([prfl_radial_evol_stacked[c, i] for i in axes(prfl_radial_evol_stacked, 2)], 1, :), pos=pos_radial, height=height_radial),
    ]
    core_rows_reps = Any[]
    for r in axes(prfl_axial_evol, 2)
        push!(core_rows_reps, (label="rep $(val_vars.rep[r]) axial", evol=reshape([prfl_axial_evol[c, r, i] for i in axes(prfl_axial_evol, 3)], 1, :), pos=pos_axial, height=vis_evol_prfl_core.height))
        push!(core_rows_reps, (label="rep $(val_vars.rep[r]) radial", evol=reshape([prfl_radial_evol[c, r, i] for i in axes(prfl_radial_evol, 3)], 1, :), pos=pos_radial, height=height_radial))
    end
    append!(core_rows_reps, core_stack)
    save_prfl_comparison!(core_stack, "$tag_IB axial/radial profile stack", "prfl_core_stack", tag_IB, vis_evol_prfl_core, val_t, val_istp)
    save_prfl_comparison!(core_rows_reps, "$tag_IB axial/radial profile reps and stack", "prfl_core_reps_stack", tag_IB, vis_evol_prfl_core, val_t, val_istp)
    log_done("finished profile comparison figures for $tag_IB", t_stage)

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
    path_pca = joinpath(path_output, "PCA density", tag_IB)
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

if !isnothing(modes_pca_prfl_modl)
    t_stage = log_step("building and saving modulation profile PCA figures")
    path_pca_prfl = joinpath(path_output, "PCA prfl modl", tag)
    isdir(path_pca_prfl) || mkpath(path_pca_prfl)
    fig_pca_prfl = Figure()
    for idx_mode in 1:n_pca_modes_prfl_modl
        mode = modes_pca_prfl_modl[idx_mode]
        spectra_params = [pca_spectra_prfl_modl[idx_mode, c] for c in axes(pca_spectra_prfl_modl, 2)]
        fig_pca_prfl.layout |> clean_gridlayout!
        gl_pca_prfl = GridLayout()
        fig_pca_prfl[1, 1] = gl_pca_prfl
        axs_pca_prfl = set_panel_pca_duet_params!(
            gl_pca_prfl,
            val_vars.IB;
            mode_kind=:profile1d,
            width_evol=400,
            width_spct=400,
            height_evol=120,
            height_spct=120,
        )
        gl_pca_prfl[0, 1] = Label(fig_pca_prfl, "$tag | profile PCA #$idx_mode"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
        plot_mode_evol_spct_duet_params!(
            axs_pca_prfl,
            mode,
            spectra_params,
            val_vars.IB,
            val_vars.istp;
            mode_kind=:profile1d,
            y_modl=y_modl_pca,
        )
        resize_to_layout!(fig_pca_prfl)
        fig_pca_prfl |> f -> save(joinpath(path_pca_prfl, @sprintf("%s_prfl_modl_%d.png", tag, idx_mode)), f; backend=CairoMakie)
    end
    log_done("saved modulation profile PCA figures", t_stage)
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

for spec in trend_property_specs
    for variant in spec.variants
        name_variant_file = replace(variant.label, r"[^A-Za-z0-9]+" => "_") |> s -> replace(s, r"^_|_$" => "")
        title_property = "$tag | $(spec.name) | $(variant.label)"
        fig_property, axs_IB_istp = set_axis_trend_variant_IB_istp!(
            val_vars.IB,
            val_vars.istp,
            spec,
            variant,
            title_property;
            trend_panel_per_prop_kwargs...,
        )
        for c in axes(trend_sidepeak_nvlp, 1), i in axes(trend_sidepeak_nvlp, 3)
            fit_evol = (!@isdefined(fit_evol_properties) || isnothing(fit_evol_properties)) ? nothing : get(fit_evol_properties[c, i], variant.name, nothing)
            plot_trend_variant_overlay!(
                axs_IB_istp[c, i],
                trend_sidepeak_nvlp[c, :, i],
                spec,
                variant,
                val_vars.istp[i];
                fit_evol,
                to_legend=(c == firstindex(trend_sidepeak_nvlp, 1) && i == firstindex(trend_sidepeak_nvlp, 3)),
            )
        end
        resize_to_layout!(fig_property)
        for format in ["pdf", "png"]
            fig_property |> f -> save(
                joinpath(path_output, @sprintf("trend_overlay_%s_%s_[%s].%s", spec.name, name_variant_file, tag, format)),
                f;
                backend=CairoMakie,
            )
        end
    end
end
