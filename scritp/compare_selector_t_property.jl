## Run this after data is generated
for t_min in 0:10:100, t_max in 20:20:240
    if t_max - t_min <= 30
        continue
    end
    tag_range = @sprintf("%d-%d", t_min, t_max)
    selector_t_common = t -> t_min .< t .< t_max
    selector_t_spectrum = (;
        number=selector_t_common,
        sp_weight=selector_t_common,
        sp_height=selector_t_common,
        sp_width=selector_t_common,
        sp_wavenum=selector_t_common,
        nvlp_size=selector_t_common,
    )
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
    log_done("analyzed stacked trends", t_stage)

    trend_stacked_over_rep = [
        trend_sidepeak_nvlp[c, :, i] |> mean_dict
        for c in axes(trend_sidepeak_nvlp, 1), i in axes(trend_sidepeak_nvlp, 3)
    ]

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
        fig_spectrum |> f -> save(
            joinpath(path_output, tag, @sprintf("%s_%s_spectrum_IB.png", spec.name, tag_range)),
            f;
            backend=CairoMakie,
            force=true,
        )
    end
    log_done("saved spectrum-vs-IB figures", t_stage)
end
