using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using CairoMakie: Axis
using Colors: Oklch
using LaTeXStrings

function set_panel_trend_sidepeak_nvlp!(gl::GridLayout, col::Int; extra=false)
    gl |> clean_gridlayout!
    w, h = (400, 200)
    col_freq = extra ? 3 : 2
    ax_evol_sum_dens = Axis(gl[1, 1]; width=w, height=h, ylabel="density sum")
    ax_evol_weight = Axis(gl[2, 1]; width=w, height=h, ylabel="side peak \nweight")
    ax_evol_height = Axis(gl[3, 1]; width=w, height=h, ylabel="side peak \nheight")
    ax_evol_width = Axis(gl[4, 1]; width=w, height=h, ylabel="side peak \nwidth (μm⁻¹)")
    ax_evol_wavenum = Axis(gl[5, 1]; width=w, height=h, ylabel="side peak \nwavenum (μm⁻¹)")
    ax_evol_sizes = Axis(gl[6, 1]; width=w, height=h, ylabel="envelope size (μm)")
    if extra
        ax_evol_extra_weight = Axis(gl[2, 2]; width=w, height=h, ylabel="side peak \nweight")
        ax_evol_extra_height = Axis(gl[3, 2]; width=w, height=h, ylabel="side peak \nheight")
        ax_evol_extra_width = Axis(gl[4, 2]; width=w, height=h, ylabel="side peak \nwidth (μm⁻¹)")
        ax_evol_extra_wavenum = Axis(gl[5, 2]; width=w, height=h, ylabel="side peak \nwavenum (μm⁻¹)")
    end
    ax_freq_weight = Axis(gl[2, col_freq]; width=w, height=h)
    ax_freq_height = Axis(gl[3, col_freq]; width=w, height=h)
    ax_freq_width = Axis(gl[4, col_freq]; width=w, height=h)
    ax_freq_wavenum = Axis(gl[5, col_freq]; width=w, height=h)
    ax_freq_sizes = Axis(gl[6, col_freq]; width=w, height=h)
    rowgap!(gl, 4)
    colgap!(gl, 4)
    dict_axs = Dict(
        "evol-dens-sum" => ax_evol_sum_dens,
        "evol-weight" => ax_evol_weight,
        "evol-height" => ax_evol_height,
        "evol-width" => ax_evol_width,
        "evol-wavenum" => ax_evol_wavenum,
        "evol-sizes" => ax_evol_sizes,
        "freq-weight" => ax_freq_weight,
        "freq-height" => ax_freq_height,
        "freq-width" => ax_freq_width,
        "freq-wavenum" => ax_freq_wavenum,
        "freq-sizes" => ax_freq_sizes,
    )
    if extra
        dict_axs["evol-extra-weight"] = ax_evol_extra_weight
        dict_axs["evol-extra-height"] = ax_evol_extra_height
        dict_axs["evol-extra-width"] = ax_evol_extra_width
        dict_axs["evol-extra-wavenum"] = ax_evol_extra_wavenum
    end
    for ax in values(dict_axs)
        hideydecorations!(ax; label=true, ticklabels=false, ticks=false, grid=false, minorticks=false, minorgrid=false)
        hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
        if col == 1
            ax.ylabelvisible = true
        end
    end
    ax_evol_sizes.xticklabelsvisible = true
    ax_freq_sizes.xticklabelsvisible = true
    ax_evol_sizes.xlabelvisible = true
    ax_freq_sizes.xlabelvisible = true
    ax_evol_sizes.xlabel = "t hold (ms)"
    ax_freq_sizes.xlabel = "freq (Hz)"
    return dict_axs
end

function plot_mode_evol_freq_duet!(axs::Dict{String,Axis}, mode::ModeWeight, val_t::AbstractVector)
    ndims(mode.profile) == 3 && size(mode.profile, 1) == 2 || throw(ArgumentError("mode.profile must be a 3D array with size[1]==2."))
    clrmap = gen_clrmap_posneg(0.60 * 360, 0.96 * 360)
    c = maximum(abs, mode.profile)
    heatmap!(axs["l"], mode.profile[1, :, :]'; colormap=clrmap, colorrange=(-c, c))
    heatmap!(axs["r"], mode.profile[2, :, :]'; colormap=clrmap, colorrange=(-c, c))
    axs["l"].aspect = DataAspect()
    axs["r"].aspect = DataAspect()
    axs["l"] |> hidedecorations!
    axs["r"] |> hidedecorations!
    for rep = 1:size(mode.weight, 1)
        lines!(axs["evol"], val_t, mode.weight[rep, :]; color=(:black, 0.2))
    end
end

function plot_mode_evol_freq_solo!(axs::Dict{String,Axis}, mode::ModeWeight, val_t::AbstractVector)
    ndims(mode.profile) == 2 || throw(ArgumentError("mode.profile must be a 2D array. "))
    clrmap = gen_clrmap_posneg(0.60 * 360, 0.96 * 360)
    c = maximum(abs, mode.profile)
    heatmap!(axs["mode"], mode.profile[:, :]; colormap=clrmap, colorrange=(-c, c))
    axs["mode"].aspect = DataAspect()
    axs["mode"] |> hidedecorations!
    for rep = 1:size(mode.weight, 1)
        lines!(axs["evol"], val_t, mode.weight[rep, :]; color=(:black, 0.2))
    end
end

function plot_trends!(axs::Dict, trend::Dict, istp; to_clean=false, alpha=1.0, is_stacked=false, to_legend=false, to_overlay=false)
    hue_theme = hue_theme_istp[istp]
    clr_mmt = Oklch(0.52, 0.14, hue_theme)
    clr_fit = (:springgreen3, 1.0)
    clr_theme = Oklch(0.52, 0.14, hue_theme)
    clr_theme1 = Oklch(0.52, 0.14, hue_theme - 20)
    clr_theme2 = Oklch(0.52, 0.14, hue_theme + 20)
    clr_shade_selected = RGBAf(Oklch(0.95, 0.1, hue_theme), 0.2)
    if to_clean
        for (k, obj) in axs
            obj isa Axis && empty!(obj)
        end
        vspan!(axs["evol-weight"], trend["t_vec_sel_sp"][1] - 1, trend["t_vec_sel_sp"][end] + 1; color=clr_shade_selected)
        vspan!(axs["evol-height"], trend["t_vec_sel_sp"][1] - 1, trend["t_vec_sel_sp"][end] + 1; color=clr_shade_selected)
        vspan!(axs["evol-width"], trend["t_vec_sel_sp"][1] - 1, trend["t_vec_sel_sp"][end] + 1; color=clr_shade_selected)
        vspan!(axs["evol-wavenum"], trend["t_vec_sel_sp"][1] - 1, trend["t_vec_sel_sp"][end] + 1; color=clr_shade_selected)
        vspan!(axs["evol-sizes"], trend["t_vec_sel_nvlp"][1] - 1, trend["t_vec_sel_nvlp"][end] + 1; color=clr_shade_selected)
    end
    lines!(axs["evol-dens-sum"], trend["t_vec"], trend["evol-all-dens-sum"]; color=(clr_theme, alpha), label="sum")
    lines!(axs["evol-weight"], trend["t_vec"], trend["evol-all-fit-weight"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-height"], trend["t_vec"], trend["evol-all-fit-height"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-width"], trend["t_vec"], trend["evol-all-fit-width"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-wavenum"], trend["t_vec"], trend["evol-all-fit-wavenum"]; color=(clr_fit, alpha), label="fit")
    if to_overlay
        lines!(axs["evol-extra-weight"], trend["t_vec"], trend["evol-all-moment-weight"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-extra-height"], trend["t_vec"], trend["evol-all-moment-height"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-extra-width"], trend["t_vec"], trend["evol-all-moment-width"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-extra-wavenum"], trend["t_vec"], trend["evol-all-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    else
        lines!(axs["evol-weight"], trend["t_vec"], trend["evol-all-moment-weight"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-height"], trend["t_vec"], trend["evol-all-moment-height"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-width"], trend["t_vec"], trend["evol-all-moment-width"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-wavenum"], trend["t_vec"], trend["evol-all-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    end
    lines!(axs["evol-sizes"], trend["t_vec"], trend["evol-all-fit-size-x"]; color=(clr_theme1, alpha), label="fit")
    lines!(axs["evol-sizes"], trend["t_vec"], trend["evol-all-fit-size-y"]; color=(clr_theme2, alpha), label="fit")
    lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-fit-weight"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-fit-height"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-fit-width"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-fit-wavenum"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-moment-weight"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-moment-height"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-moment-width"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-sizes"], trend["freq_query"], trend["freq-sel-fit-size-x"]; color=(clr_theme1, alpha), label="fit size x")
    lines!(axs["freq-sizes"], trend["freq_query"], trend["freq-sel-fit-size-y"]; color=(clr_theme2, alpha), label="fit size y")
    ylims!(axs["evol-weight"], -0.02, 0.22)
    ylims!(axs["evol-height"], -0.1, 1.1)
    ylims!(axs["evol-width"], 0.02, 0.205)
    ylims!(axs["evol-wavenum"], 0.22, 0.38)
    ylims!(axs["evol-sizes"], 1, 11)
    if to_legend
        axislegend(axs["evol-dens-sum"]; position=:rt, framevisible=false, labelsize=14)
        axislegend(axs["freq-weight"]; position=:lt, framevisible=false, labelsize=14)
        axislegend(axs["freq-sizes"]; position=:lt, framevisible=false, labelsize=14)
    end
    if to_overlay
        ylims!(axs["evol-extra-weight"], -0.02, 0.22)
        ylims!(axs["evol-extra-height"], -0.1, 1.1)
        ylims!(axs["evol-extra-width"], 0.02, 0.205)
        ylims!(axs["evol-extra-wavenum"], 0.22, 0.38)
    end
end

function plot_trend_all!(axs_trend::Dict, trend_reps::AbstractVector, trend_stacked_over_rep::Dict, istp)
    function set_tick_grid!(axs)
        for ax in [axs["evol-weight"], axs["evol-height"], axs["evol-width"], axs["evol-wavenum"], axs["evol-sizes"]]
            ax.xticks = 0:50:200
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(5)
        end
        for ax in [axs["freq-weight"], axs["freq-height"], axs["freq-width"], axs["freq-wavenum"], axs["freq-sizes"]]
            ax.xticks = 0:10:100
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(2)
        end
    end
    for r = axes(trend_reps, 1)
        trend = trend_reps[r]
        # plot on both the individual reps and all reps combined
        for (a, axs) in enumerate([axs_trend["repeats"][r], axs_trend["all"]])
            alpha = a == 1 ? 1.0 : 0.5
            # shade the selected time points on all reps combined only once
            to_clean = r == 1 || a == 1
            to_legend = a == 1
            to_overlay = a > 1
            axs |> set_tick_grid!
            plot_trends!(axs, trend, istp; to_clean, alpha, to_legend, to_overlay)
        end
    end
    # plot on the stacked axes
    axs = axs_trend["stacked"]
    plot_trends!(axs, trend_stacked_over_rep, istp; to_clean=true, alpha=1.0, to_legend=true)
    axs |> set_tick_grid!
end
