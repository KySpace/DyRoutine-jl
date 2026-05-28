using GLMakie: Grid
using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using CairoMakie: Axis
using Colors: Oklch
using LaTeXStrings

function gen_clrmap_posneg(hue_pos, hue_neg)
    return [Oklch(1 - abs(t), 0.24 * abs(t), t > 0 ? hue_pos : hue_neg) |> c -> RGBAf(c) for t in range(-1, 1; length=256)]
end

function plot_num_stat_evo!(
    ax::Axis,
    val_t::AbstractVector,
    stat_number::AbstractVector,
    val_istp;
    hue_shift=0.0,
    label=nothing
)
    hue = hue_theme_istp[val_istp] + hue_shift
    val_mean = map(s -> s[1], stat_number)
    val_err = map(s -> s[2], stat_number)
    clr_line = Oklch(0.62, 0.18, hue) |> c -> RGBAf(c, 0.95)
    clr_bar = Oklch(0.62, 0.18, hue) |> c -> RGBAf(c, 0.40)

    errorbars!(ax, val_t, val_mean, val_err; color=clr_bar)
    lines!(ax, val_t, val_mean; color=clr_line, linewidth=2.6, label=val_istp)
    scatter!(ax, val_t, val_mean; color=clr_bar, markersize=9)
    # return path_plot
end

function set_axis_sidepeak_nvlp!(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer}, panel_setter::Function, runinfo)
    length(n_dim_vars) == 3 || throw(ArgumentError("n_dim_vars must be a 3-tuple"))
    fig = Figure()
    axs_repeats = Array{Dict}(undef, n_dim_vars[1])
    fig[0, 1] = Label(fig, text="$(runinfo.date) $(@sprintf("run%02d", runinfo.runid)) IB=$(@sprintf("%.3f", runinfo.IB))A $(runinfo.tag_head)"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    for r in 1:n_dim_vars[1]
        fig[1, r] = Label(fig, text="repeat $r"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
        print("\r\033[2K\rbuilding axes for side peak trend for repeat $r")
        gl = GridLayout()
        fig[2, r] = gl
        axs_repeats[r] = panel_setter(gl, r)
    end
    println("\r\033[2K\raxes built for trends.")
    fig[2, n_dim_vars[1]+1] |> Box
    colsize!(fig.layout, n_dim_vars[1] + 1, Fixed(2))
    gl = GridLayout()
    fig[1, n_dim_vars[1]+2] = Label(fig, text="Processed after stacked"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+2] = gl
    axs_stacked = panel_setter(gl, 1)
    gl = GridLayout()
    fig[1, n_dim_vars[1]+3] = Label(fig, text="Reps overlayed"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+3] = gl
    axs_all = panel_setter(gl, 1; extra=true)
    return fig, Dict("repeats" => axs_repeats, "stacked" => axs_stacked, "all" => axs_all)
end

function set_axes_2axes!(vals::NamedTuple, panel_setter::Function, runinfo)
    fig = Figure()
    n_dim_vars = vals |> vs -> map(length, vs) |> Tuple
    name_vars = vals |> propertynames |> ns -> string.(ns)
    length(n_dim_vars) == 2 || throw(ArgumentError("Select only 2 vars for to make this table."))
    axs = Array{Dict}(undef, n_dim_vars)
    for (r, v_r) in enumerate(vals[1]), (c, v_c) in enumerate(vals[2])
        print("\r\033[2K\rbuilding axes for $r-$c")
        gl_rc = GridLayout()
        fig[r, c] = gl_rc
        Label(gl_rc[0, 1], text="$(name_vars[1])=$(v_r) | $(name_vars[2])=$(v_c)"; tellwidth=false, tellheight=false, halign=:center, valign=:top)
        gl = GridLayout()
        gl_rc[1, 1] = gl
        axs[r, c] = panel_setter(gl; row_cmpl=n_dim_vars[1] - r, col=c)
        gl_rc |> l -> colgap!(l, 0)
        gl_rc |> l -> rowgap!(l, 0)
        gl_rc |> l -> rowsize!(l, 0, 20)
    end
    fig.layout |> l -> colgap!(l, 8)
    fig.layout |> l -> rowgap!(l, 4)
    println("\r\033[2K\raxes built for the trends")
    return fig, axs
end

function set_panel_single_axis(gl::GridLayout; row_cmpl=0, col=1)
    gl |> clean_gridlayout!
    ax = Axis(gl[1, 1]; width=400, height=200)
    hidedecorations!(ax; label=false, ticklabels=false, ticks=true, grid=true, minorticks=true, minorgrid=false)
    if col == 1
        ax.yticklabelsvisible = true
        ax.ylabelvisible = true
    end
    if row_cmpl == 0
        ax.xticklabelsvisible = true
        ax.xlabelvisible = true
    end
    return Dict("ax" => ax)
end

function set_axis_stack_all!(_, panel_setter::Function, runinfo)
    fig = Figure()
    fig[0, 1] = Label(fig, text="$(runinfo.date) $(@sprintf("run%02d", runinfo.runid)) IB=$(@sprintf("%.3f", runinfo.IB))A $(runinfo.tag_head)"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    gl = GridLayout()
    fig[1, 1] = Label(fig, text="Processed after stacked"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
    fig[2, 1] = gl
    axs_stacked = panel_setter(gl, 1)
    gl = GridLayout()
    fig[1, 2] = Label(fig, text="Reps overlayed"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
    fig[2, 2] = gl
    axs_all = panel_setter(gl, 2)
    return fig, Dict("stacked" => axs_stacked, "all" => axs_all)
end

function set_panel_trend_sidepeak_nvlp!(gl::GridLayout, col::Int; extra=false)
    gl |> clean_gridlayout!
    w, h = (400, 200)
    kwargs_wh = (width=w, height=h, yticklabelspace=40.0)
    col_freq = extra ? 3 : 2
    ax_evol_sum_dens = Axis(gl[1, 1]; kwargs_wh..., ylabel="density sum")
    ax_evol_weight = Axis(gl[2, 1]; kwargs_wh..., ylabel="side peak \nweight")
    ax_evol_height = Axis(gl[3, 1]; kwargs_wh..., ylabel="side peak \nheight")
    ax_evol_width = Axis(gl[4, 1]; kwargs_wh..., ylabel="side peak \nwidth (μm⁻¹)")
    ax_evol_wavenum = Axis(gl[5, 1]; kwargs_wh..., ylabel="side peak \nwavenum (μm⁻¹)")
    ax_evol_sizes = Axis(gl[6, extra ? 2 : 1]; kwargs_wh..., ylabel="envelope size (μm)")
    if extra
        ax_evol_extra_weight = Axis(gl[2, 2]; kwargs_wh..., ylabel="side peak \nweight")
        ax_evol_extra_height = Axis(gl[3, 2]; kwargs_wh..., ylabel="side peak \nheight")
        ax_evol_extra_width = Axis(gl[4, 2]; kwargs_wh..., ylabel="side peak \nwidth (μm⁻¹)")
        ax_evol_extra_wavenum = Axis(gl[5, 2]; kwargs_wh..., ylabel="side peak \nwavenum (μm⁻¹)")
    end
    ax_freq_weight = Axis(gl[2, col_freq]; kwargs_wh...)
    ax_freq_height = Axis(gl[3, col_freq]; kwargs_wh...)
    ax_freq_width = Axis(gl[4, col_freq]; kwargs_wh...)
    ax_freq_wavenum = Axis(gl[5, col_freq]; kwargs_wh...)
    ax_freq_sizes = Axis(gl[6, col_freq]; kwargs_wh...)
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

function set_panel_trend_nvlp!(gl::GridLayout; col::Int=1, row_cmpl::Int=0, extra=false)
    gl |> clean_gridlayout!
    w, h = (400, 200)
    kwargs_wh = (width=w, height=h, yticklabelspace=40.0)
    col_freq = extra ? 3 : 2
    ax_evol_sizes = Axis(gl[1, extra ? 2 : 1]; kwargs_wh..., ylabel="envelope size (μm)")
    ax_freq_sizes = Axis(gl[1, col_freq]; kwargs_wh...)
    rowgap!(gl, 4)
    colgap!(gl, 4)
    dict_axs = Dict(
        "evol-sizes" => ax_evol_sizes,
        "freq-sizes" => ax_freq_sizes,
    )
    for ax in values(dict_axs)
        hideydecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
        hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
        if col == 1
            ax.ylabelvisible = true
        end
        if row_cmpl == 0
            ax.xlabelvisible = true
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

function plot_shade_range!(axs, sample, clr)
    step = sample |> diff |> s -> filter(!iszero, s) |> minimum
    head = minimum(sample) - step / 2
    tail = maximum(sample) + step / 2
    for ax in axs
        vspan!(ax, head, tail; color=clr)
    end
end

function plot_trends_sidepeak!(axs::Dict, trend::Dict, istp; to_clean=false, alpha=1.0, is_stacked=false, to_legend=false, to_overlay=false)
    hue_theme = hue_theme_istp[istp]
    clr_mmt = Oklch(0.52, 0.14, hue_theme)
    clr_fit = (:springgreen3, 1.0)
    clr_theme = Oklch(0.52, 0.14, hue_theme)
    clr_shade_selected = RGBAf(Oklch(0.95, 0.1, hue_theme), 0.2)
    axs_sidepeaks_evol = axs |> a -> matching_axes(a, r"(evol(-extra)?)-(weight|width|height|wavenum)")
    axs_sidepeaks_freq = axs |> a -> matching_axes(a, r"freq-(weight|width|height|wavenum)")
    if to_clean
        axs |> a -> matching_axes(a, r"(freq|(evol(-extra)?))-(dens-sum|weight|width|height|wavenum)") |> clear_axes!
        plot_shade_range!(axs_sidepeaks_evol, trend["t_vec_sel_sp"], clr_shade_selected)
    end
    lines!(axs["evol-dens-sum"], trend["t_vec"], trend["evol-all-dens-sum"]; color=(clr_theme, alpha), label="sum")
    lines!(axs["evol-weight"], trend["t_vec"], trend["evol-all-fit-weight"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-height"], trend["t_vec"], trend["evol-all-fit-height"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-width"], trend["t_vec"], trend["evol-all-fit-width"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-wavenum"], trend["t_vec"], trend["evol-all-fit-wavenum"]; color=(clr_fit, alpha), label="fit")
    if to_overlay # overlay plots are plotted separately for fit and moment
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
    lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-fit-weight"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-fit-height"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-fit-width"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-fit-wavenum"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-moment-weight"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-moment-height"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-moment-width"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    ylims!(axs["evol-weight"], -0.02, 0.22)
    ylims!(axs["evol-height"], -0.1, 1.1)
    ylims!(axs["evol-width"], 0.02, 0.205)
    ylims!(axs["evol-wavenum"], 0.22, 0.38)
    if to_legend
        axislegend(axs["evol-dens-sum"]; position=:rt, framevisible=false, labelsize=14)
        axislegend(axs["freq-weight"]; position=:lt, framevisible=false, labelsize=14)
    end
    if to_overlay
        ylims!(axs["evol-extra-weight"], -0.02, 0.22)
        ylims!(axs["evol-extra-height"], -0.1, 1.1)
        ylims!(axs["evol-extra-width"], 0.02, 0.205)
        ylims!(axs["evol-extra-wavenum"], 0.22, 0.38)
    end
end

function plot_trends_nvlp!(axs::Dict, trend::Dict, istp; to_clean=false, alpha=1.0, is_stacked=false, to_legend=false, to_overlay=false)
    hue_theme = hue_theme_istp[istp]
    clr_theme = Oklch(0.52, 0.14, hue_theme)
    clr_shade_selected = RGBAf(Oklch(0.95, 0.1, hue_theme), 0.2)
    clr_theme1 = Oklch(0.52, 0.14, hue_theme - 20)
    clr_theme2 = Oklch(0.52, 0.14, hue_theme + 20)
    if to_clean
        axs |> a -> matching_axes(a, r"(freq|evol)-(sizes)") |> clear_axes!
        plot_shade_range!([axs["evol-sizes"]], trend["t_vec_sel_nvlp"], clr_shade_selected)
    end
    lines!(axs["evol-sizes"], trend["t_vec"], trend["evol-all-fit-size-x"]; color=(clr_theme1, alpha), label="fit")
    lines!(axs["evol-sizes"], trend["t_vec"], trend["evol-all-fit-size-y"]; color=(clr_theme2, alpha), label="fit")
    lines!(axs["freq-sizes"], trend["freq_query"], trend["freq-sel-fit-size-x"]; color=(clr_theme1, alpha), label="fit size x")
    lines!(axs["freq-sizes"], trend["freq_query"], trend["freq-sel-fit-size-y"]; color=(clr_theme2, alpha), label="fit size y")
    ylims!(axs["evol-sizes"], 1, 11)
    if to_legend
        axislegend(axs["freq-sizes"]; position=:lt, framevisible=false, labelsize=14)
    end
end

function plot_trend_all!(axs_trend::Dict, trend_reps::AbstractVector, trend_stacked_over_rep::Dict, istp)
    function set_tick_grid!(axs)
        axs_sidepeaks_evol = axs |> a -> matching_axes(a, r"(evol(-extra)?)-(weight|width|height|wavenum|dens-sum|sizes)")
        axs_sidepeaks_freq = axs |> a -> matching_axes(a, r"freq-(weight|width|height|wavenum|dens-sum|sizes)")
        for ax in axs_sidepeaks_evol
            ax.xticks = 0:50:200
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(5)
        end
        for ax in axs_sidepeaks_freq
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
            plot_trends_sidepeak!(axs, trend, istp; to_clean, alpha, to_legend, to_overlay)
            plot_trends_nvlp!(axs, trend, istp; to_clean, alpha, to_legend, to_overlay)
        end
    end
    # plot on the stacked axes
    axs = axs_trend["stacked"]
    plot_trends_sidepeak!(axs, trend_stacked_over_rep, istp; to_clean=true, alpha=1.0, to_legend=true)
    plot_trends_nvlp!(axs, trend_stacked_over_rep, istp; to_clean=true, alpha=1.0, to_legend=true)
    axs |> set_tick_grid!
end

function plot_trend_nvlp!(axs_trend::Dict, trend_reps::AbstractVector, trend_stacked_over_rep::Dict, istp)
    function set_tick_grid!(axs)
        axs_sidepeaks_evol = axs |> a -> matching_axes(a, r"(evol(-extra)?)-(weight|width|height|wavenum|dens-sum|sizes)")
        axs_sidepeaks_freq = axs |> a -> matching_axes(a, r"freq-(weight|width|height|wavenum|dens-sum|sizes)")
        for ax in axs_sidepeaks_evol
            ax.xticks = 0:50:200
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(5)
        end
        for ax in axs_sidepeaks_freq
            ax.xticks = 0:10:100
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(2)
        end
    end
    axs = axs_trend["all"]
    for r = axes(trend_reps, 1)
        trend = trend_reps[r]
        # plot on both the individual reps and all reps combined
        alpha = 0.5
        # shade the selected time points on all reps combined only once
        to_clean = r == 1
        to_legend = r == 1
        axs |> set_tick_grid!
        plot_trends_nvlp!(axs, trend, istp; to_clean, alpha, to_legend)
    end
    # plot on the stacked axes
    axs = axs_trend["stacked"]
    plot_trends_nvlp!(axs, trend_stacked_over_rep, istp; to_clean=true, alpha=1.0, to_legend=true)
    axs |> set_tick_grid!
end
