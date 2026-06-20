using Printf
using Colors: Oklch
using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using Colors: Oklch
using LaTeXStrings

function set_axis_full(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer}, panel_setter::Function; to_plot_stacked=true)
    CairoMakie.activate!()
    CairoMakie.activate!()
    fig = Figure()
    length(n_dim_vars) == 3 || throw(ArgumentError("n_dim_vars must be a 3-tuple"))
    axs_solo = Array{Dict}(undef, n_dim_vars)
    axs_stacked = Array{Dict}(undef, n_dim_vars[2:end])
    for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
        print("\r\033[2Kbuilding solo axis for rep $r, $t")
        gl = GridLayout()
        # fig[1, 1][t, (r-1)*n_dim_vars[3]+i] = gl
        fig[1, 3*(i-1)+1][t, r] = gl
        axs_solo[r, t, i] = panel_setter(gl)
    end
    if to_plot_stacked
        for t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
            print("\r\033[2Kbuilding stack axis for $t")
            gl = GridLayout()
            fig[1, 3*(i-1)+2][t, 1] = gl
            axs_stacked[t, i] = panel_setter(gl)
        end
    end
    colsize!(fig.layout, 3, Fixed(2))
    return fig, axs_solo, axs_stacked
end

function set_panel_solo_essn_2d!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_dens = Axis(gl[1:2, 1])
    ax_modl = Axis(gl[1, 2])
    ax_prfl_ft = Axis(gl[2, 2])
    colsize!(gl, 1, Fixed(200))
    colsize!(gl, 2, Fixed(300))
    rowsize!(gl, 1, Fixed(200))
    rowsize!(gl, 2, Fixed(100))
    return Dict("dens" => ax_dens, "modl" => ax_modl, "prfl_ft" => ax_prfl_ft)
end

function set_panel_solo_modl!(gl::GridLayout)
    gl |> clean_gridlayout!
    label = Label(gl[0, 1:6]; tellwidth=false, tellheight=false, halign=:left, valign=:bottom)
    ax_prfl_ft_sideway = Axis(gl[1, 1], width=150, height=180)
    ax_dens = Axis(gl[1, 2], width=90, height=180)
    ax_dens_core = Axis(gl[1, 3], width=90, height=180)
    ax_dens_core_masked = Axis(gl[1, 4], width=90, height=180)
    ax_modl = Axis(gl[1, 5], width=90, height=180)
    ax_prfl_ft_upright = Axis(gl[1, 6], width=240, height=180)
    colgap!(gl, 5)
    rowgap!(gl, 2)
    rowsize!(gl, 0, 4)
    return Dict("dens" => ax_dens, "dens_core" => ax_dens_core, "dens_core_masked" => ax_dens_core_masked, "modl" => ax_modl, "upright" => ax_prfl_ft_upright, "sideway" => ax_prfl_ft_sideway, "label" => label)
end

function to_masked_clr(dens, mask, hue; sat_max=0.24, max=16, thres_alpha=0.1, l_max=1.0, l_min=0.0, alpha_base=0.1)
    size(dens) == size(mask) || throw(DimensionMismatch("dens size $(size(dens)) does not match mask size $(size(mask))."))
    dens_norm = clamp.(dens, 0, max) ./ max
    alpha = (n, m) -> m ? (thres_alpha <= 0 ? (n > 0 ? 1.0 : alpha_base) : (n > thres_alpha ? 1.0 : (n / thres_alpha * (1 - alpha_base) + alpha_base))) : 0.0
    shader = (n, m) -> Oklch(l_max - (l_max - l_min) * abs(n), sat_max * abs(n), hue) |> c -> RGBAf(c, alpha(n, m))
    return [shader(dens_norm[x, y], mask[x, y]) for x in 1:size(dens, 1), y in 1:size(dens, 2)]
end

function set_core_density_grid!(ax::Axis, lims_core)
    xlims, ylims = lims_core
    major_ticks(lims) = (ceil(Int, lims[1] / 10) * 10):10:(floor(Int, lims[2] / 10) * 10)
    minor_ticks(lims) = (ceil(Int, lims[1] / 2) * 2):2:(floor(Int, lims[2] / 2) * 2)
    ax.xticks = major_ticks(xlims)
    ax.yticks = major_ticks(ylims)
    ax.xminorticks = minor_ticks(xlims)
    ax.yminorticks = minor_ticks(ylims)
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    ax.xgridcolor = (:black, 0.22)
    ax.ygridcolor = (:black, 0.22)
    ax.xminorgridcolor = (:black, 0.10)
    ax.yminorgridcolor = (:black, 0.10)
    ax.xgridwidth = 0.6
    ax.ygridwidth = 0.6
    ax.xminorgridwidth = 0.35
    ax.yminorgridwidth = 0.35
    return ax
end

function draw_solo_modl!(
    axs::Dict{String},
    extr::SoloExtract,
    info_solo;
    dens_max=16.0,
    peak_height_max=2,
    fit_tailess_model::Function=fit_prfl_modl_sidepeak_1d_model,
)
    isnothing(extr.envelope) && return
    isnothing(extr.sidepeak) && return

    foreach(a -> a isa Axis && empty!(a), values(axs))
    essn = extr.essentials
    modl2d_norm = essn.modl2d |> m -> m ./ (sum(m) * (essn.step_modl[2] / 2)^2)
    masks = essn.prfl_modl.mask
    x_modl, y_modl = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_modl
    x_posi, y_posi = essn.smwh |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    x_posi_core, y_posi_core = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    y_modl_sm = (0:1:essn.smwh_core[2]) * essn.step_modl[2]
    hue_theme = hue_theme_istp[info_solo["istp"]]
    clrmap_dens = gen_clrmap_solo(hue_theme; thres_alpha=0.05, alpha_base=0.05)
    clr_mark_nvlp = RGBAf(Oklch(0.52, 0.10, hue_theme + 90), 1.0)
    clr_moments = Oklch(0.52, 0.14, hue_theme)
    hue_fringe = 105
    hue_center = 154

    lims_full = (essn.smwh .+ 0.5) .* essn.step_posi |> l -> map(a -> [-a, a], l)
    lims_core = (essn.smwh_core .+ 0.5) .* essn.step_posi |> l -> map(a -> [-a, a], l)

    nvlp = extr.envelope.params_asymm
    mask_hann = gen_win_hann_2d(essn.smwh_core)
    calc_dens_mask = mask -> (essn.dens2d_core .* mask_hann |> fft |> fftshift |> m -> m .* mask |> ifftshift |> ifft |> d -> real.(d))
    mask_dens_core = @. !(masks.fringe | masks.center)
    dens_core_masked = calc_dens_mask(mask_dens_core)

    heatmap!(axs["dens"], x_posi, y_posi, essn.dens2d'; colorrange=(0, dens_max), colormap=clrmap_dens, rasterize=true)
    heatmap!(axs["dens_core"], x_posi_core, y_posi_core, essn.dens2d_core'; colorrange=(0, dens_max), colormap=clrmap_dens, rasterize=true)
    heatmap!(axs["dens_core_masked"], x_posi_core, y_posi_core, dens_core_masked'; colorrange=(0, dens_max / 2), colormap=clrmap_dens, rasterize=true)
    draw_bound!(axs["dens"], essn.offset_cent_core, essn.smwh_core, essn.step_posi; color=:black, linewidth=0.5)
    draw_rotated_ellipse_corners!(axs["dens"], nvlp.cent, nvlp.size, nvlp.rotation; color=:white, linewidth=4)
    draw_rotated_ellipse_corners!(axs["dens"], nvlp.cent, nvlp.size, nvlp.rotation; color=clr_mark_nvlp, linewidth=2)

    clr_max_modl = dens_max * 5 / 8
    clr_modl_rest = to_masked_clr(modl2d_norm, (@. !(masks.side | masks.fringe | masks.center)), hue_theme; max=clr_max_modl, thres_alpha=0.0, sat_max=0)
    clr_modl_fringe = to_masked_clr(modl2d_norm, masks.fringe, hue_fringe; max=clr_max_modl, thres_alpha=0.0)
    clr_modl_center = to_masked_clr(modl2d_norm, masks.center, hue_center; max=clr_max_modl, thres_alpha=0.0)
    clr_modl_prfl = to_masked_clr(modl2d_norm, masks.side, hue_theme; max=clr_max_modl, thres_alpha=0.0)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_rest; rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_fringe; rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_center; rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_prfl; rasterize=true)
    band!(
        axs["upright"],
        y_modl_sm,
        0,
        fit_tailess_model(y_modl_sm, extr.sidepeak.fit_tailess.params);
        color=(:darkseagreen1, 0.28),
    ) |> b -> translate!(b, 0, 0, -1)
    lines!(axs["upright"], y_modl_sm, essn.prfl_modl.main.normed_px[essn.smwh_core[2]+1:end], color=(:black, 0.35), linestyle=:dash, linewidth=1)
    lines!(axs["sideway"], essn.prfl_modl.main.normed_px[essn.smwh_core[2]+1:end], y_modl_sm, color=(:black, 0.35), linestyle=:dash, linewidth=1)
    lines!(axs["upright"], y_modl_sm, essn.prfl_modl.side.normed_px[essn.smwh_core[2]+1:end], color=(:black, 0.65), linewidth=0.8)
    lines!(axs["sideway"], essn.prfl_modl.side.normed_px[essn.smwh_core[2]+1:end], y_modl_sm, color=(:black, 0.65), linewidth=0.8)
    lines!(axs["upright"], y_modl_sm, extr.sidepeak.prfl_norm_tailess_px[essn.smwh_core[2]+1:end], color=:black, linewidth=1.4)
    lines!(axs["sideway"], extr.sidepeak.prfl_norm_tailess_px[essn.smwh_core[2]+1:end], y_modl_sm, color=:black, linewidth=1.4)
    axs["sideway"].yreversed = true
    axs["sideway"] |> hidedecorations!
    axs["modl"] |> hidedecorations!
    axs["dens"] |> hidedecorations!
    hidedecorations!(axs["dens_core"]; grid=false, minorgrid=false)
    hidedecorations!(axs["dens_core_masked"]; grid=false, minorgrid=false)
    axs["upright"].yticklabelsvisible = false
    axs["upright"].xticklabelsvisible = false
    axs["dens"].aspect = DataAspect()
    axs["dens_core"].aspect = DataAspect()
    axs["dens_core_masked"].aspect = DataAspect()
    limits!(axs["dens"], lims_full...)
    limits!(axs["dens_core"], lims_core...)
    limits!(axs["dens_core_masked"], lims_core...)
    set_core_density_grid!(axs["dens_core"], lims_core)
    set_core_density_grid!(axs["dens_core_masked"], lims_core)
    xlims!(axs["upright"], 0, 0.6)
    limits!(axs["modl"], (0, 0.6), (-0.6, 0.6))
    ylims!(axs["upright"], -0.2, peak_height_max + 0.2)
    ylims!(axs["sideway"], 0.15, 0.45)
    xlims!(axs["sideway"], 0.0, peak_height_max)
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.2))
    vlines!(axs["upright"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    hlines!(axs["upright"], 0.0; color=(:darkseagreen1, 0.5))
    vlines!(axs["upright"], extr.sidepeak.params_tailess.wavenum; color=(:mediumspringgreen, 1.0))
    vlines!(axs["sideway"], extr.sidepeak.params_tailess.height; color=(:mediumspringgreen, 1.0))
    mmt = extr.sidepeak.moments
    sp = extr.sidepeak.params_tailess
    mmt_coor_min = mmt.coor |> minimum
    mmt_coor_max = mmt.coor |> maximum
    make_dual = a -> [a, a]
    errorbars!(axs["upright"], [mmt.wavenum], [1.7], [mmt.width], [mmt.width]; direction=:x, color=clr_moments, whiskerwidth=8)
    lines!(axs["sideway"], [mmt.height, mmt.height], [mmt_coor_min, mmt_coor_max]; color=(clr_moments, 1.0))
    band!(axs["sideway"], [0, mmt.height], mmt_coor_min |> make_dual, mmt_coor_max |> make_dual, color=(clr_moments, 0.1))
    band!(axs["upright"], [mmt_coor_min, mmt_coor_max], 0 |> make_dual, -0.02 |> make_dual; color=(clr_moments, 1))

    sprint2f = (x) -> @sprintf("%.2f", x)
    axs["label"].text = @sprintf("%.03f A | %.01f ms | rep %s", info_solo["IB"], info_solo["t_hold"], info_solo["repeat"])
    text!(axs["dens"], 0.05, 0.05; text="[$(nvlp.size[1] |> sprint2f), $(nvlp.size[2] |> sprint2f)] μm \nrss/sum: $(nvlp.rel_residue |> sprint2f)", space=:relative, color=clr_mark_nvlp, strokewidth=0.5, strokecolor=:white, font=:bold, fontsize=11, align=(:left, :bottom))
    text!(axs["sideway"], 0.98, 0.78; text="fit: $(sp.height |> sprint2f), $(sp.weight |> sprint2f)", space=:relative, color=:springgreen3, fontsize=14, align=(:right, :top))
    text!(axs["sideway"], 0.98, 0.88; text="μ₀: $(mmt.height |> sprint2f), $(mmt.weight |> sprint2f)", space=:relative, color=clr_moments, fontsize=14, align=(:right, :top))
    text!(axs["sideway"], 0.98, 0.98; text="$(mmt_coor_min |> sprint2f)-$(mmt_coor_max |> sprint2f) μm⁻¹", space=:relative, color=clr_moments, fontsize=14, align=(:right, :top))
    text!(axs["upright"], 0.98, 0.78; text="rss/sum: $(sp.rel_residue |> sprint2f)", space=:relative, color=:springgreen3, fontsize=14, align=(:right, :top))
    text!(axs["upright"], 0.98, 0.88; text="fit: $(sp.wavenum |> sprint2f) ± $(sp.width |> sprint2f)", space=:relative, color=:springgreen3, fontsize=14, align=(:right, :top))
    text!(axs["upright"], 0.98, 0.98; text="μ₁, μ₂: $(mmt.wavenum |> sprint2f) ± $(mmt.width |> sprint2f)", space=:relative, color=clr_moments, fontsize=14, align=(:right, :top))
end

function draw_solo_essn_2d!(axs::Dict{String,Axis}, essn::SoloEssentials, info_solo; dens_max=16.0, peak_height_max=2)
    foreach(empty!, values(axs))
    modl2d_norm = essn.modl2d |> m -> m ./ (sum(m) * (essn.step_modl[2] / 2)^2)
    x_modl, y_modl = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_modl
    x_posi, y_posi = essn.smwh |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    y_modl_sm = (0:1:essn.smwh_core[2]) * essn.step_modl[2]
    hue_theme = hue_theme_istp[info_solo["istp"]]
    clrmap_dens = gen_clrmap_solo(hue_theme; thres_alpha=0.05, alpha_base=0.05)
    clr_moments = Oklch(0.52, 0.14, hue_theme)

    heatmap!(axs["dens"], x_posi, y_posi, essn.dens2d'; colorrange=(0, dens_max), colormap=clrmap_dens, rasterize=true)

    heatmap!(axs["modl"], y_modl_sm, x_modl, modl2d_norm[essn.smwh_core[2]+1:end, :]; colorrange=(0, dens_max * 5 / 8), colormap=clrmap_dens, rasterize=true)
    lines!(axs["upright"], y_modl_sm, essn.prfl_modl.side.normed_px[essn.smwh_core[2]+1:end], color=(:black, 1.0), linewidth=1)
    lines!(axs["sideway"], essn.prfl_modl.side.normed_px[essn.smwh_core[2]+1:end], y_modl_sm, color=(:black, 0.4), linewidth=1)
    axs["sideway"].yreversed = true
    axs["dens"].aspect = DataAspect()
    xlims!(axs["upright"], 0, 0.6)
    xlims!(axs["modl"], 0, 0.6)
    xlims!(axs["dens"], -5, 5)
    ylims!(axs["dens"], -10, 10)
    ylims!(axs["upright"], -0.2, peak_height_max + 0.2)
    ylims!(axs["modl"], (-10.5, 10.5) .* essn.step_modl[1])
    ylims!(axs["sideway"], 0.15, 0.45)
    xlims!(axs["sideway"], 0.0, peak_height_max)
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.2))
    # vlines!(axs["upright"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    # hlines!(axs["upright"], 0.0; color=(:darkseagreen1, 0.5))

    sprint2f = (x) -> @sprintf("%.2f", x)
    # text!(axs["modl"], 0.35, -0.16; text="$(@sprintf("%.1f", info_solo["t_hold"])) ms | rep $(info_solo["repeat"])", color=:black, strokewidth=0.6, strokecolor=:white, fontsize=24, align=(:center, :top))
    # text!(axs["sideway"], 1.45, 0.38; text="$(mmt_coor_min |> sprint2f)-$(mmt_coor_max |> sprint2f) μm⁻¹", color=clr_moments, fontsize=14, align=(:right, :top))
end
