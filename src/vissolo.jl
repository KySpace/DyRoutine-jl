using Printf
using Colors: Oklch
using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using Colors: Oklch
using LaTeXStrings

function set_axis_full(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer}, panel_setter::Function)
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
    for t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
        print("\r\033[2Kbuilding stack axis for $t")
        gl = GridLayout()
        fig[1, 3*(i-1)+2][t, 1] = gl
        axs_stacked[t, i] = panel_setter(gl)
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
    ax_modl = Axis(gl[1, 3])
    ax_dens = Axis(gl[1, 2])
    ax_prfl_ft_upright = Axis(gl[1, 4])
    ax_prfl_ft_sideway = Axis(gl[1, 1])
    colsize!(gl, 1, Fixed(200))
    colsize!(gl, 2, Fixed(80))
    colsize!(gl, 3, Fixed(240))
    colsize!(gl, 4, Fixed(240))
    colgap!(gl, 5)
    rowsize!(gl, 1, Fixed(160))
    return Dict("dens" => ax_dens, "modl" => ax_modl, "upright" => ax_prfl_ft_upright, "sideway" => ax_prfl_ft_sideway)
end

function draw_solo_modl!(axs::Dict{String,Axis}, extr::SoloExtract, info_solo)
    isnothing(extr.envelope) && return
    isnothing(extr.sidepeak) && return

    foreach(empty!, values(axs))
    essn = extr.essentials
    modl2d_norm = essn.modl2d |> m -> m ./ (sum(m) * (essn.step_modl / 2)^2)
    x, y = essn.smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x, y) .* essn.step_posi
    x_modl, y_modl = (x, y) .* essn.step_modl
    y_modl_sm = (0:1:essn.smwh[2]) * essn.step_modl
    hue_theme = hue_theme_istp[info_solo["istp"]]
    clrmap = gen_clrmap_solo(hue_theme)
    clr_mark_nvlp = RGBAf(Oklch(0.52, 0.10, hue_theme + 90), 1.0)
    clr_moments = Oklch(0.52, 0.14, hue_theme)

    nvlp = extr.envelope.params_asymm
    shade_mainpeak = extr.sidepeak.fit_tailess.fitfn_main(y_modl_sm)
    shade_peaks = extr.sidepeak.fit_tailess.fitfn(y_modl_sm)
    band!(axs["upright"], y_modl_sm, 0, shade_mainpeak, color=(:gray, 0.1))
    band!(axs["upright"], y_modl_sm, shade_mainpeak, shade_peaks, color=(:darkseagreen1, 0.5))

    heatmap!(axs["dens"], x_posi, y_posi, essn.dens2d'; colorrange=(0, 16.0), colormap=clrmap, rasterize=true)
    draw_rotated_ellipse_corners!(axs["dens"], nvlp.cent, nvlp.size, nvlp.rotation; color=:white, linewidth=4)
    draw_rotated_ellipse_corners!(axs["dens"], nvlp.cent, nvlp.size, nvlp.rotation; color=clr_mark_nvlp, linewidth=2)

    heatmap!(axs["modl"], y_modl_sm, x_modl, modl2d_norm[essn.smwh[2]+1:end, :]; colorrange=(0, 10.0), colormap=clrmap, rasterize=true)
    lines!(axs["upright"], y_modl_sm, essn.prfl_modl_norm_px[essn.smwh[2]+1:end], color=(:black, 0.4), linewidth=1)
    lines!(axs["sideway"], essn.prfl_modl_norm_px[essn.smwh[2]+1:end], y_modl_sm, color=(:black, 0.4), linewidth=1)
    lines!(axs["upright"], y_modl_sm, extr.sidepeak.prfl_norm_tailess_px[essn.smwh[2]+1:end], color=:black, linewidth=1)
    lines!(axs["sideway"], extr.sidepeak.prfl_norm_tailess_px[essn.smwh[2]+1:end], y_modl_sm, color=:black, linewidth=1)
    axs["sideway"].yreversed = true
    axs["sideway"] |> hidedecorations!
    axs["modl"] |> hidedecorations!
    axs["dens"] |> hidedecorations!
    axs["upright"].yticklabelsvisible = false
    axs["upright"].xticklabelsvisible = false
    axs["dens"].aspect = DataAspect()
    xlims!(axs["upright"], 0, 0.6)
    xlims!(axs["modl"], 0, 0.6)
    xlims!(axs["dens"], -5, 5)
    ylims!(axs["dens"], -10, 10)
    ylims!(axs["upright"], -0.2, 1.8)
    ylims!(axs["modl"], (-10.5, 10.5) .* essn.step_modl)
    ylims!(axs["sideway"], 0.15, 0.45)
    xlims!(axs["sideway"], 0.0, 1.5)
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.2))
    vlines!(axs["upright"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    hlines!(axs["upright"], 0.0; color=(:darkseagreen1, 0.5))
    vlines!(axs["upright"], extr.sidepeak.params_tailess.wavenum; color=(:mediumspringgreen, 1.0))
    vlines!(axs["sideway"], extr.sidepeak.params_tailess.height; color=(:mediumspringgreen, 1.0))
    mmt = extr.sidepeak.moments
    sp = extr.sidepeak.params_tailess
    mmt_coor_min = mmt.coor |> minimum
    mmt_coor_max = mmt.coor |> maximum
    errorbars!(axs["upright"], [mmt.wavenum], [1.7], [mmt.width], [mmt.width]; direction=:x, color=clr_moments, whiskerwidth=8)
    lines!(axs["sideway"], [mmt.height, mmt.height], [mmt_coor_min, mmt_coor_max]; color=(clr_moments, 1.0))
    band!(axs["sideway"], [0, mmt.height], mmt_coor_min |> a -> [a, a], mmt_coor_max |> a -> [a, a]; color=(clr_moments, 0.1))

    sprint2f = (x) -> @sprintf("%.2f", x)
    text!(axs["modl"], 0.35, -0.16; text="$(info_solo["t_hold"]) ms | rep $(info_solo["repeat"])", color=:black, strokewidth=0.6, strokecolor=:white, fontsize=24, align=(:center, :top))
    text!(axs["dens"], -4.8, 9.8; text="[$(nvlp.size[1] |> sprint2f), $(nvlp.size[2] |> sprint2f)] μm \nrss/sum: $(nvlp.rel_residue |> sprint2f)", color=clr_mark_nvlp, strokewidth=0.5, strokecolor=:white, font=:bold, fontsize=11, align=(:left, :top))
    text!(axs["sideway"], 1.45, 0.44; text="fit: $(sp.height |> sprint2f), $(sp.weight |> sprint2f)", color=:springgreen3, fontsize=14, align=(:right, :top))
    text!(axs["sideway"], 1.45, 0.41; text="μ₀: $(mmt.height |> sprint2f), $(mmt.weight |> sprint2f)", color=clr_moments, fontsize=14, align=(:right, :top))
    text!(axs["sideway"], 1.45, 0.38; text="$(mmt_coor_min |> sprint2f)-$(mmt_coor_max |> sprint2f) μm⁻¹", color=clr_moments, fontsize=14, align=(:right, :top))
    text!(axs["upright"], 0.58, 1.4; text="fit: $(sp.wavenum |> sprint2f) ± $(sp.width |> sprint2f)", color=:springgreen3, fontsize=14, align=(:right, :top))
    text!(axs["upright"], 0.58, 1.2; text="rss/sum: $(sp.rel_residue |> sprint2f)", color=:springgreen3, fontsize=14, align=(:right, :top))
    text!(axs["upright"], 0.58, 1.6; text="μ₁, μ₂: $(mmt.wavenum |> sprint2f) ± $(mmt.width |> sprint2f)", color=clr_moments, fontsize=14, align=(:right, :top))
end

function draw_solo_essn_2d!(axs::Dict{String,Axis}, essn::SoloEssentials, info_solo)
    foreach(empty!, values(axs))
    modl2d_norm = essn.modl2d |> m -> m ./ (sum(m) * (essn.step_modl / 2)^2)
    x, y = essn.smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x, y) .* essn.step_posi
    x_modl, y_modl = (x, y) .* essn.step_modl
    y_modl_sm = (0:1:essn.smwh[2]) * essn.step_modl
    clrmap = gen_clrmap_solo(hue_theme_istp[info_solo["istp"]])
    heatmap!(axs["dens"], x_posi, y_posi, essn.dens2d'; colorrange=(0, 16.0), colormap=clrmap)
    heatmap!(axs["modl"], y_modl_sm, x_modl, modl2d_norm[essn.smwh[2]+1:end, :]; colorrange=(0, 10.0), colormap=:binary)
    axs["dens"].aspect = DataAspect()
    # axs["modl"].aspect = DataAspect()
    ylims!(axs["prfl_ft"], 0, 2.5)
    xlims!(axs["prfl_ft"], 0, 0.8)
    xlims!(axs["modl"], 0, 0.8)
    ylims!(axs["modl"], -0.5, 0.5)
    axs["prfl_ft"].yticksvisible = false
    axs["prfl_ft"].yticklabelsvisible = false
    axs["modl"] |> hidedecorations!
    # axs["dens"] |> hidedecorations!
    lines!(axs["prfl_ft"], y_modl_sm, essn.prfl_modl_norm_px |> fold_symmetric; color=:black)
    vlines!(axs["prfl_ft"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    hlines!(axs["modl"], [-10.5, 10.5] .* essn.step_modl; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    text!(axs["dens"], 0, 14; text=@sprintf("%i ms | rep %i", info_solo["t_hold"], info_solo["repeat"]), color=:black, fontsize=24, align=(:center, :bottom))
end
