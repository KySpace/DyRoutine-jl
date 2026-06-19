include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
## function definitions
function set_panel_solo_modl_masks!(gl::GridLayout)
    gl |> clean_gridlayout!
    label = Label(gl[0, 1:5]; tellwidth=false, tellheight=false, halign=:left, valign=:bottom)
    ax_dens_core = Axis(gl[1, 1], width=120, height=240)
    ax_dens_core_masked_1 = Axis(gl[1, 2], width=120, height=240)
    ax_dens_core_masked_2 = Axis(gl[1, 3], width=120, height=240)
    ax_dens_core_masked_in = Axis(gl[1, 4], width=120, height=240)
    ax_modl = Axis(gl[1, 5], width=120, height=240)
    ax_prfl_ft_upright = Axis(gl[1, 6], width=240, height=240)
    colgap!(gl, 5)
    rowgap!(gl, 2)
    rowsize!(gl, 0, 4)
    return Dict(
        "dens_core" => ax_dens_core,
        "dens_masked_1" => ax_dens_core_masked_1,
        "dens_masked_2" => ax_dens_core_masked_2,
        "dens_core_masked_in" => ax_dens_core_masked_in,
        "modl" => ax_modl,
        "upright" => ax_prfl_ft_upright,
        "label" => label)
end

function to_masked_clr(dens, mask, hue; sat_max=0.24, max=16, thres_alpha=0.1, l_max=0.8, l_min=0.0, alpha_base=0.1)
    size(dens) == size(mask) || throw(ArgumentError("dens and mask must have the same size"))
    norm = d -> clamp.(d, 0, max) / max
    dens_norm = dens |> norm
    alpha = (n, m) -> m ?
        (n > thres_alpha ? 1.0 : (n / thres_alpha * (1 - alpha_base) + alpha_base)) : 0
    shader = (n, m) -> Oklch(l_max - (l_max - l_min) * abs(n), sat_max * abs(n), hue) |> c -> RGBAf(c, alpha(n, m))
    return [shader(dens_norm[x, y], mask[x, y]) for x in 1:size(dens, 1), y in 1:size(dens, 2)]
end

function draw_solo_modl_mask!(axs::Dict{String}, extr::SoloExtract, info_solo, mask1, mask2; dens_max=16.0, peak_height_max=2)

    isnothing(extr.envelope) && return
    isnothing(extr.sidepeak) && return
    hue_1 = 105 # fringe
    hue_2 = 154 # envelope
    foreach(a -> a isa Axis && empty!(a), values(axs))
    essn = extr.essentials
    mask_hann = gen_win_hann_2d(essn.smwh_core)
    norm_modl = m -> map(p -> clamp(p, 0, Inf), m) |> m -> m ./ (sum(m) * (essn.step_modl[2] / 2)^2)
    modl2d_norm = essn.modl2d |> norm_modl
    calc_dens_mask = mask -> (essn.dens2d_core .* mask_hann |> fft |> fftshift |> m -> m .* mask |> ifftshift |> ifft |> d -> real.(d) |> norm_modl)
    dens_masked_1 = calc_dens_mask(mask1)
    dens_masked_2 = calc_dens_mask(mask2)
    x_modl, y_modl = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_modl
    x_posi, y_posi = essn.smwh |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    x_posi_core, y_posi_core = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    y_modl_sm = (0:1:essn.smwh_core[2]) * essn.step_modl[2]
    hue_theme = hue_theme_istp[info_solo["istp"]]
    clrmap = gen_clrmap_solo(hue_theme)
    clrmap_1 = gen_clrmap_solo(hue_1)
    clrmap_2 = gen_clrmap_solo(hue_2)
    clr_mark_nvlp = RGBAf(Oklch(0.52, 0.10, hue_theme + 90), 1.0)
    clr_moments = Oklch(0.52, 0.14, hue_theme)

    lims_core = (essn.smwh_core .+ 0.5) .* essn.step_posi |> l -> map(a -> [-a, a], l)
    lims_modl = ((0, 0.6), (-0.6, 0.6))

    nvlp = extr.envelope.params_asymm
    shade_mainpeak = extr.sidepeak.fit_tailess.fitfn_main(y_modl_sm)
    shade_peaks = extr.sidepeak.fit_tailess.fitfn(y_modl_sm)
    band!(axs["upright"], y_modl_sm, 0, shade_mainpeak, color=(:gray, 0.1))
    band!(axs["upright"], y_modl_sm, shade_mainpeak, shade_peaks, color=(:darkseagreen1, 0.5))

    heatmap!(axs["dens_core"], x_posi_core, y_posi_core, (essn.dens2d_core .* gen_win_hann_2d(essn.smwh_core))'; colorrange=(0, dens_max), colormap=clrmap, rasterize=true)
    heatmap!(axs["dens_masked_1"], x_posi_core, y_posi_core, dens_masked_1'; colorrange=(0, dens_max), colormap=clrmap_1, rasterize=true)
    # heatmap!(axs["dens_masked_2"], x_posi_core, y_posi_core, dens_masked_2'; colorrange=(0, dens_max), colormap=clrmap_2, rasterize=true)

    clr_max_modl = dens_max * 5 / 8
    clr_modl_mask1 = to_masked_clr(modl2d_norm, mask1, hue_1; max=clr_max_modl, thres_alpha=0.0)
    clr_modl_mask2 = to_masked_clr(modl2d_norm, mask2, hue_2; max=clr_max_modl, thres_alpha=0.0)
    heatmap!(axs["modl"], y_modl, x_modl, modl2d_norm; colorrange=(0, clr_max_modl), colormap=clrmap, rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_mask1; rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_mask2; rasterize=true)
    hlines!(axs["modl"], (essn.smw_modl+0.5)*essn.step_modl[1]; color=(:black, 0.4), linewidth=1)
    hlines!(axs["modl"], -(essn.smw_modl+0.5)*essn.step_modl[1]; color=(:black, 0.4), linewidth=1)
    lines!(axs["upright"], y_modl_sm, essn.prfl_modl_norm_px[essn.smwh_core[2]+1:end], color=(:black, 0.4), linewidth=1)
    lines!(axs["upright"], y_modl_sm, extr.sidepeak.prfl_norm_tailess_px[essn.smwh_core[2]+1:end], color=:black, linewidth=1)
    # axs["modl"] |> hidedecorations!
    axs["dens_core"] |> hidedecorations!
    axs["upright"].yticklabelsvisible = false
    axs["upright"].xticklabelsvisible = false
    limits!(axs["dens_core"], lims_core...)
    limits!(axs["dens_core_masked_in"], lims_core...)
    limits!(axs["dens_masked_1"], lims_core...)
    limits!(axs["dens_masked_2"], lims_core...)
    limits!(axs["modl"], lims_modl...)
    limits!(axs["upright"], (0, 0.6), (-0.2, peak_height_max + 0.2))
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.2))
    vlines!(axs["upright"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    hlines!(axs["upright"], 0.0; color=(:darkseagreen1, 0.5))
    vlines!(axs["upright"], extr.sidepeak.params_tailess.wavenum; color=(:mediumspringgreen, 1.0))
    # sp = extr.sidepeak.params_tailess
    # make_dual = a -> [a, a]
    # band!(axs["upright"], [mmt_coor_min, mmt_coor_max], 0 |> make_dual, -0.02 |> make_dual; color=(clr_moments, 1))
    linkaxes!(axs["dens_core"], axs["dens_core_masked_in"], axs["dens_masked_1"], axs["dens_masked_2"])
    sprint2f = (x) -> @sprintf("%.2f", x)
    axs["label"].text = "$(@sprintf("%.1f", info_solo["t_hold"])) ms | rep $(info_solo["repeat"])"
end

function copy_symmetric_2d(modl::AbstractMatrix{<:Real})
    size(modl, 1) |> isodd || error("Expect array of odd height")
    smh_modl = (size(modl, 1) + 1) / 2 |> Int
    modl_symm = modl .+ reverse(modl; dims=(1,2))
    modl_symm[smh_modl, :] = modl_symm[smh_modl, :] ./ 2
    return modl_symm
end
function copy_symmetric_2d(modl::AbstractMatrix{<:Bool})
    size(modl, 1) |> isodd || error("Expect array of odd height")
    smh_modl = (size(modl, 1) + 1) / 2 |> Int
    modl_symm = modl .| reverse(modl; dims=(1,2))
    return modl_symm
end
# Want to see how the fringes or features in density corresponds to peaks in modulation
#



## Loading
## Set a mask and show the results

c, r, t, i = (4, 1, 10, 2)
extr = extr_fmt[c, r, t, i]
info = info_fmt[c, r, t, i]
essn = extr.essentials
x_modl, y_modl = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_modl
x_posi, y_posi = essn.smwh |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
x_posi_core, y_posi_core = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi

fig = Figure()
gl = GridLayout(fig[1, 1])

axs = set_panel_solo_modl_masks!(gl)
foreach(a -> a isa Axis && empty!(a), values(axs))

mask_fringe = [((x - 0.37)/0.14)^2 + ((y - 0.32)/0.1)^2 < 1 for y in y_modl, x in x_modl] |> copy_symmetric_2d
mask_envelop = fill(false, length(y_modl), length(x_modl))
draw_solo_modl_mask!(axs, extr, info, mask_fringe, mask_envelop; dens_max=16.0, peak_height_max=2)

# norm_modl = m -> m ./ (sum(m) * (essn.step_modl[2] / 2)^2)
# modl_fringe = essn.modl2d .* mask_fringe |> norm_modl
# mask_hann = gen_win_hann_2d(smwh_core)
# dens2d_fringe = essn.dens2d_core .* mask_hann |> fft |> fftshift |> m -> m .* mask_fringe |> ifftshift |> ifft |> d -> abs.(d)
# axs["modl"] |> empty!
# clr_mask_fringe = to_masked_clr(essn.modl2d, mask_fringe, 105; max= 1600, thres_alpha=0.1)
# lims_modl = ((0, 0.6), (-0.6, 0.6))
# limits!(axs["modl"], lims_modl...)
# heatmap!(axs["dens_masked_1"], x_posi, y_posi, dens2d_fringe')
# heatmap!(axs["modl"], x_modl, y_modl, clr_mask_fringe')

fig |> resize_to_layout!
fig |> display
