include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
## function definitions
function set_panel_solo_modl_masks!(gl::GridLayout)
    gl |> clean_gridlayout!
    kwargs_label = (;tellwidth=false, tellheight=false, halign=:left, valign=:bottom)
    label = Label(gl[0, 1:5]; tellwidth=false, tellheight=false, halign=:left, valign=:bottom)
    Label(gl[1, 1]; text="density", kwargs_label...)
    Label(gl[1, 2]; text="fringe", kwargs_label...)
    Label(gl[1, 3]; text="envelop", kwargs_label...)
    Label(gl[1, 4]; text="masked modulation", kwargs_label...)
    Label(gl[1, 5]; text="FT 2D", kwargs_label...)
    Label(gl[1, 6]; text="FT profile", kwargs_label...)
    ax_dens_core = Axis(gl[2, 1], width=120, height=240)
    ax_dens_core_masked_1 = Axis(gl[2, 2], width=120, height=240)
    ax_dens_core_masked_2 = Axis(gl[2, 3], width=120, height=240)
    ax_dens_core_masked_in = Axis(gl[2, 4], width=120, height=240)
    ax_modl = Axis(gl[2, 5], width=120, height=240)
    ax_prfl_ft_upright = Axis(gl[2, 6], width=240, height=240)
    colgap!(gl, 5)
    rowgap!(gl, 2)
    rowsize!(gl, 0, 16)
    rowsize!(gl, 1, 16)
    return Dict(
        "dens_core" => ax_dens_core,
        "dens_masked_1" => ax_dens_core_masked_1,
        "dens_masked_2" => ax_dens_core_masked_2,
        "dens_core_masked_in" => ax_dens_core_masked_in,
        "modl" => ax_modl,
        "upright" => ax_prfl_ft_upright,
        "label" => label)
end

function to_masked_clr(dens, mask, hue; sat_max=0.24, max=16, thres_alpha=0.1, l_max=1.0, l_min=0.0, alpha_base=0.1)
    size(dens) == size(mask) || throw(ArgumentError("dens and mask must have the same size"))
    norm = d -> clamp.(d, 0, max) / max
    dens_norm = dens |> norm
    alpha = (n, m) -> m ?
        (n > thres_alpha ? 1.0 : (n / thres_alpha * (1 - alpha_base) + alpha_base)) : 0
    shader = (n, m) -> Oklch(l_max - (l_max - l_min) * abs(n), sat_max * abs(n), hue) |> c -> RGBAf(c, alpha(n, m))
    return [shader(dens_norm[x, y], mask[x, y]) for x in 1:size(dens, 1), y in 1:size(dens, 2)]
end

function draw_solo_modl_mask!(axs::Dict{String}, extr::SoloExtract, info_solo, mask1, mask2, mask_prfl, prfls_masked, fit_tail; dens_max=16.0, peak_height_max=2)

    isnothing(extr.envelope) && return
    isnothing(extr.sidepeak) && return
    hue_1 = 105 # fringe
    hue_2 = 154 # envelope
    foreach(a -> a isa Axis && empty!(a), values(axs))
    essn = extr.essentials
    mask_hann = gen_win_hann_2d(essn.smwh_core)
    norm_modl = m -> m ./ (sum(essn.modl2d) * (essn.step_modl[2] / 2)^2)
    # norm_modl_offset = m -> clamp(m, 0, Inf) |> m -> m ./ (sum(m) * (essn.step_modl[2] / 2)^2)
    modl2d_norm = essn.modl2d |> norm_modl
    calc_dens_mask = mask -> (essn.dens2d_core .* mask_hann |> fft |> fftshift |> m -> m .* mask |> ifftshift |> ifft |> d -> real.(d))
    dens_masked_1 = mask1 |> calc_dens_mask |> m -> clamp.(m, 0, Inf)
    dens_masked_2 = mask2 |> calc_dens_mask |> m -> clamp.(m, 0, Inf)
    dens_masked_cleaned = .!(mask1 .| mask2) |> calc_dens_mask
    x_modl, y_modl = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_modl
    x_posi, y_posi = essn.smwh |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    x_posi_core, y_posi_core = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    y_modl_sm = (0:1:essn.smwh_core[2]) * essn.step_modl[2]
    hue_theme = hue_theme_istp[info_solo["istp"]]
    clrmap = gen_clrmap_solo(hue_theme; thres_alpha=0.05, alpha_base=0.05)
    clrmap_1 = gen_clrmap_solo(hue_1; thres_alpha=0.05, alpha_base=0.05)
    clrmap_2 = gen_clrmap_solo(hue_2; thres_alpha=0.05, alpha_base=0.05)
    clr_mark_nvlp = RGBAf(Oklch(0.52, 0.10, hue_theme + 90), 1.0)
    clr_moments = Oklch(0.52, 0.14, hue_theme)
    clr_prfl_masked = Oklch(0.52, 0.14, hue_theme)

    lims_core = (essn.smwh_core .+ 0.5) .* essn.step_posi |> l -> map(a -> [-a, a], l)
    lims_modl = ((0, 0.6), (-0.6, 0.6))

    nvlp = extr.envelope.params_asymm
    shade_mainpeak = zeros(length(y_modl_sm))
    shade_peaks = fit_prfl_modl_sidepeak_1d_model(y_modl_sm, extr.sidepeak.fit_tailess.params)
    band!(axs["upright"], y_modl_sm, 0, shade_mainpeak, color=(:gray, 0.1))
    band!(axs["upright"], y_modl_sm, shade_mainpeak, shade_peaks, color=(:darkseagreen1, 0.5))

    heatmap!(axs["dens_core"], x_posi_core, y_posi_core, (essn.dens2d_core .* gen_win_hann_2d(essn.smwh_core))'; colorrange=(0, dens_max), colormap=clrmap, rasterize=true)
    heatmap!(axs["dens_masked_1"], x_posi_core, y_posi_core, dens_masked_1'; colorrange=(0, dens_max), colormap=clrmap_1, rasterize=true)
    heatmap!(axs["dens_masked_2"], x_posi_core, y_posi_core, dens_masked_2'; colorrange=(0, dens_max), colormap=clrmap_2, rasterize=true)
    heatmap!(axs["dens_core_masked_in"], x_posi_core, y_posi_core, dens_masked_cleaned'; colorrange=(0, dens_max/2), colormap=clrmap, rasterize=true)

    # prfl_masked = essn.modl2d .* mask_prfl |>
    clr_max_modl = dens_max * 5 / 8
    clr_modl_mask1 = to_masked_clr(modl2d_norm, mask1, hue_1; max=clr_max_modl, thres_alpha=0.0)
    clr_modl_mask2 = to_masked_clr(modl2d_norm, mask2, hue_2; max=clr_max_modl, thres_alpha=0.0)
    clr_modl_mask_prfl = to_masked_clr(modl2d_norm, mask_prfl, hue_theme; max=clr_max_modl, thres_alpha=0.0)
    clr_modl_mask_rest = to_masked_clr(modl2d_norm, (@. !(mask_prfl | mask1 | mask2)), hue_theme; max=clr_max_modl, thres_alpha=0.0, sat_max=0)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_mask_rest; rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_mask1; rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_mask2; rasterize=true)
    heatmap!(axs["modl"], y_modl, x_modl, clr_modl_mask_prfl; rasterize=true)
    rng_prfl = essn.smwh_core[2] |> h -> h+1:(2*h + 1)
    lines!(axs["upright"], y_modl, essn.prfl_modl.side.normed_px, color=(:black, 0.4), linewidth=1)
    lines!(axs["upright"], y_modl, extr.sidepeak.prfl_norm_tailess_px, color=:black, linewidth=1)
    lines!(axs["upright"], y_modl, prfls_masked.prfl_main_normed_px, color=:black, linestyle=:dash, linewidth=1)
    lines!(axs["upright"], y_modl, prfls_masked.prfl_sidepeak_normed_px, color=clr_prfl_masked, linewidth=0.6)
    lines!(axs["upright"], y_modl, prfls_masked.prfl_sidepeak_normed_px - fit_tail(y_modl), color=clr_prfl_masked, linewidth=1.4)
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
    linkaxes!(axs["dens_core"], axs["dens_core_masked_in"], axs["dens_masked_1"], axs["dens_masked_2"], axs["dens_core_masked_in"])
    for ax in [axs["dens_core_masked_in"], axs["dens_masked_1"], axs["dens_masked_2"], axs["dens_core_masked_in"]]
        ax.yticklabelsvisible = false
        ax.xlabel = "μm"
    end
    axs["dens_core"].xlabel = "μm"
    for ax in [axs["modl"], axs["upright"]]
        ax.xlabel = "μm⁻¹"
    end
    sprint2f = (x) -> @sprintf("%.2f", x)
    axs["label"].text = @sprintf("%.03f A | %.01f ms | rep %s", info_solo["IB"], info_solo["t_hold"], info_solo["repeat"])
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

## investigate tail

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
    e -> fit_prfl_modl_twinpeak_decay_1d(y_modl, e.prfl_modl.side.normed_px, selector_tail_stack(y_modl); fit_stack_kwargs...)
    for c in axes(essn_2d_fmt, 1), t in axes(essn_2d_fmt, 3), i in axes(essn_2d_fmt, 4)
]
log_done("fit stacked modulation tails", t_stage)
tail_params_D = fit_prfl_modl_over_rep_1d |> fs -> map(f -> f.params[6], fs)
tail_params_λ = fit_prfl_modl_over_rep_1d |> fs -> map(f -> f.params[7], fs)

fig_tail = Figure()
ax_D = Axis(fig_tail[1,1]; title="D")
ax_λ = Axis(fig_tail[1,2]; title="λ")
c = 1
for c in axes(essn_2d_fmt, 1), i in axes(essn_2d_fmt, 4)
    # [ax_D ax_λ] |> clear_axes!
    clr_theme = Oklch(0.52, 0.14, hue_theme_istp[val_vars.istp[i]] + 5 * c)
    lines!(ax_D, val_vars.t_hold, tail_params_D[c,:,i]; color=RGBAf(clr_theme, 0.4))
    lines!(ax_λ, val_vars.t_hold, tail_params_λ[c,:,i]; color=RGBAf(clr_theme, 0.4))
end
fig_tail |> display

## Set a mask and show the results
using GLMakie
GLMakie.activate!()

crti_samples = [
    (4, 2, 19, 2, "strong fringe")
    (4, 2,  1, 2, "messy envelop")
    (4, 1, 30, 2, "small modulation peak")
    (3, 2,  4, 2, "")
    (3, 2,  3, 2, "small spacing")
    (5, 3, 34, 2, "large spacing side peak")
]
fig = Figure()
gl = GridLayout(fig[1, 1])
axs = set_panel_solo_modl_masks!(gl)
foreach(a -> a isa Axis && empty!(a), values(axs))

c, r, t, i = (3, 2,  3, 2)
extr = extr_fmt[c, r, t, i]
info = info_fmt[c, r, t, i]
essn = extr.essentials
fit_tail = y -> fit_prfl_modl_twinpeak_decay_1d_tail(y, fit_prfl_modl_over_rep_1d[c, t, i].params)
x_modl, y_modl = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_modl
x_posi, y_posi = essn.smwh |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
x_posi_core, y_posi_core = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi

mask_fringe = [((x - 0.37)/0.16)^2 + ((y - 0.32)/0.1)^2 < 1 for y in y_modl, x in x_modl] |> copy_symmetric_2d
mask_envelop = [
    begin
        θ = -20 * π / 180
        rot = [
            cos(θ)  -sin(θ)
            sin(θ)   cos(θ)
        ]
        x_r, y_r = rot * [x, y]
        ((x_r/0.4)^2 + (y_r/0.15)^2 < 1) | ((x/0.55)^2 + (y/0.1)^2 < 1)
    end
    for y in y_modl, x in x_modl] |> copy_symmetric_2d
mask_sidepeak = [@. abs(x) < 0.51 for y in y_modl, x in x_modl] |> m -> @. (m & !(mask_fringe | mask_envelop))
mask_main = [@. (abs(x) < 0.51) for y in y_modl, x in x_modl] |> m -> @. (m & !mask_fringe)
calc_prfl_masked = (modl, mask) -> sum(modl .* mask; dims=2) ./ sum(Int.(mask); dims=2) |> vec
function calc_prfl_norm_px_masked(modl2d, mask_sp, mask_main, step_modl)
    prfl_main = calc_prfl_masked(modl2d, mask_main)
    prfl_sidepeak = calc_prfl_masked(modl2d, mask_sp)
    norm_prfl_main = prfl_main |> sum |> s -> s * step_modl[2] / 2
    prfl_main_normed_px = prfl_main ./ norm_prfl_main
    prfl_sidepeak_normed_px = prfl_sidepeak ./ norm_prfl_main
    (;norm_prfl_main, prfl_main_normed_px, prfl_sidepeak_normed_px)
end
prfls_masked = calc_prfl_norm_px_masked(essn.modl2d, mask_sidepeak, mask_main, essn.step_modl)
draw_solo_modl_mask!(axs, extr, info, mask_fringe, mask_envelop, mask_sidepeak, prfls_masked, fit_tail; dens_max=16.0, peak_height_max=2)
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
