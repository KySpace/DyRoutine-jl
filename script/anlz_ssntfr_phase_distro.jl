using GLMakie
using FFTW: fft, fftshift, ifftshift
using HDF5
using ImageFiltering
using JLD2
using LinearAlgebra: norm
using LsqFit: curve_fit, residuals, stderror
using Printf
using Statistics
import CairoMakie

GLMakie.activate!()

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
path_data = joinpath(path_root, "0204_interference", "result", "data.h5")

# commit 57d2e69f9017be9957b38674e01e6fbc3aae013d
# Slightly larger width but with more tukey parameter
path_output = joinpath(path_root, "AnlzRoutine", "42.MeanAbsl2D.Narrow")
path_fit_jld2 = joinpath(path_output, "SSNTFR_phase_distro_fit.jld2")

tag = "SSNTFR"
filename_plot_phase_distro = "$(tag)_phase_distro_table"
filename_plot_cohr_modl = "$(tag)_cohr_modl"
val_istp = ["162", "164"]
val_IB_ref = [
    5.310,
    5.312,
    5.314,
    5.316,
    5.317,
    5.318,
    5.319,
    5.320,
    5.322,
    5.324,
    5.326,
    5.328,
    5.330,
    5.332,
    5.334,
    5.338,
    5.342,
]
smwh = (150, 150)
smwh_dens_ft = (50, 50)
alpha_tukey = (0.2, 1.0)
mag = 22.06
pixsz = 6.5
bin = 1
sigma_center_filter = 5
use_common_xy_center = :free # :free, :fixed_x, or :fixed_xy
x_max_fit_peak = 10 # μm
x_max_fit_modl = 4 # μm
x_fit_offset = 0.0 # μm
smh_dens_strip = 30
y_strip_offset = -0.0 # μm
amp_gauss_init = 6.0
x0_gauss_init = 0.0
sigma_gauss_init = 12.0
bg_gauss_init = 0.5
# [amp_gauss, x0_gauss, sigma_gauss, bg_gauss]
fit_lower_gauss = [0.0, -10.0, 6.0, -2.0]
fit_upper_gauss = [25.0, 10.0, 20.0, 1.0]
amp_modl_init = 1.0
slope_modl_init = 0.0
quad_modl_init = 0.0
lambda_modl_init = 5.0
phi_modl_init = (0.0, pi)
# [amp_modl, slope_modl, quad_modl, lambda_modl, phi_modl]
fit_lower_modl = [0.0, -0.20, -0.050, 4.0, -2pi]
fit_upper_modl = [8.0, 0.20, 0.050, 7.5, 2pi]
lambda_hue_min = 4.0
lambda_hue_max = 7.5
lambda_hue_span = 260.0
polar_lightness = 0.74
polar_chroma = 0.12
polar_alpha = 1.0
polar_lightness_rss_bad = 0.92
polar_chroma_rss_bad = 0.04
polar_alpha_rss_bad = 0.1
rss_rel_ramp = (0.4, 0.5)
markersize_fit = 7
markersize_fit_selected = 13

# live inspector selections
ib, istp, idx_rep = (5, 1, 1)
y_row = 0.0
ylims_profile = (-1.0, 15.0)
hue_scheme = :lambda
phase_mode = :phi0 # :phi0 from modulation fit, or :phip referenced to fitted peak x0

function gauss_1d_model(x, p)
    (A, x0, σ, bg) = p
    @. A * exp(-((x - x0)/σ)^2) + bg
end
name_gauss_params(p) = (; A=p[1], x0=p[2], σ=p[3], bg=p[4])

function modl_vary_1d_model(x, p)
    (M, a, b, λ, φ) = p
    @. M * (1 + a * x + b * x^2) * cos(2π * x/λ - φ)
end
name_modl_params(p) = (; M=p[1], a=p[2], b=p[3], λ=p[4], φ=p[5])

function load_density_payload(path_data::AbstractString, val_istp::AbstractVector{<:AbstractString})
    name_dataset_by_istp = Dict(
        "162" => "im64us",
        "164" => "im62us",
    )

    h5open(path_data, "r") do file
        dens_loaded = map(val_istp) do istp
            read(file[name_dataset_by_istp[istp]])
        end
        _, _, n_rep, n_IB = size(first(dens_loaded))
        dens_raw = Array{Matrix{Float64}}(undef, n_IB, length(val_istp), n_rep)
        for idx_IB in 1:n_IB, idx_istp in eachindex(val_istp), idx_rep in 1:n_rep
            dens_raw[idx_IB, idx_istp, idx_rep] = Float64.(copy(@view dens_loaded[idx_istp][:, :, idx_rep, idx_IB]))
        end
        return dens_raw
    end
end

function gaussian_offset_1d(x, p)
    return @. p[1] * exp(-((x - p[2])^2) / (2 * p[3]^2)) + p[4]
end

##
function shader_cmpx(ampl, phase; l_max=0.7438, c_max=0.1255, hue_offset=0, prescale=(t -> t), alpha_base=0.1, thres_alpha=0.05)
    l = @pipe ampl |> prescale |> clamp(_, 0, 1) |> 1 - _
    alpha = 1 - l |> u -> u > thres_alpha ? 1.0 : (u / thres_alpha * (1 - alpha_base) + alpha_base)
    c = l < l_max ? (l / l_max * c_max) : ((1 - l) / (1 - l_max) * c_max)
    hue = @pipe mod(phase, 2pi)/2pi*360 + hue_offset |> mod(_, 360)
    RGBAf(Oklch(l, c, hue), alpha)
end
function draw_profile_inspector!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    y_dens::AbstractVector{<:Real},
    dens_core::AbstractMatrix,
    dens_core_ft_masked::AbstractMatrix,
    dens_core_ft_amp::AbstractMatrix,
    dens_core_ft_phs::AbstractMatrix,
    dens_core_ft_cmpx_mean::AbstractMatrix,
    dens_core_ft_absl_mean::AbstractMatrix,
    x_dens_ft::AbstractVector{<:Real},
    y_dens_ft::AbstractVector{<:Real},
    kx_ft::AbstractVector{<:Real},
    ky_ft::AbstractVector{<:Real},
    fit_peak::AbstractMatrix,
    xy_center::AbstractArray{<:Tuple{Int,Int},3},
    val_istp::AbstractVector;
    ib::Integer,
    istp::Integer,
    idx_rep::Integer,
    y_row::Real,
    smidx_mean_profile::Integer,
    y_strip_offset::Real,
    ylims_profile::Tuple{<:Real,<:Real},
    x_max_fit_peak::Real,
    x_max_fit_modl::Real,
    x_fit_offset::Real,
    hue_scheme::Symbol,
    lambda_hue_min::Real,
    lambda_hue_max::Real,
    lambda_hue_span::Real,
    polar_lightness::Real,
    polar_chroma::Real,
    polar_alpha::Real,
    polar_lightness_rss_bad::Real,
    polar_chroma_rss_bad::Real,
    polar_alpha_rss_bad::Real,
    rss_rel_ramp::Tuple{<:Real,<:Real},
    clrrng_ft2d_absl_mean::Tuple{<:Real,<:Real},
    clrrng_ft2d_cmpx_mean::Tuple{<:Real,<:Real},
    markersize_fit::Real,
    markersize_fit_selected::Real,
    x_center_px0::Real,
    cohr::AbstractArray{<:NamedTuple},
)
    ib in axes(dens_core, 1) || throw(ArgumentError("ib must be in $(axes(dens_core, 1)), got $ib."))
    istp in axes(dens_core, 2) || throw(ArgumentError("istp must be in $(axes(dens_core, 2)), got $istp."))
    size(dens_core, 2) == length(val_istp) || throw(DimensionMismatch(
        "dens_core second dimension $(size(dens_core, 2)) must match length(val_istp) $(length(val_istp)).",
    ))
    size(fit_peak) == size(dens_core) || throw(DimensionMismatch(
        "fit_peak size $(size(fit_peak)) must match dens_core size $(size(dens_core)).",
    ))
    size(dens_core_ft_masked) == size(dens_core) || throw(DimensionMismatch(
        "dens_core_ft_masked size $(size(dens_core_ft_masked)) must match dens_core size $(size(dens_core)).",
    ))
    size(dens_core_ft_amp) == size(dens_core) || throw(DimensionMismatch(
        "dens_core_ft_amp size $(size(dens_core_ft_amp)) must match dens_core size $(size(dens_core)).",
    ))
    size(dens_core_ft_phs) == size(dens_core) || throw(DimensionMismatch(
        "dens_core_ft_phs size $(size(dens_core_ft_phs)) must match dens_core size $(size(dens_core)).",
    ))
    size(dens_core_ft_cmpx_mean) == size(dens_core) || throw(DimensionMismatch(
        "dens_core_ft_cmpx_mean size $(size(dens_core_ft_cmpx_mean)) must match dens_core size $(size(dens_core)).",
    ))
    size(dens_core_ft_absl_mean) == size(dens_core) || throw(DimensionMismatch(
        "dens_core_ft_absl_mean size $(size(dens_core_ft_absl_mean)) must match dens_core size $(size(dens_core)).",
    ))
    size(cohr) == size(dens_core) || throw(DimensionMismatch(
        "cohr size $(size(cohr)) must match dens_core size $(size(dens_core)).",
    ))
    size(xy_center) == (size(dens_core, 1), size(dens_core, 2), length(first(dens_core))) || throw(DimensionMismatch(
        "xy_center size $(size(xy_center)) must match (IB, istp, rep) $((size(dens_core, 1), size(dens_core, 2), length(first(dens_core)))).",
    ))
    for idx in CartesianIndices(dens_core)
        isempty(dens_core[idx]) && continue
        length(fit_peak[idx]) == length(dens_core[idx]) || throw(DimensionMismatch(
            "fit_peak[$(Tuple(idx)...)] length $(length(fit_peak[idx])) must match dens_core length $(length(dens_core[idx])).",
        ))
        length(dens_core_ft_masked[idx]) == length(dens_core[idx]) || throw(DimensionMismatch(
            "dens_core_ft_masked[$(Tuple(idx)...)] length $(length(dens_core_ft_masked[idx])) must match dens_core length $(length(dens_core[idx])).",
        ))
        length(dens_core_ft_amp[idx]) == length(dens_core[idx]) || throw(DimensionMismatch(
            "dens_core_ft_amp[$(Tuple(idx)...)] length $(length(dens_core_ft_amp[idx])) must match dens_core length $(length(dens_core[idx])).",
        ))
        length(dens_core_ft_phs[idx]) == length(dens_core[idx]) || throw(DimensionMismatch(
            "dens_core_ft_phs[$(Tuple(idx)...)] length $(length(dens_core_ft_phs[idx])) must match dens_core length $(length(dens_core[idx])).",
        ))
        size(first(dens_core[idx])) == (length(y_dens), length(x_dens)) || throw(DimensionMismatch(
            "dens_core[$(Tuple(idx)...)] crop size $(size(first(dens_core[idx]))) must match " *
            "(length(y_dens), length(x_dens)) $((length(y_dens), length(x_dens))).",
        ))
        size(first(dens_core_ft_masked[idx])) == (length(y_dens_ft), length(x_dens_ft)) || throw(DimensionMismatch(
            "dens_core_ft_masked[$(Tuple(idx)...)] masked density size $(size(first(dens_core_ft_masked[idx]))) must match " *
            "(length(y_dens_ft), length(x_dens_ft)) $((length(y_dens_ft), length(x_dens_ft))).",
        ))
        size(first(dens_core_ft_amp[idx])) == (length(ky_ft), length(kx_ft)) || throw(DimensionMismatch(
            "dens_core_ft_amp[$(Tuple(idx)...)] crop FT size $(size(first(dens_core_ft_amp[idx]))) must match " *
            "(length(ky_ft), length(kx_ft)) $((length(ky_ft), length(kx_ft))).",
        ))
        size(first(dens_core_ft_phs[idx])) == (length(ky_ft), length(kx_ft)) || throw(DimensionMismatch(
            "dens_core_ft_phs[$(Tuple(idx)...)] crop FT size $(size(first(dens_core_ft_phs[idx]))) must match " *
            "(length(ky_ft), length(kx_ft)) $((length(ky_ft), length(kx_ft))).",
        ))
        size(dens_core_ft_cmpx_mean[idx]) == (length(ky_ft), length(kx_ft)) || throw(DimensionMismatch(
            "dens_core_ft_cmpx_mean[$(Tuple(idx)...)] size $(size(dens_core_ft_cmpx_mean[idx])) must match " *
            "(length(ky_ft), length(kx_ft)) $((length(ky_ft), length(kx_ft))).",
        ))
        size(dens_core_ft_absl_mean[idx]) == (length(ky_ft), length(kx_ft)) || throw(DimensionMismatch(
            "dens_core_ft_absl_mean[$(Tuple(idx)...)] size $(size(dens_core_ft_absl_mean[idx])) must match " *
            "(length(ky_ft), length(kx_ft)) $((length(ky_ft), length(kx_ft))).",
        ))
    end

    dens_vec = dens_core[ib, istp]
    dens_ft_masked_vec = dens_core_ft_masked[ib, istp]
    dens_ft_amp_vec = dens_core_ft_amp[ib, istp]
    dens_ft_phs_vec = dens_core_ft_phs[ib, istp]
    isempty(dens_vec) && throw(ArgumentError("dens_core[$ib, $istp] has no selected crops."))
    n_rep_profile = length(dens_vec)
    idx_rep = mod1(idx_rep, length(dens_vec))
    idx_row = argmin(abs.(y_dens .- y_row))
    idx_strip_center = argmin(abs.(y_dens .- y_strip_offset))
    idxs_center = max(1, idx_strip_center - smidx_mean_profile):min(length(y_dens), idx_strip_center + smidx_mean_profile)
    mask_fit_peak_plot = abs.(x_dens .- x_fit_offset) .<= x_max_fit_peak
    dens2d = dens_vec[idx_rep]
    fit_info = fit_peak[ib, istp][idx_rep]
    obs_rss_rel_ramp = Observable((Float64(rss_rel_ramp[1]), Float64(rss_rel_ramp[2])))

    gen_theme_clr(idx_istp::Integer, alpha::Real) =
        RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
    gen_theme_clrmap(idx_istp::Integer) =
        gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
    gen_hue_fit(fit_info, hue_scheme_live::Symbol) =
        if hue_scheme_live == :lambda
            lambda_norm = clamp((fit_info.fit_modl.params.λ - lambda_hue_min) / (lambda_hue_max - lambda_hue_min), 0, 1)
            lambda_hue_span * (1 - lambda_norm)
        elseif hue_scheme_live == :rep
            n_rep_profile > 1 ? 360 * (fit_info.idx_rep - 1) / (n_rep_profile - 1) : 0.0
        else
            throw(ArgumentError("Unknown hue_scheme $hue_scheme_live. Expected :lambda or :rep."))
        end
    gen_fit_color(fit_info, hue_scheme_live::Symbol) = begin
        ramp_start, ramp_end = obs_rss_rel_ramp[]
        rss_norm = clamp(
            (fit_info.fit_modl.rss_rel - ramp_start) / max(ramp_end - ramp_start, eps(Float64)),
            0,
            1,
        )
        lightness = polar_lightness + rss_norm * (polar_lightness_rss_bad - polar_lightness)
        chroma = polar_chroma + rss_norm * (polar_chroma_rss_bad - polar_chroma)
        alpha = polar_alpha + rss_norm * (polar_alpha_rss_bad - polar_alpha)
        RGBAf(Oklch(lightness, chroma, gen_hue_fit(fit_info, hue_scheme_live)), alpha)
    end
    gen_fit_phase(fit_info, phase_mode_live::Symbol) =
        if phase_mode_live == :phi0
            fit_info.fit_modl.params.φ
        elseif phase_mode_live == :phip
            fit_info.fit_modl.params.φ - 2pi * fit_info.fit_gauss.params.x0 / fit_info.fit_modl.params.λ
        else
            throw(ArgumentError("Unknown phase_mode $phase_mode_live. Expected :phi0 or :phip."))
        end
    gen_fit_polar_payload(
        idx_IB::Integer,
        idx_istp::Integer,
        idx_rep_selected::Integer,
        hue_scheme_live::Symbol,
        phase_mode_live::Symbol,
    ) = begin
        fits = fit_peak[idx_IB, idx_istp]
        ids_success = findall(f -> f.success, fits)
        theta = [mod(gen_fit_phase(fits[idx], phase_mode_live), 2pi) for idx in ids_success]
        radius = [abs(fits[idx].fit_modl.params.M) for idx in ids_success]
        radius_outer = fill(polar_radius_outer, length(theta))
        color = map(ids_success) do idx
            gen_fit_color(fits[idx], hue_scheme_live)
        end
        markersize = [idx == idx_rep_selected ? markersize_fit_selected : markersize_fit for idx in ids_success]
        return (; theta, radius, radius_outer, color, markersize)
    end
    gen_moment_payload(idx_IB::Integer, idx_istp::Integer) = begin
        moment = cohr[idx_IB, idx_istp]
        theta = [mod(moment.moment_angel, 2pi)]
        radius = [clamp(moment.moment_length, 0, 1) * moment_radius_scale]
        return (; theta, radius)
    end

    clr_mean = RGBAf(0.35, 0.35, 0.35, 0.62)
    clr_strip = RGBAf(0.86, 0.86, 0.86, 0.50)
    clr_fit_peak_span = RGBAf(0.86, 0.86, 0.86, 0.18)
    clr_fit_modl_span = RGBAf(0.86, 0.86, 0.86, 0.42)
    clr_fit = RGBAf(Oklch(0.60, 0.17, 145), 0.95)
    clr_moment = RGBAf(0.02, 0.02, 0.02, 0.92)
    polar_radius_max = maximum([
        abs(fit.fit_modl.params.M)
        for fits in fit_peak
        for fit in fits
        if fit.success
    ])
    polar_radius_max = max(polar_radius_max, eps(Float64))
    moment_radius_scale = polar_radius_max
    polar_radius_outer = 1.08 * polar_radius_max
    polar_radius_limit = 1.20 * polar_radius_max
    polar_rticks = (0:2:6, string.(0:2:6))
    step_x = median(diff(x_dens))
    step_y = median(diff(y_dens))
    y_strip_min = y_dens[first(idxs_center)] - step_y / 2
    y_strip_max = y_dens[last(idxs_center)] + step_y / 2
    x_fit_peak_min, x_fit_peak_max = (x_fit_offset - x_max_fit_peak, x_fit_offset + x_max_fit_peak)
    x_fit_modl_min, x_fit_modl_max = (x_fit_offset - x_max_fit_modl, x_fit_offset + x_max_fit_modl)
    mask_sidepeak = [sqrt(x^2 + y^2) > 0.1 for y in ky_ft, x in kx_ft]
    gen_x0_peak(idx_IB::Integer, idx_istp::Integer) =
        [fit_peak[idx_IB, idx_istp][idx].fit_gauss.params.x0 for idx in 1:n_rep_profile]
    gen_ft_rgba(amp, phs; rng_amp=clrrng_ft2d_cmpx_mean) = begin
        amp_max = isnothing(rng_amp) ? 
            max(amp[mask_sidepeak] |> vec |> maximum, eps(Float64)) :
            rng_amp[2] / 2
        rgba = shader_cmpx.(amp ./ amp_max, phs; prescale=ampl_prescaler)
        phs_sample = range(-π, π, 256)
        amp_sample = range(0, 2amp_max, 256)
        clrmp = [shader_cmpx.(a ./ amp_max, φ; prescale=ampl_prescaler) for a in amp_sample, φ in phs_sample]
        (rgba, (;amp=amp_sample, phs=phs_sample, colormap=clrmp))
    end
    gen_ft_rgba(ft::AbstractMatrix{<:Complex}; kwargs...) = gen_ft_rgba(abs.(ft), angle.(ft); kwargs...)
    gen_ft_rgba_img(args...; kwargs...) = first(gen_ft_rgba(args...; kwargs...))
    gen_ft_rgba_clr(args...; kwargs...) = last(gen_ft_rgba(args...; kwargs...))
    gen_fit_gauss_text(fit_info_live) = @sprintf(
        "A=%.2f\nx0=%.2f\nσ=%.1f\nbg=%.3f",
        fit_info_live.fit_gauss.params.A,
        fit_info_live.fit_gauss.params.x0,
        fit_info_live.fit_gauss.params.σ,
        fit_info_live.fit_gauss.params.bg,
    )
    gen_fit_modl_text_left(fit_info_live) = @sprintf(
        "M=%g\na=%g\nb=%g",
        fit_info_live.fit_modl.params.M,
        fit_info_live.fit_modl.params.a,
        fit_info_live.fit_modl.params.b,
    )
    gen_fit_modl_text_right(fit_info_live) = @sprintf(
        "λ=%.2f\nφ0=%.2fπ | %.1f°\nφp=%.2fπ | %.1f°",
        fit_info_live.fit_modl.params.λ,
        fit_info_live.fit_modl.params.φ / 2,
        fit_info_live.fit_modl.params.φ / 2pi * 360,
        gen_fit_phase(fit_info_live, :phip) / 2,
        gen_fit_phase(fit_info_live, :phip) / 2pi * 360,
    )
    gen_fit_modl_text_rss(fit_info_live) = @sprintf("rss=%.2f", fit_info_live.fit_modl.rss_rel)

    obs_idx_IB = Observable(ib)
    obs_idx_istp = Observable(istp)
    obs_idx_rep = Observable(idx_rep)
    obs_hue_scheme = Observable(hue_scheme)
    obs_phase_mode = Observable(phase_mode)
    obs_ft_cmpx_view = Observable(:complex)
    obs_idx_row = Observable(idx_row)
    obs_val_row = Observable(y_dens[idx_row])
    obs_dens2d = Observable(dens2d)
    obs_dens2d_hm = lift(ds -> ds', obs_dens2d)
    obs_dens2d_ft_masked = Observable(dens_ft_masked_vec[idx_rep]')
    ft_clr = gen_ft_rgba_clr(dens_ft_amp_vec[idx_rep], dens_ft_phs_vec[idx_rep]; rng_amp=nothing)
    obs_dens2d_ft = Observable(gen_ft_rgba_img(dens_ft_amp_vec[idx_rep], dens_ft_phs_vec[idx_rep]; rng_amp=nothing)')
    obs_dens2d_ft_amp = Observable(dens_ft_amp_vec[idx_rep]')
    obs_dens2d_ft_clr = Observable(ft_clr.colormap')
    ft_cmpx_mean = dens_core_ft_cmpx_mean[ib, istp]
    ft_cmpx_mean_clr = gen_ft_rgba_clr(ft_cmpx_mean)
    obs_dens2d_ft_cmpx_mean = Observable(gen_ft_rgba_img(ft_cmpx_mean)')
    obs_dens2d_ft_cmpx_mean_amp = Observable(abs.(ft_cmpx_mean)')
    obs_dens2d_ft_cmpx_mean_clr = Observable(ft_cmpx_mean_clr.colormap')
    obs_dens2d_ft_cmpx_mean_clr_amp = Observable(collect(ft_cmpx_mean_clr.amp))
    obs_dens2d_ft_cmpx_mean_clr_phs = Observable(collect(ft_cmpx_mean_clr.phs))
    obs_dens2d_ft_cmpx_mean_clr_amp_max = Observable(maximum(ft_cmpx_mean_clr.amp))
    obs_dens2d_ft_absl_mean = Observable(dens_core_ft_absl_mean[ib, istp]')
    ft_absl_cb_x = [0.0, 1.0]
    ft_absl_cb_y = collect(range(clrrng_ft2d_absl_mean...; length=256))
    ft_absl_cb = [y for _ in ft_absl_cb_x, y in ft_absl_cb_y]
    obs_ft_cmpx_visible = lift(view -> view == :complex, obs_ft_cmpx_view)
    obs_ft_amp_visible = lift(view -> view == :amp, obs_ft_cmpx_view)
    obs_colorrange = Observable((0.0, maximum(dens2d)))
    obs_clrmap = Observable(gen_theme_clrmap(istp))
    obs_clr_theme = Observable(gen_theme_clr(istp, 0.3))
    obs_profile_row = Observable(vec(@view dens2d[idx_row, :]))
    obs_profile_row_mean = Observable(vec(mean(@view(dens2d[idxs_center, :]); dims=1)))
    obs_profile_modl = Observable(fit_info.profile_modl[mask_fit_peak_plot])
    obs_fit_gauss = Observable(fit_info.fit_gauss.fit[mask_fit_peak_plot])
    obs_fit_modl = Observable(fit_info.fit_modl.fit[mask_fit_peak_plot])
    obs_x0_peak = Observable(gen_x0_peak(ib, istp))
    obs_x0_peak_current = Observable(fit_info.fit_gauss.params.x0)
    payload_fit_polar = gen_fit_polar_payload(ib, istp, idx_rep, hue_scheme, phase_mode)
    obs_fit_theta = Observable(payload_fit_polar.theta)
    obs_fit_eta = Observable(payload_fit_polar.radius)
    obs_fit_outer_radius = Observable(payload_fit_polar.radius_outer)
    obs_fit_color = Observable(payload_fit_polar.color)
    obs_fit_markersize = Observable(payload_fit_polar.markersize)
    payload_moment = gen_moment_payload(ib, istp)
    obs_moment_theta = Observable(payload_moment.theta)
    obs_moment_radius = Observable(payload_moment.radius)
    obs_fit_gauss_text = Observable(gen_fit_gauss_text(fit_info))
    obs_fit_modl_text_left = Observable(gen_fit_modl_text_left(fit_info))
    obs_fit_modl_text_right = Observable(gen_fit_modl_text_right(fit_info))
    obs_fit_modl_text_rss = Observable(gen_fit_modl_text_rss(fit_info))
    obs_title = lift(obs_idx_IB, obs_idx_istp, obs_idx_rep, obs_val_row) do idx_IB_live, idx_istp_live, idx_rep_live, val_row_live
        @sprintf(
            "IB idx=%d, istp=%s, rep=%d/%d, y_row=%.3f μm, strip_y=%.3f μm",
            idx_IB_live,
            string(val_istp[idx_istp_live]),
            idx_rep_live,
            length(dens_core[idx_IB_live, idx_istp_live]),
            val_row_live,
            y_dens[idx_strip_center],
        )
    end
    obs_title_row = lift(obs_val_row) do val_row_live
        @sprintf("y_row=%.3f μm", val_row_live)
    end
    obs_title_fit_polar = lift(obs_idx_IB, obs_idx_istp, obs_hue_scheme, obs_phase_mode) do idx_IB_live, idx_istp_live, hue_scheme_live, phase_mode_live
        @sprintf("fit %s, M: IB=%.3f, istp=%s, hue=%s", string(phase_mode_live), val_IB[idx_IB_live], string(val_istp[idx_istp_live]), string(hue_scheme_live))
    end
    obs_title_x0_peak = lift(obs_idx_IB, obs_idx_istp) do idx_IB_live, idx_istp_live
        @sprintf("x0 peak: IB=%.3f, istp=%s", val_IB[idx_IB_live], string(val_istp[idx_istp_live]))
    end

    Label(fig[0, 1:2]; text=obs_title, tellwidth=false, halign=:left)

    ax_hm = Axis(
        fig[1, 1];
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect(),
        xgridvisible=true,
        ygridvisible=true,
    )
    try
        deregister_interaction!(ax_hm, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    hspan!(ax_hm, y_strip_min, y_strip_max; color=clr_strip)
    vspan!(ax_hm, x_fit_peak_min, x_fit_peak_max; color=clr_fit_peak_span)
    vspan!(ax_hm, x_fit_modl_min, x_fit_modl_max; color=clr_fit_modl_span)
    hm = heatmap!(ax_hm, x_dens, y_dens, obs_dens2d_hm; colormap=obs_clrmap, colorrange=obs_colorrange, rasterize=true)
    hlines!(ax_hm, lift(x -> [x], obs_val_row); color=obs_clr_theme, linewidth=0.9)

    gl_ft_live = GridLayout(fig[2, 1:2])
    ax_ft_masked = Axis(
        gl_ft_live[1, 1];
        xlabel="x (μm)",
        ylabel="y (μm)",
        title="masked dens",
        aspect=DataAspect(),
        titlesize=11,
    )
    heatmap!(
        ax_ft_masked,
        x_dens_ft,
        y_dens_ft,
        obs_dens2d_ft_masked;
        colormap=obs_clrmap,
        colorrange=obs_colorrange,
        rasterize=true,
    )

    ax_ft = Axis(
        gl_ft_live[1, 2];
        xlabel="kx (μm⁻¹)",
        ylabel="ky (μm⁻¹)",
        title="shot FT",
        aspect=DataAspect(),
        titlesize=11,
    )
    heatmap!(
        ax_ft,
        kx_ft,
        ky_ft,
        obs_dens2d_ft;
        visible=obs_ft_cmpx_visible,
        rasterize=true,
    )
    heatmap!(
        ax_ft,
        kx_ft,
        ky_ft,
        obs_dens2d_ft_amp;
        colormap=obs_clrmap,
        colorrange=clrrng_ft2d_absl_mean,
        visible=obs_ft_amp_visible,
        rasterize=true,
    )
    xlims!(ax_ft, 0, 0.5)
    ylims!(ax_ft, -0.2, 0.2)

    ax_ft_cmpx_mean = Axis(
        gl_ft_live[1, 3];
        xlabel="kx (μm⁻¹)",
        ylabel="ky (μm⁻¹)",
        title="mean FT",
        aspect=DataAspect(),
        titlesize=11,
    )
    heatmap!(
        ax_ft_cmpx_mean,
        kx_ft,
        ky_ft,
        obs_dens2d_ft_cmpx_mean;
        visible=obs_ft_cmpx_visible,
        rasterize=true,
    )
    heatmap!(
        ax_ft_cmpx_mean,
        kx_ft,
        ky_ft,
        obs_dens2d_ft_cmpx_mean_amp;
        colormap=obs_clrmap,
        colorrange=clrrng_ft2d_absl_mean,
        visible=obs_ft_amp_visible,
        rasterize=true,
    )
    xlims!(ax_ft_cmpx_mean, 0, 0.5)
    ylims!(ax_ft_cmpx_mean, -0.2, 0.2)

    ax_ft_cmpx_clr = Axis(
        gl_ft_live[1, 4];
        xlabel="phase",
        ylabel="amp",
        title="FT shader",
        titlesize=11,
    )
    heatmap!(
        ax_ft_cmpx_clr,
        obs_dens2d_ft_cmpx_mean_clr_phs,
        obs_dens2d_ft_cmpx_mean_clr_amp,
        obs_dens2d_ft_cmpx_mean_clr;
        visible=obs_ft_cmpx_visible,
        rasterize=true,
    )
    heatmap!(
        ax_ft_cmpx_clr,
        ft_absl_cb_x,
        ft_absl_cb_y,
        ft_absl_cb;
        colormap=obs_clrmap,
        colorrange=clrrng_ft2d_absl_mean,
        visible=obs_ft_amp_visible,
        rasterize=true,
    )
    xlims!(ax_ft_cmpx_clr, -π, π)
    ylims!(ax_ft_cmpx_clr, clrrng_ft2d_cmpx_mean)
    ft_cmpx_clr_ylim_handler = on(obs_dens2d_ft_cmpx_mean_clr_amp_max) do amp_max_live
        ylims!(ax_ft_cmpx_clr, 0, amp_max_live)
    end

    ax_ft_absl_mean = Axis(
        gl_ft_live[1, 5];
        xlabel="kx (μm⁻¹)",
        ylabel="ky (μm⁻¹)",
        title="mean |FT|",
        aspect=DataAspect(),
        titlesize=11,
    )
    heatmap!(
        ax_ft_absl_mean,
        kx_ft,
        ky_ft,
        obs_dens2d_ft_absl_mean;
        colormap=obs_clrmap,
        colorrange=clrrng_ft2d_absl_mean,
        rasterize=true,
    )
    xlims!(ax_ft_absl_mean, 0, 0.5)
    ylims!(ax_ft_absl_mean, -0.2, 0.2)

    ax_ft_absl_clr = Axis(
        gl_ft_live[1, 6];
        xlabel="",
        ylabel="|FT|",
        title="|FT| scale",
        titlesize=11,
    )
    heatmap!(
        ax_ft_absl_clr,
        ft_absl_cb_x,
        ft_absl_cb_y,
        ft_absl_cb;
        colormap=obs_clrmap,
        colorrange=clrrng_ft2d_absl_mean,
        rasterize=true,
    )
    hidexdecorations!(ax_ft_absl_clr; grid=false)
    colgap!(gl_ft_live, 6)
    rowgap!(gl_ft_live, 0)
    colsize!(gl_ft_live, 4, Fixed(150))
    colsize!(gl_ft_live, 6, Fixed(52))

    ax_row = Axis(
        fig[3, 1];
        xlabel="x (μm)",
        ylabel="density",
        title=obs_title_row,
    )
    try
        deregister_interaction!(ax_row, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    vspan!(ax_row, x_fit_modl_min, x_fit_modl_max; color=clr_fit_modl_span)
    lines!(ax_row, x_dens, obs_profile_row_mean; color=clr_mean, linewidth=2.5)
    lines!(ax_row, x_dens, obs_profile_row; color=obs_clr_theme, linewidth=1.7)
    lines!(ax_row, x_dens[mask_fit_peak_plot], obs_fit_gauss; color=clr_fit, linewidth=1.0)
    text!(
        ax_row,
        0.98,
        0.96;
        text=obs_fit_gauss_text,
        space=:relative,
        align=(:right, :top),
        color=clr_fit,
        font="Consolas",
        fontsize=13,
    )
    xlims!(ax_row, extrema(x_dens))
    ylims!(ax_row, ylims_profile)

    ax_modl = Axis(
        fig[4, 1];
        xlabel="x (μm)",
        ylabel="profile - gaussian",
        title="modulation residual",
    )
    try
        deregister_interaction!(ax_modl, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    vspan!(ax_modl, x_fit_modl_min, x_fit_modl_max; color=clr_fit_modl_span)
    lines!(ax_modl, x_dens[mask_fit_peak_plot], obs_profile_modl; color=clr_mean, linewidth=2.5)
    lines!(ax_modl, x_dens[mask_fit_peak_plot], obs_fit_modl; color=clr_fit, linewidth=1.0)
    vlines!(ax_modl, lift(x0 -> [x0], obs_x0_peak_current); color=clr_fit, linewidth=1.0, linestyle=:dash)
    text!(
        ax_modl,
        0.02,
        0.96;
        text=obs_fit_modl_text_left,
        space=:relative,
        align=(:left, :top),
        color=clr_fit,
        font="Consolas",
        fontsize=13,
    )
    text!(
        ax_modl,
        0.98,
        0.96;
        text=obs_fit_modl_text_right,
        space=:relative,
        align=(:right, :top),
        color=clr_fit,
        font="Consolas",
        fontsize=13,
    )
    text!(
        ax_modl,
        0.98,
        0.04;
        text=obs_fit_modl_text_rss,
        space=:relative,
        align=(:right, :bottom),
        color=clr_fit,
        font="Consolas",
        fontsize=13,
    )
    xlims!(ax_modl, extrema(x_dens[mask_fit_peak_plot]))
    ylims!(ax_modl, (-5, 5))

    ax_fit_polar = PolarAxis(
        fig[1, 2];
        title=obs_title_fit_polar,
        thetaticklabelsize=9,
        rticklabelsize=9,
        rticks=polar_rticks,
    )

    scatter!(
        ax_fit_polar,
        obs_fit_theta,
        obs_fit_eta;
        color=obs_fit_color,
        markersize=obs_fit_markersize,
        strokecolor=(:black, 0.40),
        strokewidth=0.35,
    )
    scatter!(
        ax_fit_polar,
        obs_fit_theta,
        obs_fit_outer_radius;
        color=obs_fit_color,
        markersize=obs_fit_markersize,
        strokecolor=(:black, 0.40),
        strokewidth=0.35,
    )
    scatter!(
        ax_fit_polar,
        obs_moment_theta,
        obs_moment_radius;
        color=clr_moment,
        marker=:diamond,
        markersize=13,
        strokecolor=:white,
        strokewidth=0.8,
    )
    rlims!(ax_fit_polar, 0, polar_radius_limit)

    ax_center = Axis(
        fig[4, 2];
        xlabel="rep",
        ylabel="x0 (μm)",
        title=obs_title_x0_peak,
    )
    lines!(ax_center, 1:n_rep_profile, obs_x0_peak; color=clr_mean, linewidth=1.5)
    vlines!(ax_center, lift(idx -> [idx], obs_idx_rep); color=clr_fit, linewidth=1.2)
    xlims!(ax_center, 1, n_rep_profile)

    function update_profiles!()
        dens_vec_live = dens_core[obs_idx_IB[], obs_idx_istp[]]
        dens_ft_masked_vec_live = dens_core_ft_masked[obs_idx_IB[], obs_idx_istp[]]
        dens_ft_amp_vec_live = dens_core_ft_amp[obs_idx_IB[], obs_idx_istp[]]
        dens_ft_phs_vec_live = dens_core_ft_phs[obs_idx_IB[], obs_idx_istp[]]
        ft_cmpx_mean_live = dens_core_ft_cmpx_mean[obs_idx_IB[], obs_idx_istp[]]
        ft_cmpx_mean_clr_live = gen_ft_rgba_clr(ft_cmpx_mean_live)
        isempty(dens_vec_live) && return nothing
        obs_idx_rep[] = mod1(obs_idx_rep[], length(dens_vec_live))
        dens2d_live = dens_vec_live[obs_idx_rep[]]
        dens2d_ft_masked_live = dens_ft_masked_vec_live[obs_idx_rep[]]
        dens2d_ft_amp_live = dens_ft_amp_vec_live[obs_idx_rep[]]
        dens2d_ft_phs_live = dens_ft_phs_vec_live[obs_idx_rep[]]
        fit_info_live = fit_peak[obs_idx_IB[], obs_idx_istp[]][obs_idx_rep[]]
        obs_dens2d[] = dens2d_live
        obs_dens2d_ft_masked[] = dens2d_ft_masked_live'
        obs_dens2d_ft[] = gen_ft_rgba_img(dens2d_ft_amp_live, dens2d_ft_phs_live; rng_amp=nothing)'
        obs_dens2d_ft_amp[] = dens2d_ft_amp_live'
        obs_dens2d_ft_cmpx_mean[] = gen_ft_rgba_img(ft_cmpx_mean_live)'
        obs_dens2d_ft_cmpx_mean_amp[] = abs.(ft_cmpx_mean_live)'
        obs_dens2d_ft_cmpx_mean_clr[] = ft_cmpx_mean_clr_live.colormap'
        obs_dens2d_ft_cmpx_mean_clr_amp[] = collect(ft_cmpx_mean_clr_live.amp)
        obs_dens2d_ft_cmpx_mean_clr_phs[] = collect(ft_cmpx_mean_clr_live.phs)
        obs_dens2d_ft_cmpx_mean_clr_amp_max[] = maximum(ft_cmpx_mean_clr_live.amp)
        obs_dens2d_ft_absl_mean[] = dens_core_ft_absl_mean[obs_idx_IB[], obs_idx_istp[]]'
        obs_colorrange[] = (0.0, maximum(dens2d_live))
        obs_clrmap[] = gen_theme_clrmap(obs_idx_istp[])
        obs_clr_theme[] = gen_theme_clr(obs_idx_istp[], 0.3)
        obs_profile_row[] = vec(@view dens2d_live[obs_idx_row[], :])
        obs_profile_row_mean[] = vec(mean(@view(dens2d_live[idxs_center, :]); dims=1))
        obs_profile_modl[] = fit_info_live.profile_modl[mask_fit_peak_plot]
        obs_fit_gauss[] = fit_info_live.fit_gauss.fit[mask_fit_peak_plot]
        obs_fit_modl[] = fit_info_live.fit_modl.fit[mask_fit_peak_plot]
        obs_x0_peak[] = gen_x0_peak(obs_idx_IB[], obs_idx_istp[])
        obs_x0_peak_current[] = fit_info_live.fit_gauss.params.x0
        payload_fit_polar_live = gen_fit_polar_payload(obs_idx_IB[], obs_idx_istp[], obs_idx_rep[], obs_hue_scheme[], obs_phase_mode[])
        obs_fit_theta[] = payload_fit_polar_live.theta
        obs_fit_eta[] = payload_fit_polar_live.radius
        obs_fit_outer_radius[] = payload_fit_polar_live.radius_outer
        obs_fit_color[] = payload_fit_polar_live.color
        obs_fit_markersize[] = payload_fit_polar_live.markersize
        payload_moment_live = gen_moment_payload(obs_idx_IB[], obs_idx_istp[])
        obs_moment_theta[] = payload_moment_live.theta
        obs_moment_radius[] = payload_moment_live.radius
        obs_fit_gauss_text[] = gen_fit_gauss_text(fit_info_live)
        obs_fit_modl_text_left[] = gen_fit_modl_text_left(fit_info_live)
        obs_fit_modl_text_right[] = gen_fit_modl_text_right(fit_info_live)
        obs_fit_modl_text_rss[] = gen_fit_modl_text_rss(fit_info_live)
        return nothing
    end

    function update_cut_profiles!(x_click::Real, y_click::Real)
        idx_row_live = argmin(abs.(y_dens .- y_click))
        obs_idx_row[] = idx_row_live
        obs_val_row[] = y_dens[idx_row_live]
        update_profiles!()
        return nothing
    end

    function update_data_index!(step_IB::Integer, step_istp::Integer, step_profile::Integer)
        obs_idx_IB[] = mod1(obs_idx_IB[] + step_IB, size(dens_core, 1))
        obs_idx_istp[] = mod1(obs_idx_istp[] + step_istp, size(dens_core, 2))
        dens_vec_live = dens_core[obs_idx_IB[], obs_idx_istp[]]
        isempty(dens_vec_live) && return nothing
        obs_idx_rep[] = mod1(obs_idx_rep[] + step_profile, length(dens_vec_live))
        update_profiles!()
        return nothing
    end

    function cycle_hue_scheme!()
        obs_hue_scheme[] = obs_hue_scheme[] == :lambda ? :rep : :lambda
        update_profiles!()
        return nothing
    end

    function cycle_phase_mode!()
        obs_phase_mode[] = obs_phase_mode[] == :phi0 ? :phip : :phi0
        update_profiles!()
        return nothing
    end

    function cycle_ft_cmpx_view!()
        obs_ft_cmpx_view[] = obs_ft_cmpx_view[] == :complex ? :amp : :complex
        return nothing
    end

    click_handler = on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press && is_mouseinside(ax_hm.scene)
            xy_click = mouseposition(ax_hm)
            update_cut_profiles!(xy_click[1], xy_click[2])
        end
        return Consume(false)
    end

    gl_ctrl = GridLayout(fig[3, 2])
    gl_ctrl_cycle = GridLayout(gl_ctrl[1, 1])
    gl_ctrl_tune = GridLayout(gl_ctrl[1, 2])
    labels = ("IB", "istp", "rep")
    steps = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    button_handlers = map(enumerate(labels)) do (idx_ctrl, label_ctrl)
        step = steps[idx_ctrl]
        btn_prev = Button(gl_ctrl_cycle[idx_ctrl, 1]; label="←", width=34, height=30)
        Label(gl_ctrl_cycle[idx_ctrl, 2]; text=label_ctrl, tellwidth=true, tellheight=false, halign=:center, valign=:center)
        btn_next = Button(gl_ctrl_cycle[idx_ctrl, 3]; label="→", width=34, height=30)
        (
            on(btn_prev.clicks) do _
                update_data_index!((-step[1]), (-step[2]), (-step[3]))
            end,
            on(btn_next.clicks) do _
                update_data_index!(step...)
            end,
        )
    end
    btn_hue = Button(gl_ctrl_tune[1, 1]; label=lift(s -> "hue: $(s)", obs_hue_scheme), height=30)
    hue_handler = on(btn_hue.clicks) do _
        cycle_hue_scheme!()
    end
    btn_phase = Button(gl_ctrl_tune[2, 1]; label=lift(s -> "phase: $(s)", obs_phase_mode), height=30)
    phase_handler = on(btn_phase.clicks) do _
        cycle_phase_mode!()
    end
    btn_ft_cmpx_view = Button(gl_ctrl_tune[3, 1]; label=lift(s -> "FT: $(s)", obs_ft_cmpx_view), height=30)
    ft_cmpx_view_handler = on(btn_ft_cmpx_view.clicks) do _
        cycle_ft_cmpx_view!()
    end
    Label(
        gl_ctrl_tune[4, 1];
        text=lift(r -> @sprintf("rss ramp %.2f..%.2f", r[1], r[2]), obs_rss_rel_ramp),
        tellwidth=false,
        tellheight=true,
        halign=:center,
    )
    slider_rss_rel_ramp = IntervalSlider(
        gl_ctrl_tune[5, 1];
        range=0.0:0.01:2.0,
        startvalues=obs_rss_rel_ramp[],
    )
    rss_rel_ramp_handler = on(slider_rss_rel_ramp.interval) do interval
        obs_rss_rel_ramp[] = (Float64(interval[1]), Float64(interval[2]))
        update_profiles!()
    end

    colsize!(fig.layout, 1, Fixed(560))
    colsize!(fig.layout, 2, Fixed(420))
    rowsize!(fig.layout, 1, Fixed(380))
    rowsize!(fig.layout, 2, Fixed(150))
    rowsize!(fig.layout, 3, Fixed(130))
    rowsize!(fig.layout, 4, Fixed(220))
    resize_to_layout!(fig)
    return (;
        ax_hm,
        ax_ft_masked,
        ax_ft,
        ax_ft_cmpx_mean,
        ax_ft_cmpx_clr,
        ax_ft_absl_mean,
        ax_ft_absl_clr,
        ax_row,
        ax_modl,
        ax_fit_polar,
        ax_center,
        hm,
        idx_IB=obs_idx_IB,
        idx_istp=obs_idx_istp,
        idx_rep=obs_idx_rep,
        idx_row=obs_idx_row,
        y_row=obs_val_row,
        click_handler,
        button_handlers,
        hue_handler,
        phase_handler,
        ft_cmpx_view_handler,
        rss_rel_ramp_handler,
        ft_cmpx_clr_ylim_handler,
    )
end
##
function draw_phase_distro_table!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    y_dens::AbstractVector{<:Real},
    kx_ft::AbstractVector{<:Real},
    ky_ft::AbstractVector{<:Real},
    val_IB::AbstractVector,
    val_istp::AbstractVector,
    ntfr2d_mean::AbstractMatrix,
    dens_core_ft_cmpx_mean::AbstractMatrix,
    dens_core_ft_absl_mean::AbstractMatrix,
    fit_peak::AbstractMatrix,
    xy_center::AbstractArray{<:Tuple{Int,Int},3};
    x_max_fit_peak::Real,
    x_max_fit_modl::Real,
    x_fit_offset::Real,
    smidx_mean_profile::Integer,
    y_strip_offset::Real,
    x_center_px0::Real,
    lambda_hue_min::Real,
    lambda_hue_max::Real,
    lambda_hue_span::Real,
    polar_lightness::Real,
    polar_chroma::Real,
    polar_alpha::Real,
    polar_lightness_rss_bad::Real,
    polar_chroma_rss_bad::Real,
    polar_alpha_rss_bad::Real,
    rss_rel_ramp::Tuple{<:Real,<:Real},
    clrrng_ft2d_absl_mean::Tuple{<:Real,<:Real},
    clrrng_ft2d_cmpx_mean::Tuple{<:Real,<:Real},
    markersize_fit::Real,
    cohr::AbstractArray{<:NamedTuple},
)
    n_IB = length(val_IB)
    n_istp = length(val_istp)
    n_rep = size(xy_center, 3)
    size(ntfr2d_mean) == (n_IB, n_istp) || throw(DimensionMismatch(
        "ntfr2d_mean size $(size(ntfr2d_mean)) must match (IB, istp) $((n_IB, n_istp)).",
    ))
    size(dens_core_ft_cmpx_mean) == (n_IB, n_istp) || throw(DimensionMismatch(
        "dens_core_ft_cmpx_mean size $(size(dens_core_ft_cmpx_mean)) must match (IB, istp) $((n_IB, n_istp)).",
    ))
    size(dens_core_ft_absl_mean) == (n_IB, n_istp) || throw(DimensionMismatch(
        "dens_core_ft_absl_mean size $(size(dens_core_ft_absl_mean)) must match (IB, istp) $((n_IB, n_istp)).",
    ))
    size(fit_peak) == (n_IB, n_istp) || throw(DimensionMismatch(
        "fit_peak size $(size(fit_peak)) must match (IB, istp) $((n_IB, n_istp)).",
    ))
    size(cohr) == (n_IB, n_istp) || throw(DimensionMismatch(
        "cohr size $(size(cohr)) must match (IB, istp) $((n_IB, n_istp)).",
    ))
    size(xy_center) == (n_IB, n_istp, n_rep) || throw(DimensionMismatch(
        "xy_center size $(size(xy_center)) must match (IB, istp, rep) $((n_IB, n_istp, n_rep)).",
    ))
    for idx in CartesianIndices(dens_core_ft_cmpx_mean)
        size(dens_core_ft_cmpx_mean[idx]) == (length(ky_ft), length(kx_ft)) || throw(DimensionMismatch(
            "dens_core_ft_cmpx_mean[$(Tuple(idx)...)] size $(size(dens_core_ft_cmpx_mean[idx])) must match " *
            "(length(ky_ft), length(kx_ft)) $((length(ky_ft), length(kx_ft))).",
        ))
        size(dens_core_ft_absl_mean[idx]) == (length(ky_ft), length(kx_ft)) || throw(DimensionMismatch(
            "dens_core_ft_absl_mean[$(Tuple(idx)...)] size $(size(dens_core_ft_absl_mean[idx])) must match " *
            "(length(ky_ft), length(kx_ft)) $((length(ky_ft), length(kx_ft))).",
        ))
    end

    step_x = median(diff(x_dens))
    step_y = median(diff(y_dens))
    idx_strip_center = argmin(abs.(y_dens .- y_strip_offset))
    idxs_center = max(1, idx_strip_center - smidx_mean_profile):min(length(y_dens), idx_strip_center + smidx_mean_profile)
    y_strip_min = y_dens[first(idxs_center)] - step_y / 2
    y_strip_max = y_dens[last(idxs_center)] + step_y / 2
    x_fit_peak_min, x_fit_peak_max = (x_fit_offset - x_max_fit_peak, x_fit_offset + x_max_fit_peak)
    x_fit_modl_min, x_fit_modl_max = (x_fit_offset - x_max_fit_modl, x_fit_offset + x_max_fit_modl)
    clr_strip = RGBAf(0.86, 0.86, 0.86, 0.34)
    clr_fit_peak_span = RGBAf(0.86, 0.86, 0.86, 0.14)
    clr_fit_modl_span = RGBAf(0.86, 0.86, 0.86, 0.32)
    colorrange_dens = (0.0, maximum(maximum, ntfr2d_mean))
    radius_max = maximum([
        abs(fit.fit_modl.params.M)
        for fits in fit_peak
        for fit in fits
        if fit.success
    ])
    radius_max = max(radius_max, eps(Float64))
    moment_radius_scale = radius_max
    radius_outer = 1.08 * radius_max
    radius_limit = 1.20 * radius_max
    polar_rticks = (0:2:6, string.(0:2:6))
    center_vals = [
        fit_peak[idx_IB, idx_istp][idx_rep].fit_gauss.params.x0
        for idx_IB in 1:n_IB, idx_istp in 1:n_istp, idx_rep in 1:n_rep
    ]
    center_ylim = extrema(center_vals)
    center_pad = max(0.5, 0.08 * (center_ylim[2] - center_ylim[1]))
    center_ylim = (center_ylim[1] - center_pad, center_ylim[2] + center_pad)

    gen_theme_clr(idx_istp::Integer, alpha::Real) =
        RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
    gen_fit_hue(fit_info, hue_scheme::Symbol) =
        if hue_scheme == :lambda
            lambda_norm = clamp((fit_info.fit_modl.params.λ - lambda_hue_min) / (lambda_hue_max - lambda_hue_min), 0, 1)
            lambda_hue_span * (1 - lambda_norm)
        elseif hue_scheme == :rep
            n_rep > 1 ? 360 * (fit_info.idx_rep - 1) / (n_rep - 1) : 0.0
        else
            throw(ArgumentError("Unknown hue_scheme $hue_scheme."))
        end
    gen_fit_color(fit_info, hue_scheme::Symbol) = begin
        ramp_start, ramp_end = rss_rel_ramp
        rss_norm = clamp(
            (fit_info.fit_modl.rss_rel - ramp_start) / max(ramp_end - ramp_start, eps(Float64)),
            0,
            1,
        )
        lightness = polar_lightness + rss_norm * (polar_lightness_rss_bad - polar_lightness)
        chroma = polar_chroma + rss_norm * (polar_chroma_rss_bad - polar_chroma)
        alpha = polar_alpha + rss_norm * (polar_alpha_rss_bad - polar_alpha)
        RGBAf(Oklch(lightness, chroma, gen_fit_hue(fit_info, hue_scheme)), alpha)
    end
    gen_fit_phase(fit_info, phase_mode::Symbol) =
        if phase_mode == :phi0
            fit_info.fit_modl.params.φ
        elseif phase_mode == :phip
            fit_info.fit_modl.params.φ - 2pi * fit_info.fit_gauss.params.x0 / fit_info.fit_modl.params.λ
        else
            throw(ArgumentError("Unknown phase_mode $phase_mode. Expected :phi0 or :phip."))
        end
    gen_polar_payload(idx_IB::Integer, idx_istp::Integer, hue_scheme::Symbol, phase_mode::Symbol) = begin
        fits = fit_peak[idx_IB, idx_istp]
        ids_success = findall(f -> f.success, fits)
        theta = [mod(gen_fit_phase(fits[idx], phase_mode), 2pi) for idx in ids_success]
        radius = [abs(fits[idx].fit_modl.params.M) for idx in ids_success]
        color = [
            gen_fit_color(fits[idx], hue_scheme)
            for idx in ids_success
        ]
        return (; theta, radius, color)
    end
    gen_moment_payload(idx_IB::Integer, idx_istp::Integer) = begin
        moment = cohr[idx_IB, idx_istp]
        theta = [mod(moment.moment_angel, 2pi)]
        radius = [clamp(moment.moment_length, 0, 1) * moment_radius_scale]
        return (; theta, radius)
    end
    mask_sidepeak = [sqrt(x^2 + y^2) > 0.1 for y in ky_ft, x in kx_ft]
    gen_ft_rgba_table(ft::AbstractMatrix{<:Complex}; rng_amp=clrrng_ft2d_cmpx_mean) = begin
        amp = abs.(ft)
        phs = angle.(ft)
        # amp_max = max(maximum(vec(amp[mask_sidepeak])), eps(Float64))
        amp_max = rng_amp[2] / 2
        rgba = shader_cmpx.(amp ./ amp_max, phs; prescale=ampl_prescaler)
        phs_sample = range(-π, π, 256)
        amp_sample = range(rng_amp..., 256)
        clrmp = [shader_cmpx.(a ./ amp_max, φ; prescale=ampl_prescaler) for a in amp_sample, φ in phs_sample]
        return (rgba, (; amp=collect(amp_sample), phs=collect(phs_sample), colormap=clrmp))
    end
    gen_ft_rgba_table_img(ft::AbstractMatrix{<:Complex}) = first(gen_ft_rgba_table(ft))
    gen_ft_rgba_table_clr(ft::AbstractMatrix{<:Complex}) = last(gen_ft_rgba_table(ft))
    ft_absl_cb_x = [0.0, 1.0]
    ft_absl_cb_y = collect(range(clrrng_ft2d_absl_mean...; length=256))
    ft_absl_cb = [y for _ in ft_absl_cb_x, y in ft_absl_cb_y]

    idx_col_dens = 1
    idx_col_phi0_lambda = 2
    idx_col_phi0_rep = 3
    idx_col_separator = 4
    idx_col_phip_lambda = 5
    idx_col_phip_rep = 6
    idx_col_peak = 7
    idx_col_lambda_hist = 8
    idx_col_ft_cmpx = 9
    idx_col_ft_cmpx_amp = 10
    idx_col_ft_absl = 11
    idx_col_IB_right = 12
    lambda_hist_min, lambda_hist_max = fit_lower_modl[4], fit_upper_modl[4]
    lambda_hist_bins = range(lambda_hist_min, lambda_hist_max; length=24)

    Label(
        fig[0, 1:idx_col_IB_right];
        text=@sprintf(
            "%s phase distro: mean density, φ0/φp fit polar distributions, λ histograms, and x0; peak fit %.1f..%.1f μm, mod fit %.1f..%.1f μm, y strip %.1f..%.1f μm",
            tag,
            x_fit_peak_min,
            x_fit_peak_max,
            x_fit_modl_min,
            x_fit_modl_max,
            y_strip_min,
            y_strip_max,
        ),
        tellwidth=false,
        tellheight=true,
        halign=:left,
    )
    Label(fig[1, idx_col_dens]; text="dens 2D", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_phi0_lambda]; text="φ0 hue λ", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_phi0_rep]; text="φ0 hue rep", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_phip_lambda]; text="φp hue λ", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_phip_rep]; text="φp hue rep", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_peak]; text="peak position", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_lambda_hist]; text="λ hist", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_ft_cmpx]; text="mean FT complex", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_ft_cmpx_amp]; text="mean FT amp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_ft_absl]; text="mean |FT|", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_IB_right]; text="IB", tellwidth=true, tellheight=true, halign=:left, font=:bold)
    Box(
        fig[2:(n_IB + 1), idx_col_separator];
        color=:black,
        strokewidth=0,
    )

    function prep_nested_group!(gl::GridLayout, width_inner::Real; height_label::Real=28, gap::Real=4)
        rowsize!(gl, 1, Fixed(height_label))
        colgap!(gl, gap)
        return gl
    end

    function draw_polar_group!(row::Integer, idx_col::Integer, idx_IB::Integer, IB, phase_mode_tbl::Symbol, hue_scheme::Symbol)
        gl_polar = GridLayout(fig[row, idx_col])
        prep_nested_group!(gl_polar, 340; height_label=34)
        for idx_istp in 1:n_istp
            Label(
                gl_polar[1, idx_istp];
                text=idx_istp == 1 ? @sprintf("IB=%.3f A", IB) : "",
                tellwidth=false,
                tellheight=false,
                halign=:left,
                valign=:bottom,
                fontsize=9,
                color=:black,
            )
            ax_polar = PolarAxis(
                gl_polar[2, idx_istp];
                thetaticklabelsize=7,
                rticklabelsize=7,
                rticks=polar_rticks,
            )
            payload = gen_polar_payload(idx_IB, idx_istp, hue_scheme, phase_mode_tbl)
            scatter!(
                ax_polar,
                payload.theta,
                payload.radius;
                color=payload.color,
                markersize=markersize_fit,
                strokecolor=(:black, 0.36),
                strokewidth=0.25,
            )
            scatter!(
                ax_polar,
                payload.theta,
                fill(radius_outer, length(payload.theta));
                color=payload.color,
                markersize=markersize_fit,
                strokecolor=(:black, 0.36),
                strokewidth=0.25,
            )
            payload_moment = gen_moment_payload(idx_IB, idx_istp)
            scatter!(
                ax_polar,
                payload_moment.theta,
                payload_moment.radius;
                color=RGBAf(0.02, 0.02, 0.02, 0.92),
                marker=:diamond,
                markersize=9,
                strokecolor=:white,
                strokewidth=0.55,
            )
            rlims!(ax_polar, 0, radius_limit)
        end
        return gl_polar
    end

    function draw_ft_group!(row::Integer, idx_col::Integer, idx_IB::Integer, title_left::AbstractString, kind::Symbol)
        gl_ft = GridLayout(fig[row, idx_col])
        prep_nested_group!(gl_ft, 230; height_label=34)
        for idx_istp in 1:n_istp
            Label(
                gl_ft[1, idx_istp];
                text=idx_istp == 1 ? title_left : "",
                tellwidth=false,
                tellheight=false,
                halign=:left,
                valign=:bottom,
                fontsize=9,
                color=:black,
            )
            ax_ft = Axis(
                gl_ft[2, idx_istp];
                xlabel="",
                ylabel="",
                aspect=DataAspect(),
            )
            if kind == :complex
                heatmap!(
                    ax_ft,
                    kx_ft,
                    ky_ft,
                    gen_ft_rgba_table_img(dens_core_ft_cmpx_mean[idx_IB, idx_istp])';
                    rasterize=true,
                )
            elseif kind == :complex_amp
                clrmap = gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
                heatmap!(
                    ax_ft,
                    kx_ft,
                    ky_ft,
                    abs.(dens_core_ft_cmpx_mean[idx_IB, idx_istp])';
                    colormap=clrmap,
                    colorrange=clrrng_ft2d_absl_mean,
                    rasterize=true,
                )
            elseif kind == :absolute
                clrmap = gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
                heatmap!(
                    ax_ft,
                    kx_ft,
                    ky_ft,
                    dens_core_ft_absl_mean[idx_IB, idx_istp]';
                    colormap=clrmap,
                    colorrange=clrrng_ft2d_absl_mean,
                    rasterize=true,
                )
            else
                throw(ArgumentError("Unknown FT table kind $kind."))
            end
            xlims!(ax_ft, extrema(kx_ft))
            ylims!(ax_ft, extrema(ky_ft))
            hideydecorations!(ax_ft; grid=false)
        end
        if kind == :complex
            ft_clr = gen_ft_rgba_table_clr(dens_core_ft_cmpx_mean[idx_IB, 1])
            ax_ft_clr = Axis(
                gl_ft[2, n_istp + 1];
                xlabel="φ",
                ylabel="amp",
                titlesize=8,
            )
            heatmap!(
                ax_ft_clr,
                ft_clr.phs,
                ft_clr.amp,
                ft_clr.colormap';
                rasterize=true,
            )
            xlims!(ax_ft_clr, -π, π)
            ylims!(ax_ft_clr, clrrng_ft2d_cmpx_mean)
            ax_ft_clr.xticks = ([-π, 0, π], ["-π", "0", "+π"])
            colsize!(gl_ft, n_istp + 1, Fixed(80))
        elseif kind == :absolute
            gl_ft_cb = GridLayout(gl_ft[2, n_istp + 1])
            for idx_istp in 1:n_istp
                clrmap = gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
                ax_ft_cb = Axis(
                    gl_ft_cb[1, idx_istp];
                    xlabel="",
                    ylabel=idx_istp == 1 ? "|FT|" : "",
                )
                heatmap!(
                    ax_ft_cb,
                    ft_absl_cb_x,
                    ft_absl_cb_y,
                    ft_absl_cb;
                    colormap=clrmap,
                    colorrange=clrrng_ft2d_absl_mean,
                    rasterize=true,
                )
                hidexdecorations!(ax_ft_cb)
                hideydecorations!(ax_ft_cb; grid=false, label=idx_istp != 1, ticklabels=idx_istp != 1, ticks=idx_istp != 1)
            end
            colgap!(gl_ft_cb, 2)
            colsize!(gl_ft, n_istp + 1, Fixed(78))
        elseif kind == :complex_amp
            colgap!(gl_ft, 4)
        end
        return gl_ft
    end

    for (idx_IB, IB) in enumerate(val_IB)
        row = idx_IB + 1
        is_bottom_row = idx_IB == n_IB
        Label(fig[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:right)
        Label(fig[row, idx_col_IB_right]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:left)

        gl_dens = GridLayout(fig[row, idx_col_dens])
        prep_nested_group!(gl_dens, 360)
        for idx_istp in 1:n_istp
            Label(gl_dens[1, idx_istp]; text=string(val_istp[idx_istp]), tellwidth=false, tellheight=false, halign=:center, valign=:bottom, fontsize=9)
            ax_dens = Axis(
                gl_dens[2, idx_istp];
                xlabel=is_bottom_row ? "x (μm)" : "",
                ylabel=idx_istp == 1 ? "y (μm)" : "",
                aspect=DataAspect(),
            )
            clrmap = gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
            hspan!(ax_dens, y_strip_min, y_strip_max; color=clr_strip)
            vspan!(ax_dens, x_fit_peak_min, x_fit_peak_max; color=clr_fit_peak_span)
            vspan!(ax_dens, x_fit_modl_min, x_fit_modl_max; color=clr_fit_modl_span)
            heatmap!(ax_dens, x_dens, y_dens, ntfr2d_mean[idx_IB, idx_istp]'; colormap=clrmap, colorrange=colorrange_dens, rasterize=true)
            hideydecorations!(ax_dens; label=idx_istp != 1, ticklabels=true, ticks=true, grid=false)

        end

        draw_polar_group!(row, idx_col_phi0_lambda, idx_IB, IB, :phi0, :lambda)
        draw_polar_group!(row, idx_col_phi0_rep, idx_IB, IB, :phi0, :rep)
        draw_polar_group!(row, idx_col_phip_lambda, idx_IB, IB, :phip, :lambda)
        draw_polar_group!(row, idx_col_phip_rep, idx_IB, IB, :phip, :rep)

        gl_peak = GridLayout(fig[row, idx_col_peak])
        prep_nested_group!(gl_peak, 360)
        Label(gl_peak[1, 1]; text=@sprintf("IB=%.3f A", IB), tellwidth=false, tellheight=false, halign=:left, valign=:bottom, fontsize=9)
        ax_peak = Axis(
            gl_peak[2, 1];
            xlabel=is_bottom_row ? "rep" : "",
            ylabel="x0 (μm)",
        )
        for idx_istp in 1:n_istp
            center_x = [fit_peak[idx_IB, idx_istp][idx_rep].fit_gauss.params.x0 for idx_rep in 1:n_rep]
            lines!(ax_peak, 1:n_rep, center_x; color=gen_theme_clr(idx_istp, 0.88), linewidth=1.0)
        end
        xlims!(ax_peak, 1, n_rep)
        ylims!(ax_peak, center_ylim)
        hideydecorations!(ax_peak; label=false, ticklabels=false, ticks=false, grid=false)

        gl_lambda_hist = GridLayout(fig[row, idx_col_lambda_hist])
        prep_nested_group!(gl_lambda_hist, 150)
        Label(gl_lambda_hist[1, 1:2]; text="162 + 164", tellwidth=false, tellheight=false, halign=:center, valign=:bottom, fontsize=9)
        ax_lambda_hist = Axis(gl_lambda_hist[2, 1:2]; xlabel=is_bottom_row ? "λ (μm)" : "", ylabel="count")
        for idx_istp in 1:n_istp
            lambdas = [
                fit.fit_modl.params.λ
                for fit in fit_peak[idx_IB, idx_istp]
                if fit.success
            ]
            hist!(
                ax_lambda_hist,
                lambdas;
                bins=lambda_hist_bins,
                color=gen_theme_clr(idx_istp, 0.38),
                strokecolor=gen_theme_clr(idx_istp, 0.85),
                strokewidth=0.7,
            )
        end
        xlims!(ax_lambda_hist, lambda_hist_min, lambda_hist_max)
        !is_bottom_row && hidexdecorations!(ax_lambda_hist; grid=false)
        hideydecorations!(ax_lambda_hist; label=false, ticklabels=false, ticks=false, grid=false)

        draw_ft_group!(row, idx_col_ft_cmpx, idx_IB, @sprintf("IB=%.3f A", IB), :complex)
        draw_ft_group!(row, idx_col_ft_cmpx_amp, idx_IB, "", :complex_amp)
        draw_ft_group!(row, idx_col_ft_absl, idx_IB, "", :absolute)

        rowsize!(fig.layout, row, Fixed(360))
    end

    colsize!(fig.layout, idx_col_dens, Fixed(730))
    colsize!(fig.layout, idx_col_phi0_lambda, Fixed(690))
    colsize!(fig.layout, idx_col_phi0_rep, Fixed(690))
    colsize!(fig.layout, idx_col_separator, Fixed(2))
    colsize!(fig.layout, idx_col_phip_lambda, Fixed(690))
    colsize!(fig.layout, idx_col_phip_rep, Fixed(690))
    colsize!(fig.layout, idx_col_peak, Fixed(370))
    colsize!(fig.layout, idx_col_lambda_hist, Fixed(310))
    colsize!(fig.layout, idx_col_ft_cmpx, Fixed(560))
    colsize!(fig.layout, idx_col_ft_cmpx_amp, Fixed(470))
    colsize!(fig.layout, idx_col_ft_absl, Fixed(560))
    colsize!(fig.layout, idx_col_IB_right, Fixed(55))
    colgap!(fig.layout, idx_col_dens, 14)
    colgap!(fig.layout, idx_col_phi0_lambda, 10)
    colgap!(fig.layout, idx_col_phi0_rep, 6)
    colgap!(fig.layout, idx_col_separator - 1, 6)
    colgap!(fig.layout, idx_col_separator, 6)
    colgap!(fig.layout, idx_col_phip_lambda, 10)
    colgap!(fig.layout, idx_col_phip_rep, 14)
    colgap!(fig.layout, idx_col_peak, 14)
    colgap!(fig.layout, idx_col_lambda_hist, 14)
    colgap!(fig.layout, idx_col_ft_cmpx, 10)
    colgap!(fig.layout, idx_col_ft_cmpx_amp, 10)
    colgap!(fig.layout, idx_col_ft_absl, 14)
    rowgap!(fig.layout, 0)
    resize_to_layout!(fig)
    return fig
end

##
println("  [$tag] loading densities from $path_data")
dens_raw_fmt = load_density_payload(path_data, val_istp)
n_IB, n_istp, n_rep = size(dens_raw_fmt)
wh_raw = size(dens_raw_fmt[1, 1, 1])
(x_center_px0, y_center_px0) = (wh_raw .+ 1) ./ 2
println("  [$tag] formatted densities as (IB, istp, rep)=$(size(dens_raw_fmt)), image size=$wh_raw")
length(val_IB_ref) == n_IB || throw(DimensionMismatch("val_IB_ref length $(length(val_IB_ref)) must match IB count $n_IB."))
length(val_istp) == n_istp || throw(DimensionMismatch("val_istp length $(length(val_istp)) must match istp count $n_istp."))

step_dens = pixsz * bin / mag
step_ft = 1 ./ (2 .* smwh_dens_ft .* step_dens)
x_dens = step_dens .* collect(-smwh[2]:smwh[2])
y_dens = step_dens .* collect(-smwh[1]:smwh[1])
val_IB = copy(val_IB_ref)
use_common_xy_center in (:free, :fixed_x, :fixed_xy) || throw(ArgumentError(
    "use_common_xy_center must be :free, :fixed_x, or :fixed_xy, got $use_common_xy_center.",
))

num = Array{Float64}(undef, n_IB, n_istp, n_rep)
xy_center_nvlp_px = Array{Tuple{Float64, Float64}}(undef, n_IB, n_istp, n_rep)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp, idx_rep in 1:n_rep
    dens = dens_raw_fmt[idx_IB, idx_istp, idx_rep]
    dens_smooth = imfilter(dens, Kernel.gaussian(sigma_center_filter))

    prfl_x = vec(sum(dens_smooth; dims=1))
    x_fit = collect(1.0:length(prfl_x))
    p0_x = [maximum(prfl_x), (length(prfl_x) + 1) / 2, length(prfl_x) / 10, minimum(prfl_x)]
    x_center = curve_fit(gaussian_offset_1d, x_fit, Float64.(prfl_x), p0_x).param[2]

    prfl_y = vec(sum(dens_smooth; dims=2))
    y_fit = collect(1.0:length(prfl_y))
    p0_y = [maximum(prfl_y), (length(prfl_y) + 1) / 2, length(prfl_y) / 10, minimum(prfl_y)]
    y_center = curve_fit(gaussian_offset_1d, y_fit, Float64.(prfl_y), p0_y).param[2]

    num[idx_IB, idx_istp, idx_rep] = sum(dens)
    xy_center_nvlp_px[idx_IB, idx_istp, idx_rep] = (x_center, y_center)
end
xy_center_shift = Array{Tuple{Int,Int}}(undef, n_IB, n_istp, n_rep)
if use_common_xy_center == :fixed_xy
    for idx_IB in 1:n_IB, idx_istp in 1:n_istp
        x_common = round(Int, mean(first.(xy_center_nvlp_px[idx_IB, idx_istp, :])))
        y_common = round(Int, mean(last.(xy_center_nvlp_px[idx_IB, idx_istp, :])))
        xy_center_shift[idx_IB, idx_istp, :] .= Ref((x_common, y_common))
    end
    println("  [$tag] using fixed xy_center repeated over reps for each IB, istp")
elseif use_common_xy_center == :fixed_x
    for idx_IB in 1:n_IB, idx_istp in 1:n_istp
        x_common = round(Int, mean(first.(xy_center_nvlp_px[idx_IB, idx_istp, :])))
        for idx_rep in 1:n_rep
            _, y_fit = xy_center_nvlp_px[idx_IB, idx_istp, idx_rep]
            xy_center_shift[idx_IB, idx_istp, idx_rep] = (x_common, round(Int, y_fit))
        end
    end
    println("  [$tag] using fixed x center and fitted y centers")
else
    xy_center_shift = map(c -> round.(Int, c), xy_center_nvlp_px)
    println("  [$tag] using fitted xy center")
end

xy_peak_nvlp = map(c -> (c .- x_center_px0) .* step_dens, xy_center_nvlp_px)

mask_valid_duet = trues(n_IB, n_rep)

count_profile_shot = vec(sum(mask_valid_duet; dims=2))
println("  [$tag] using all duet counts per IB=$(count_profile_shot)")

ids_rep_valid = [findall(@view mask_valid_duet[idx_IB, :]) for idx_IB in 1:n_IB]
dens_core = Array{Vector{Matrix{Float64}}}(undef, n_IB, n_istp)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    dens_core[idx_IB, idx_istp] = [
        crop_center(dens_raw_fmt[idx_IB, idx_istp, idx_rep], xy_center_shift[idx_IB, idx_istp, idx_rep], smwh) |> copy
        for idx_rep in 1:n_rep
        if mask_valid_duet[idx_IB, idx_rep]
    ]
end

idx_strip_center = argmin(abs.(y_dens .- y_strip_offset))
idxs_center = max(1, idx_strip_center - smh_dens_strip):min(length(y_dens), idx_strip_center + smh_dens_strip)
mask_fit_peak = abs.(x_dens .- x_fit_offset) .<= x_max_fit_peak
mask_fit_modl = abs.(x_dens .- x_fit_offset) .<= x_max_fit_modl
x_fit_peak = x_dens[mask_fit_peak]
x_fit_modl = x_dens[mask_fit_modl]

gen_freq_shifted(n::Integer, step::Real) =
    iseven(n) ? collect((-n ÷ 2):(n ÷ 2 - 1)) ./ (n * step) :
    collect((-(n ÷ 2)):(n ÷ 2)) ./ (n * step)
idx_ft_x_center = argmin(abs.(x_dens))
idx_ft_y_center = argmin(abs.(y_dens))
idxs_ft_x = (idx_ft_x_center - smwh_dens_ft[1]):(idx_ft_x_center + smwh_dens_ft[1])
idxs_ft_y = (idx_ft_y_center - smwh_dens_ft[2]):(idx_ft_y_center + smwh_dens_ft[2])
first(idxs_ft_x) >= firstindex(x_dens) && last(idxs_ft_x) <= lastindex(x_dens) || throw(ArgumentError(
    "smwh_dens_ft=$smwh_dens_ft exceeds x_dens bounds around x=0.",
))
first(idxs_ft_y) >= firstindex(y_dens) && last(idxs_ft_y) <= lastindex(y_dens) || throw(ArgumentError(
    "smwh_dens_ft=$smwh_dens_ft exceeds y_dens bounds around y=0.",
))
kx_ft_full = gen_freq_shifted(length(idxs_ft_x), step_dens)
ky_ft_full = gen_freq_shifted(length(idxs_ft_y), step_dens)
mask_ft_kx = (kx_ft_full .>= 0) .& (kx_ft_full .<= 0.5)
mask_ft_ky = (ky_ft_full .>= -0.2) .& (ky_ft_full .<= 0.2)
kx_ft = kx_ft_full[mask_ft_kx]
ky_ft = ky_ft_full[mask_ft_ky]
x_dens_ft = x_dens[idxs_ft_x]
y_dens_ft = y_dens[idxs_ft_y]
tukey_dens_ft = tukey1d(smwh_dens_ft[2]; alpha=alpha_tukey[2]) * tukey1d(smwh_dens_ft[1]; alpha=alpha_tukey[1])'

dens_core_ft_masked = map(dens_core) do dens_vec
    map(dens_vec) do dens2d
        @view(dens2d[idxs_ft_y, idxs_ft_x]) .* tukey_dens_ft |> copy
    end
end
dens_core_ft_cmpx = map(dens_core_ft_masked) do dens_vec_masked
    map(dens_vec_masked) do dens_masked
        @pipe dens_masked |>
            ifftshift |> fft |> fftshift |>
            _[mask_ft_ky, mask_ft_kx] |> copy |>
            ft2d -> ft2d ./ (sum(abs.(ft2d)) .* prod(step_ft) ./ 4)
    end
end
dens_core_ft_amp = map(ft_vec -> map(u -> abs.(u), ft_vec), dens_core_ft_cmpx)
dens_core_ft_phs = map(ft_vec -> map(u -> angle.(u), ft_vec), dens_core_ft_cmpx)
dens_core_ft_cmpx_mean = map(dens_core_ft_cmpx) do ft2d_cmpx_vec
    reduce(+, ft2d_cmpx_vec) ./ length(ft2d_cmpx_vec)
end
dens_core_ft_absl_mean = map(dens_core_ft_cmpx) do ft2d_cmpx_vec
    reduce(+, map(ft -> abs.(ft), ft2d_cmpx_vec)) ./ length(ft2d_cmpx_vec)
end
mask_sidepeak = [sqrt((x/0.12)^2 + (y/0.12)^2) > 1 for y in ky_ft, x in kx_ft]
mask_sidepeak = [sqrt(((abs(x)-0.2)/0.08)^2 + (y/0.08)^2) < 1 for y in ky_ft, x in kx_ft]
sum_absl_mean2d = @pipe dens_core_ft_cmpx_mean |> map(u -> abs.(u) |> w -> sum(w[mask_sidepeak]) / sum(w), _)
sum_mean_absl2d = @pipe dens_core_ft_absl_mean |> map(u ->       u |> w -> sum(w[mask_sidepeak]) / sum(w), _)

fit_peak = Array{Vector{NamedTuple}}(undef, n_IB, n_istp)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    fit_peak[idx_IB, idx_istp] = map(enumerate(dens_core[idx_IB, idx_istp])) do (idx_rep_valid, dens2d)
        profile = vec(mean(@view(dens2d[idxs_center, :]); dims=1))
        prfl_strip_mean = Float64.(profile[mask_fit_peak])
        p_init_gauss = [amp_gauss_init, x0_gauss_init, sigma_gauss_init, bg_gauss_init]
        try
            fit_gauss = curve_fit(
                gauss_1d_model,
                x_fit_peak,
                prfl_strip_mean,
                p_init_gauss;
                lower=copy(fit_lower_gauss),
                upper=copy(fit_upper_gauss),
            )
            err_gauss = try
                stderror(fit_gauss)
            catch err
                err isa SingularException || rethrow()
                fill(NaN, length(fit_gauss.param))
            end
            fit_gauss_full = gauss_1d_model(x_dens, fit_gauss.param)
            profile_modl = profile .- fit_gauss_full
            prfl_modl_mean = Float64.(profile_modl[mask_fit_modl])

            fit_trials = NamedTuple[]
            for phi_modl_seed in phi_modl_init
                p_init_modl = [amp_modl_init, slope_modl_init, quad_modl_init, lambda_modl_init, phi_modl_seed]
                try
                    fit_modl = curve_fit(
                        modl_vary_1d_model,
                        x_fit_modl,
                        prfl_modl_mean,
                        p_init_modl;
                        lower=copy(fit_lower_modl),
                        upper=copy(fit_upper_modl),
                    )
                    push!(fit_trials, (; fit=fit_modl, rss=sum(abs2, fit_modl.resid), phi_modl_init=phi_modl_seed))
                catch err
                    @warn "modl_vary_1d_model trial fit failed" idx_IB idx_istp idx_rep=ids_rep_valid[idx_IB][idx_rep_valid] phi_modl_seed err
                end
            end
            isempty(fit_trials) && error("all modulation trial fits failed")
            best_trial = fit_trials[argmin(getfield.(fit_trials, :rss))]
            fit_modl = best_trial.fit
            err_modl = try
                stderror(fit_modl)
            catch err
                err isa SingularException || rethrow()
                fill(NaN, length(fit_modl.param))
            end
            rss_rel_modl = norm(residuals(fit_modl)) / min(norm(prfl_modl_mean), norm(modl_vary_1d_model(x_fit_modl, fit_modl.param)))
            params_gauss = name_gauss_params(fit_gauss.param)
            params_modl = name_modl_params(fit_modl.param)
            (;
                idx_rep=ids_rep_valid[idx_IB][idx_rep_valid],
                success=true,
                profile,
                profile_modl,
                fit_gauss=(;
                    success=true,
                    params=params_gauss,
                    param_err=err_gauss,
                    fit=fit_gauss_full,
                    resid=copy(fit_gauss.resid),
                    rss=sum(abs2, fit_gauss.resid),
                ),
                fit_modl=(;
                    success=true,
                    params=params_modl,
                    param_err=err_modl,
                    fit=modl_vary_1d_model(x_dens, fit_modl.param),
                    resid=copy(fit_modl.resid),
                    rss=best_trial.rss,
                    rss_rel=rss_rel_modl,
                    phi_modl_init=best_trial.phi_modl_init,
                ),
            )
        catch err
            @warn "two-step phase distro fit failed" idx_IB idx_istp idx_rep=ids_rep_valid[idx_IB][idx_rep_valid] err
            (;
                idx_rep=ids_rep_valid[idx_IB][idx_rep_valid],
                success=false,
                profile,
                profile_modl=fill(NaN, length(x_dens)),
                fit_gauss=(;
                    success=false,
                    params=name_gauss_params(fill(NaN, length(fit_lower_gauss))),
                    param_err=fill(NaN, length(fit_lower_gauss)),
                    fit=fill(NaN, length(x_dens)),
                    resid=fill(NaN, length(x_fit_peak)),
                    rss=NaN,
                ),
                fit_modl=(;
                    success=false,
                    params=name_modl_params(fill(NaN, length(fit_lower_modl))),
                    param_err=fill(NaN, length(fit_lower_modl)),
                    fit=fill(NaN, length(x_dens)),
                    resid=fill(NaN, length(x_fit_modl)),
                    rss=NaN,
                    rss_rel=NaN,
                    phi_modl_init=NaN,
                ),
            )
        end
    end
end
count_fit = sum(sum(f.success for f in fits) for fits in fit_peak)
count_fit_err = sum(sum(f.success && (any(isnan, f.fit_gauss.param_err) || any(isnan, f.fit_modl.param_err)) for f in fits) for fits in fit_peak)
println("  [$tag] fitted two-step profiles for $count_fit crops; singular error estimates for $count_fit_err crops")

function calc_circ_moment(φ, weight)
    cmpn = map(φ, weight) do φ, w
        (w * cos(φ), w * sin(φ))
    end
    @pipe cmpn |>
            reduce((a, b) -> a .+ b, _) |>
            m -> (hypot(m...) / sum(weight), atan(reverse(m)...))
end

# cohr = map(CartesianIndices(fit_peak), fit_peak) do ids, fit_reps
#     φ_0 = @pipe fit_reps |> map(f -> f.fit_modl.params.φ, _)
#     M = @pipe fit_reps |> map(f -> f.fit_modl.params.M, _)
#     λ = @pipe fit_reps |> map(f -> f.fit_modl.params.λ, _)
#     x_p = @pipe fit_reps |> map(f -> f.fit_gauss.params.x0, _)
#     x_n = @pipe xy_peak_nvlp[Tuple(ids)..., :] |> map(xy -> xy[1], _)
#     φ_p = @. φ_0 - 2pi * x_p / λ
#     φ_n = @. φ_0 - 2pi * x_n / λ

#     rss_rel = @pipe fit_reps |> map(f -> f.fit_modl.rss, _)

#     ampl_modl = map(M, rss_rel) do M, rss_rel
#         rss_rel > 0.5 ? 0 : M
#     end |> skipmissing |> mean
#     moment_circ = calc_circ_moment(φ_p, Float64.(rss_rel .< 0.5))
#     (; ampl_modl, moment_length=moment_circ[1], moment_angel=moment_circ[2])
# end


# for (i, istp) in enumerate(["162", "164"])
#     lines!(axs_cohr, val_IB_ref, getfield.(cohr[:, i], :moment_length); color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.8))
#     lines!(axs_modl, val_IB_ref, getfield.(cohr[:, i], :ampl_modl);     color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.8))
# end
if !@isdefined(cohr)
    cohr = map(fit_peak) do fit_reps
        ids_good = findall(f -> f.success && isfinite(f.fit_modl.params.φ) && isfinite(f.fit_modl.params.M), fit_reps)
        if isempty(ids_good)
            (; ampl_modl=NaN, moment_length=0.0, moment_angel=0.0)
        else
            φ = [fit_reps[idx].fit_modl.params.φ for idx in ids_good]
            weight = [abs(fit_reps[idx].fit_modl.params.M) for idx in ids_good]
            weight_sum = sum(weight)
            if weight_sum <= eps(Float64)
                (; ampl_modl=0.0, moment_length=0.0, moment_angel=0.0)
            else
                moment_x = sum(w * cos(θ) for (θ, w) in zip(φ, weight))
                moment_y = sum(w * sin(θ) for (θ, w) in zip(φ, weight))
                (;
                    ampl_modl=mean(weight),
                    moment_length=hypot(moment_x, moment_y) / weight_sum,
                    moment_angel=atan(moment_y, moment_x),
                )
            end
        end
    end
end

ntfr2d_mean = map(dens_core) do ds
    isempty(ds) && throw(ArgumentError("No valid densities available for a condition."))
    dropdims(mean(stack(ds); dims=3); dims=3)
end
##

isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
ampl_prescaler = a -> a^2 * 0.4
clrrng_ft2d_absl_mean = (0, 50)
clrrng_ft2d_cmpx_mean = (0, 50)
fit_config = (;
    tag,
    path_data,
    path_output,
    smwh,
    smwh_dens_ft,
    mag,
    pixsz,
    bin,
    sigma_center_filter,
    use_common_xy_center,
    x_max_fit_peak,
    x_max_fit_modl,
    x_fit_offset,
    smh_dens_strip,
    y_strip_offset,
    x_center_px0,
    amp_gauss_init,
    x0_gauss_init,
    sigma_gauss_init,
    bg_gauss_init,
    fit_lower_gauss,
    fit_upper_gauss,
    amp_modl_init,
    slope_modl_init,
    quad_modl_init,
    lambda_modl_init,
    phi_modl_init,
    fit_lower_modl,
    fit_upper_modl,
    polar_lightness,
    polar_chroma,
    polar_alpha,
    polar_lightness_rss_bad,
    polar_chroma_rss_bad,
    polar_alpha_rss_bad,
    rss_rel_ramp,
    clrrng_ft2d_absl_mean,
    phase_mode,
)
JLD2.@save path_fit_jld2 fit_config x_dens y_dens x_dens_ft y_dens_ft kx_ft ky_ft val_IB val_istp num xy_center_nvlp_px xy_center_shift mask_valid_duet ids_rep_valid ntfr2d_mean dens_core_ft_masked dens_core_ft_amp dens_core_ft_phs dens_core_ft_cmpx_mean dens_core_ft_absl_mean fit_peak cohr
println("  [$tag] saved phase distro fit data to $path_fit_jld2")

fig_phase_distro = Figure(fontsize=12)
draw_phase_distro_table!(
    fig_phase_distro,
    x_dens,
    y_dens,
    kx_ft,
    ky_ft,
    val_IB,
    val_istp,
    ntfr2d_mean,
    dens_core_ft_cmpx_mean,
    dens_core_ft_absl_mean,
    fit_peak,
    xy_center_shift;
    x_max_fit_peak,
    x_max_fit_modl,
    x_fit_offset,
    smidx_mean_profile=smh_dens_strip,
    y_strip_offset,
    x_center_px0,
    lambda_hue_min,
    lambda_hue_max,
    lambda_hue_span,
    polar_lightness,
    polar_chroma,
    polar_alpha,
    polar_lightness_rss_bad,
    polar_chroma_rss_bad,
    polar_alpha_rss_bad,
    rss_rel_ramp,
    clrrng_ft2d_absl_mean,
    clrrng_ft2d_cmpx_mean,
    markersize_fit,
    cohr,
)
for ext in ("png", "pdf")
    save(joinpath(path_output, "$filename_plot_phase_distro.$ext"), fig_phase_distro; backend=CairoMakie)
end
println("  [$tag] saved phase distro table to $(joinpath(path_output, "$filename_plot_phase_distro.png"))")

##
fig_live = Figure(fontsize=14)
profile_axes = draw_profile_inspector!(
    fig_live,
    x_dens,
    y_dens,
    dens_core,
    dens_core_ft_masked,
    dens_core_ft_amp,
    dens_core_ft_phs,
    dens_core_ft_cmpx_mean,
    dens_core_ft_absl_mean,
    x_dens_ft,
    y_dens_ft,
    kx_ft,
    ky_ft,
    fit_peak,
    xy_center_shift,
    val_istp;
    ib,
    istp,
    idx_rep,
    y_row,
    smidx_mean_profile=smh_dens_strip,
    y_strip_offset,
    ylims_profile,
    x_max_fit_peak,
    x_max_fit_modl,
    x_fit_offset,
    hue_scheme,
    lambda_hue_min,
    lambda_hue_max,
    lambda_hue_span,
    polar_lightness,
    polar_chroma,
    polar_alpha,
    polar_lightness_rss_bad,
    polar_chroma_rss_bad,
    polar_alpha_rss_bad,
    rss_rel_ramp,
    clrrng_ft2d_absl_mean,
    clrrng_ft2d_cmpx_mean,
    markersize_fit,
    markersize_fit_selected,
    x_center_px0,
    cohr,
)
display(fig_live)

fig_cohr = Figure()
axs_cohr = Axis(fig_cohr[1, 1]; title="coherent 2D modulation spectral weight within mask")
axs_modl = Axis(fig_cohr[2, 1]; title="averaged 2D modulation spectral weight within mask")
for (i, istp) in enumerate(["162", "164"])
    lines!(axs_cohr, val_IB_ref, sum_absl_mean2d[:,i]; color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.8))
    lines!(axs_modl, val_IB_ref, sum_mean_absl2d[:,i]; color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.8))
end
for ext in ("png", "pdf")
    save(joinpath(path_output, "$filename_plot_cohr_modl.$ext"), fig_cohr; backend=CairoMakie)
end
println("  [$tag] saved phase distro table to $(joinpath(path_output, "$filename_plot_cohr_modl.png"))")
