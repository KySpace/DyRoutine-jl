using GLMakie
using HDF5
using ImageFiltering
using JLD2
using LsqFit: curve_fit, stderror
using Printf
using Statistics
import CairoMakie

GLMakie.activate!()

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
path_data = joinpath(path_root, "0204_interference", "result", "data.h5")
path_output = joinpath(path_root, "AnlzRoutine", "36.PhaseDistro")
path_fit_jld2 = joinpath(path_output, "SSNTFR_phase_distro_fit.jld2")

tag = "SSNTFR"
filename_plot_phase_distro = "$(tag)_phase_distro_table"
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
mag = 22.06
pixsz = 6.5
bin = 1
sigma_center_filter = 5
use_common_xy_center = false
x_max_fit = 10 # μm
x_fit_offset = 0.0 # μm
smh_dens_strip = 10
y_strip_offset = 0.0 # μm
amp_gauss_init = 6.0
sigma_gauss_init = 12.0
bg_gauss_init = 2.0
# [amp_gauss, sigma_gauss, bg_gauss]
fit_lower_gauss = [0.0, 10.0, 0.0]
fit_upper_gauss = [25.0, 25.0, 5.0]
amp_modl_init = 1.0
slope_modl_init = 0.0
quad_modl_init = 0.0
lambda_modl_init = 5.0
phi_modl_init = (0.0, pi)
# [amp_modl, slope_modl, quad_modl, lambda_modl, phi_modl]
fit_lower_modl = [0.0, -0.20, -0.050, 3.0, -2pi]
fit_upper_modl = [8.0, 0.20, 0.050, 6.0, 2pi]
lambda_hue_min = 4.0
lambda_hue_max = 6.0
lambda_hue_span = 260.0
polar_lightness = 0.74
polar_chroma = 0.12
markersize_fit = 7
markersize_fit_selected = 13

# live inspector selections
ib, istp, idx_rep = (5, 1, 1)
y_row = 0.0
ylims_profile = (-1.0, 15.0)
hue_scheme = :lambda

function gauss_1d_model(x, p)
    (A, σ, bg) = p
    @. A * exp(-(x/σ)^2) + bg
end

function modl_vary_1d_model(x, p)
    (M, a, b, λ, φ) = p
    @. M * (1 + a * x + b * x^2) * cos(2π * x/λ - φ)
end

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

function draw_profile_inspector!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    y_dens::AbstractVector{<:Real},
    dens_core::AbstractMatrix,
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
    x_max_fit::Real,
    x_fit_offset::Real,
    hue_scheme::Symbol,
    lambda_hue_min::Real,
    lambda_hue_max::Real,
    lambda_hue_span::Real,
    polar_lightness::Real,
    polar_chroma::Real,
    markersize_fit::Real,
    markersize_fit_selected::Real,
    x_center_px0::Real,
)
    ib in axes(dens_core, 1) || throw(ArgumentError("ib must be in $(axes(dens_core, 1)), got $ib."))
    istp in axes(dens_core, 2) || throw(ArgumentError("istp must be in $(axes(dens_core, 2)), got $istp."))
    size(dens_core, 2) == length(val_istp) || throw(DimensionMismatch(
        "dens_core second dimension $(size(dens_core, 2)) must match length(val_istp) $(length(val_istp)).",
    ))
    size(fit_peak) == size(dens_core) || throw(DimensionMismatch(
        "fit_peak size $(size(fit_peak)) must match dens_core size $(size(dens_core)).",
    ))
    size(xy_center) == (size(dens_core, 1), size(dens_core, 2), length(first(dens_core))) || throw(DimensionMismatch(
        "xy_center size $(size(xy_center)) must match (IB, istp, rep) $((size(dens_core, 1), size(dens_core, 2), length(first(dens_core)))).",
    ))
    for idx in CartesianIndices(dens_core)
        isempty(dens_core[idx]) && continue
        length(fit_peak[idx]) == length(dens_core[idx]) || throw(DimensionMismatch(
            "fit_peak[$(Tuple(idx)...)] length $(length(fit_peak[idx])) must match dens_core length $(length(dens_core[idx])).",
        ))
        size(first(dens_core[idx])) == (length(y_dens), length(x_dens)) || throw(DimensionMismatch(
            "dens_core[$(Tuple(idx)...)] crop size $(size(first(dens_core[idx]))) must match " *
            "(length(y_dens), length(x_dens)) $((length(y_dens), length(x_dens))).",
        ))
    end

    dens_vec = dens_core[ib, istp]
    isempty(dens_vec) && throw(ArgumentError("dens_core[$ib, $istp] has no selected crops."))
    n_rep_profile = length(dens_vec)
    idx_rep = mod1(idx_rep, length(dens_vec))
    idx_row = argmin(abs.(y_dens .- y_row))
    idx_strip_center = argmin(abs.(y_dens .- y_strip_offset))
    idxs_center = max(1, idx_strip_center - smidx_mean_profile):min(length(y_dens), idx_strip_center + smidx_mean_profile)
    mask_fit_plot = abs.(x_dens .- x_fit_offset) .<= x_max_fit
    dens2d = dens_vec[idx_rep]
    fit_info = fit_peak[ib, istp][idx_rep]
    dens_mean = mean(dens_vec)

    gen_theme_clr(idx_istp::Integer, alpha::Real) =
        RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
    gen_theme_clrmap(idx_istp::Integer) =
        gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
    gen_hue_fit(fit_info, hue_scheme_live::Symbol) =
        if hue_scheme_live == :lambda
            lambda_norm = clamp((fit_info.fit_modl.params[4] - lambda_hue_min) / (lambda_hue_max - lambda_hue_min), 0, 1)
            lambda_hue_span * (1 - lambda_norm)
        elseif hue_scheme_live == :rep
            n_rep_profile > 1 ? 360 * (fit_info.idx_rep - 1) / (n_rep_profile - 1) : 0.0
        else
            throw(ArgumentError("Unknown hue_scheme $hue_scheme_live. Expected :lambda or :rep."))
        end
    gen_fit_polar_payload(idx_IB::Integer, idx_istp::Integer, idx_rep_selected::Integer, hue_scheme_live::Symbol) = begin
        fits = fit_peak[idx_IB, idx_istp]
        ids_success = findall(f -> f.success, fits)
        theta = [mod(fits[idx].fit_modl.params[5], 2pi) for idx in ids_success]
        radius = [abs(fits[idx].fit_modl.params[1]) for idx in ids_success]
        color = map(ids_success) do idx
            fit_info = fits[idx]
            RGBAf(Oklch(polar_lightness, polar_chroma, gen_hue_fit(fit_info, hue_scheme_live)), 0.92)
        end
        markersize = [idx == idx_rep_selected ? markersize_fit_selected : markersize_fit for idx in ids_success]
        return (; theta, radius, color, markersize)
    end

    clr_mean = RGBAf(0.35, 0.35, 0.35, 0.62)
    clr_strip = RGBAf(0.86, 0.86, 0.86, 0.50)
    clr_fit = RGBAf(Oklch(0.60, 0.17, 145), 0.95)
    step_x = median(diff(x_dens))
    step_y = median(diff(y_dens))
    y_strip_min = y_dens[first(idxs_center)] - step_y / 2
    y_strip_max = y_dens[last(idxs_center)] + step_y / 2
    x_fit_min, x_fit_max = (x_fit_offset - x_max_fit, x_fit_offset + x_max_fit)
    gen_x_center_um(idx_IB::Integer, idx_istp::Integer) =
        [(xy_center[idx_IB, idx_istp, idx][1] - x_center_px0) * step_x for idx in 1:n_rep_profile]

    obs_idx_IB = Observable(ib)
    obs_idx_istp = Observable(istp)
    obs_idx_rep = Observable(idx_rep)
    obs_hue_scheme = Observable(hue_scheme)
    obs_idx_row = Observable(idx_row)
    obs_val_row = Observable(y_dens[idx_row])
    obs_dens2d = Observable(dens2d)
    obs_dens2d_hm = lift(ds -> ds', obs_dens2d)
    obs_colorrange = Observable((0.0, maximum(dens2d)))
    obs_clrmap = Observable(gen_theme_clrmap(istp))
    obs_clr_theme = Observable(gen_theme_clr(istp, 0.3))
    obs_profile_row = Observable(vec(@view dens2d[idx_row, :]))
    obs_profile_row_mean = Observable(vec(mean(@view(dens2d[idxs_center, :]); dims=1)))
    obs_profile_modl = Observable(fit_info.profile_modl[mask_fit_plot])
    obs_fit_gauss = Observable(fit_info.fit_gauss.fit[mask_fit_plot])
    obs_fit_modl = Observable(fit_info.fit_modl.fit[mask_fit_plot])
    obs_xy_center_x = Observable(gen_x_center_um(ib, istp))
    payload_fit_polar = gen_fit_polar_payload(ib, istp, idx_rep, hue_scheme)
    obs_fit_theta = Observable(payload_fit_polar.theta)
    obs_fit_eta = Observable(payload_fit_polar.radius)
    obs_fit_color = Observable(payload_fit_polar.color)
    obs_fit_markersize = Observable(payload_fit_polar.markersize)
    obs_fit_gauss_text = Observable(
        @sprintf(
            "A=%.3g\nσ=%.3g\nbg=%.3g",
            fit_info.fit_gauss.params...,
        ),
    )
    obs_fit_modl_text = Observable(
        @sprintf(
            "M=%.3g\na=%.3g\nb=%.3g\nλ=%.3g\nφ=%.3g",
            fit_info.fit_modl.params...,
        ),
    )
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
    obs_title_fit_polar = lift(obs_idx_IB, obs_idx_istp, obs_hue_scheme) do idx_IB_live, idx_istp_live, hue_scheme_live
        @sprintf("fit φ, M: IB=%.3f, istp=%s, hue=%s", val_IB[idx_IB_live], string(val_istp[idx_istp_live]), string(hue_scheme_live))
    end
    obs_title_center = lift(obs_idx_IB, obs_idx_istp) do idx_IB_live, idx_istp_live
        @sprintf("x center: IB=%.3f, istp=%s", val_IB[idx_IB_live], string(val_istp[idx_istp_live]))
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
    vspan!(ax_hm, x_fit_min, x_fit_max; color=clr_strip)
    hm = heatmap!(ax_hm, x_dens, y_dens, obs_dens2d_hm; colormap=obs_clrmap, colorrange=obs_colorrange, rasterize=true)
    hlines!(ax_hm, lift(x -> [x], obs_val_row); color=obs_clr_theme, linewidth=0.9)

    ax_row = Axis(
        fig[2, 1];
        xlabel="x (μm)",
        ylabel="density",
        title=obs_title_row,
    )
    try
        deregister_interaction!(ax_row, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    lines!(ax_row, x_dens, obs_profile_row_mean; color=clr_mean, linewidth=2.5)
    lines!(ax_row, x_dens, obs_profile_row; color=obs_clr_theme, linewidth=1.7)
    lines!(ax_row, x_dens[mask_fit_plot], obs_fit_gauss; color=clr_fit, linewidth=1.0)
    text!(
        ax_row,
        0.98,
        0.96;
        text=obs_fit_gauss_text,
        space=:relative,
        align=(:right, :top),
        color=clr_fit,
        fontsize=13,
    )
    xlims!(ax_row, extrema(x_dens))
    ylims!(ax_row, ylims_profile)

    ax_modl = Axis(
        fig[3, 1];
        xlabel="x (μm)",
        ylabel="profile - gaussian",
        title="modulation residual",
    )
    try
        deregister_interaction!(ax_modl, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    lines!(ax_modl, x_dens[mask_fit_plot], obs_profile_modl; color=clr_mean, linewidth=2.5)
    lines!(ax_modl, x_dens[mask_fit_plot], obs_fit_modl; color=clr_fit, linewidth=1.0)
    text!(
        ax_modl,
        0.98,
        0.96;
        text=obs_fit_modl_text,
        space=:relative,
        align=(:right, :top),
        color=clr_fit,
        fontsize=13,
    )
    xlims!(ax_modl, extrema(x_dens[mask_fit_plot]))
    ylims!(ax_modl, (-5, 5))

    ax_fit_polar = PolarAxis(
        fig[1, 2];
        title=obs_title_fit_polar,
        thetaticklabelsize=9,
        rticklabelsize=9,
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

    ax_center = Axis(
        fig[3, 2];
        xlabel="rep",
        ylabel="x center (μm)",
        title=obs_title_center,
    )
    lines!(ax_center, 1:n_rep_profile, obs_xy_center_x; color=clr_mean, linewidth=1.5)
    vlines!(ax_center, lift(idx -> [idx], obs_idx_rep); color=clr_fit, linewidth=1.2)
    xlims!(ax_center, 1, n_rep_profile)

    function update_profiles!()
        dens_vec_live = dens_core[obs_idx_IB[], obs_idx_istp[]]
        isempty(dens_vec_live) && return nothing
        obs_idx_rep[] = mod1(obs_idx_rep[], length(dens_vec_live))
        dens2d_live = dens_vec_live[obs_idx_rep[]]
        fit_info_live = fit_peak[obs_idx_IB[], obs_idx_istp[]][obs_idx_rep[]]
        dens_mean_live = mean(dens_vec_live)
        obs_dens2d[] = dens2d_live
        obs_colorrange[] = (0.0, maximum(dens2d_live))
        obs_clrmap[] = gen_theme_clrmap(obs_idx_istp[])
        obs_clr_theme[] = gen_theme_clr(obs_idx_istp[], 0.3)
        obs_profile_row[] = vec(@view dens2d_live[obs_idx_row[], :])
        obs_profile_row_mean[] = vec(mean(@view(dens2d_live[idxs_center, :]); dims=1))
        obs_profile_modl[] = fit_info_live.profile_modl[mask_fit_plot]
        obs_fit_gauss[] = fit_info_live.fit_gauss.fit[mask_fit_plot]
        obs_fit_modl[] = fit_info_live.fit_modl.fit[mask_fit_plot]
        obs_xy_center_x[] = gen_x_center_um(obs_idx_IB[], obs_idx_istp[])
        payload_fit_polar_live = gen_fit_polar_payload(obs_idx_IB[], obs_idx_istp[], obs_idx_rep[], obs_hue_scheme[])
        obs_fit_theta[] = payload_fit_polar_live.theta
        obs_fit_eta[] = payload_fit_polar_live.radius
        obs_fit_color[] = payload_fit_polar_live.color
        obs_fit_markersize[] = payload_fit_polar_live.markersize
        obs_fit_gauss_text[] = @sprintf(
            "A=%.3g\nσ=%.3g\nbg=%.3g",
            fit_info_live.fit_gauss.params...,
        )
        obs_fit_modl_text[] = @sprintf(
            "M=%.3g\na=%.3g\nb=%.3g\nλ=%.3g\nφ=%.3g",
            fit_info_live.fit_modl.params...,
        )
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

    click_handler = on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press && is_mouseinside(ax_hm.scene)
            xy_click = mouseposition(ax_hm)
            update_cut_profiles!(xy_click[1], xy_click[2])
        end
        return Consume(false)
    end

    gl_ctrl = GridLayout(fig[2, 2])
    labels = ("IB", "istp", "rep")
    steps = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    button_handlers = map(enumerate(labels)) do (idx_ctrl, label_ctrl)
        step = steps[idx_ctrl]
        btn_prev = Button(gl_ctrl[idx_ctrl, 1]; label="←", width=34, height=30)
        Label(gl_ctrl[idx_ctrl, 2]; text=label_ctrl, tellwidth=true, tellheight=false, halign=:center, valign=:center)
        btn_next = Button(gl_ctrl[idx_ctrl, 3]; label="→", width=34, height=30)
        (
            on(btn_prev.clicks) do _
                update_data_index!((-step[1]), (-step[2]), (-step[3]))
            end,
            on(btn_next.clicks) do _
                update_data_index!(step...)
            end,
        )
    end
    btn_hue = Button(gl_ctrl[4, 1:3]; label=lift(s -> "hue: $(s)", obs_hue_scheme), height=30)
    hue_handler = on(btn_hue.clicks) do _
        cycle_hue_scheme!()
    end

    colsize!(fig.layout, 1, Fixed(360))
    colsize!(fig.layout, 2, Fixed(300))
    rowsize!(fig.layout, 1, Fixed(360))
    rowsize!(fig.layout, 2, Fixed(260))
    rowsize!(fig.layout, 3, Fixed(220))
    resize_to_layout!(fig)
    return (;
        ax_hm,
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
    )
end

function draw_phase_distro_table!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    y_dens::AbstractVector{<:Real},
    val_IB::AbstractVector,
    val_istp::AbstractVector,
    ntfr2d_mean::AbstractMatrix,
    fit_peak::AbstractMatrix,
    xy_center::AbstractArray{<:Tuple{Int,Int},3};
    x_max_fit::Real,
    x_fit_offset::Real,
    smidx_mean_profile::Integer,
    y_strip_offset::Real,
    x_center_px0::Real,
    lambda_hue_min::Real,
    lambda_hue_max::Real,
    lambda_hue_span::Real,
    polar_lightness::Real,
    polar_chroma::Real,
    markersize_fit::Real,
)
    n_IB = length(val_IB)
    n_istp = length(val_istp)
    n_rep = size(xy_center, 3)
    size(ntfr2d_mean) == (n_IB, n_istp) || throw(DimensionMismatch(
        "ntfr2d_mean size $(size(ntfr2d_mean)) must match (IB, istp) $((n_IB, n_istp)).",
    ))
    size(fit_peak) == (n_IB, n_istp) || throw(DimensionMismatch(
        "fit_peak size $(size(fit_peak)) must match (IB, istp) $((n_IB, n_istp)).",
    ))
    size(xy_center) == (n_IB, n_istp, n_rep) || throw(DimensionMismatch(
        "xy_center size $(size(xy_center)) must match (IB, istp, rep) $((n_IB, n_istp, n_rep)).",
    ))

    step_x = median(diff(x_dens))
    step_y = median(diff(y_dens))
    idx_strip_center = argmin(abs.(y_dens .- y_strip_offset))
    idxs_center = max(1, idx_strip_center - smidx_mean_profile):min(length(y_dens), idx_strip_center + smidx_mean_profile)
    y_strip_min = y_dens[first(idxs_center)] - step_y / 2
    y_strip_max = y_dens[last(idxs_center)] + step_y / 2
    x_fit_min, x_fit_max = (x_fit_offset - x_max_fit, x_fit_offset + x_max_fit)
    clr_strip = RGBAf(0.86, 0.86, 0.86, 0.34)
    colorrange_dens = (0.0, maximum(maximum, ntfr2d_mean))
    radius_max = maximum([
        abs(fit.fit_modl.params[1])
        for fits in fit_peak
        for fit in fits
        if fit.success
    ])
    radius_max = max(radius_max, eps(Float64))
    center_vals = [
        (xy_center[idx_IB, idx_istp, idx_rep][1] - x_center_px0) * step_x
        for idx_IB in 1:n_IB, idx_istp in 1:n_istp, idx_rep in 1:n_rep
    ]
    center_ylim = extrema(center_vals)
    center_pad = max(0.5, 0.08 * (center_ylim[2] - center_ylim[1]))
    center_ylim = (center_ylim[1] - center_pad, center_ylim[2] + center_pad)

    gen_theme_clr(idx_istp::Integer, alpha::Real) =
        RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
    gen_fit_hue(fit_info, hue_scheme::Symbol) =
        if hue_scheme == :lambda
            lambda_norm = clamp((fit_info.fit_modl.params[4] - lambda_hue_min) / (lambda_hue_max - lambda_hue_min), 0, 1)
            lambda_hue_span * (1 - lambda_norm)
        elseif hue_scheme == :rep
            n_rep > 1 ? 360 * (fit_info.idx_rep - 1) / (n_rep - 1) : 0.0
        else
            throw(ArgumentError("Unknown hue_scheme $hue_scheme."))
        end
    gen_polar_payload(idx_IB::Integer, idx_istp::Integer, hue_scheme::Symbol) = begin
        fits = fit_peak[idx_IB, idx_istp]
        ids_success = findall(f -> f.success, fits)
        theta = [mod(fits[idx].fit_modl.params[5], 2pi) for idx in ids_success]
        radius = [abs(fits[idx].fit_modl.params[1]) for idx in ids_success]
        color = [
            RGBAf(Oklch(polar_lightness, polar_chroma, gen_fit_hue(fits[idx], hue_scheme)), 0.92)
            for idx in ids_success
        ]
        return (; theta, radius, color)
    end

    idx_col_center = 3 * n_istp + 1
    idx_col_IB_right = idx_col_center + 1

    Label(
        fig[0, 1:idx_col_IB_right];
        text=@sprintf(
            "%s phase distro: mean density, fit polar distributions, and x center; x fit %.1f..%.1f μm, y strip %.1f..%.1f μm",
            tag,
            x_fit_min,
            x_fit_max,
            y_strip_min,
            y_strip_max,
        ),
        tellwidth=false,
        tellheight=true,
        halign=:left,
    )
    for (idx_istp, istp) in enumerate(val_istp)
        Label(fig[1, idx_istp]; text="mean $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
        Label(fig[1, n_istp + idx_istp]; text="polar λ $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
        Label(fig[1, 2n_istp + idx_istp]; text="polar rep $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    end
    Label(fig[1, idx_col_center]; text="x center", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    Label(fig[1, idx_col_IB_right]; text="IB", tellwidth=true, tellheight=true, halign=:left, font=:bold)

    for (idx_IB, IB) in enumerate(val_IB)
        row = idx_IB + 1
        is_bottom_row = idx_IB == n_IB
        Label(fig[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:right)
        Label(fig[row, idx_col_IB_right]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:left)

        for idx_istp in 1:n_istp
            ax_dens = Axis(
                fig[row, idx_istp];
                xlabel=is_bottom_row ? "x (μm)" : "",
                ylabel=idx_istp == 1 ? "y (μm)" : "",
                aspect=DataAspect(),
            )
            clrmap = gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
            hspan!(ax_dens, y_strip_min, y_strip_max; color=clr_strip)
            vspan!(ax_dens, x_fit_min, x_fit_max; color=clr_strip)
            heatmap!(ax_dens, x_dens, y_dens, ntfr2d_mean[idx_IB, idx_istp]'; colormap=clrmap, colorrange=colorrange_dens, rasterize=true)
            hidexdecorations!(ax_dens; label=!is_bottom_row, ticklabels=!is_bottom_row, ticks=!is_bottom_row, grid=false)
            hideydecorations!(ax_dens; label=idx_istp != 1, ticklabels=true, ticks=true, grid=false)

            for (idx_group, hue_scheme) in enumerate((:lambda, :rep))
                idx_col = idx_group * n_istp + idx_istp
                ax_polar = PolarAxis(
                    fig[row, idx_col];
                    thetaticklabelsize=7,
                    rticklabelsize=7,
                )
                payload = gen_polar_payload(idx_IB, idx_istp, hue_scheme)
                scatter!(
                    ax_polar,
                    payload.theta,
                    payload.radius;
                    color=payload.color,
                    markersize=markersize_fit,
                    strokecolor=(:black, 0.36),
                    strokewidth=0.25,
                )
                rlims!(ax_polar, 0, radius_max)
            end
        end

        ax_center = Axis(
            fig[row, idx_col_center];
            xlabel=is_bottom_row ? "rep" : "",
            ylabel="x center (μm)",
        )
        for idx_istp in 1:n_istp
            center_x = [(xy_center[idx_IB, idx_istp, idx_rep][1] - x_center_px0) * step_x for idx_rep in 1:n_rep]
            lines!(ax_center, 1:n_rep, center_x; color=gen_theme_clr(idx_istp, 0.88), linewidth=1.0)
        end
        xlims!(ax_center, 1, n_rep)
        ylims!(ax_center, center_ylim)
        !is_bottom_row && hidexdecorations!(ax_center; grid=false)
        hideydecorations!(ax_center; label=false, ticklabels=false, ticks=false, grid=false)

        rowsize!(fig.layout, row, Fixed(360))
    end

    for idx_col in 1:(3 * n_istp)
        colsize!(fig.layout, idx_col, Fixed(idx_col <= n_istp ? 360 : 300))
    end
    colsize!(fig.layout, idx_col_center, Fixed(360))
    colsize!(fig.layout, idx_col_IB_right, Fixed(55))
    colgap!(fig.layout, n_istp, 14)
    colgap!(fig.layout, 2n_istp, 14)
    colgap!(fig.layout, 3n_istp, 14)
    rowgap!(fig.layout, 1, 4)
    resize_to_layout!(fig)
    return fig
end

println("  [$tag] loading densities from $path_data")
dens_raw_fmt = load_density_payload(path_data, val_istp)
n_IB, n_istp, n_rep = size(dens_raw_fmt)
wh_raw = size(dens_raw_fmt[1, 1, 1])
x_center_px0 = (wh_raw[2] + 1) / 2
println("  [$tag] formatted densities as (IB, istp, rep)=$(size(dens_raw_fmt)), image size=$wh_raw")
length(val_IB_ref) == n_IB || throw(DimensionMismatch("val_IB_ref length $(length(val_IB_ref)) must match IB count $n_IB."))
length(val_istp) == n_istp || throw(DimensionMismatch("val_istp length $(length(val_istp)) must match istp count $n_istp."))

step_dens = pixsz * bin / mag
x_dens = step_dens .* collect(-smwh[2]:smwh[2])
y_dens = step_dens .* collect(-smwh[1]:smwh[1])
val_IB = copy(val_IB_ref)

num = Array{Float64}(undef, n_IB, n_istp, n_rep)
xy_center = Array{Tuple{Int,Int}}(undef, n_IB, n_istp, n_rep)
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
    xy_center[idx_IB, idx_istp, idx_rep] = round.(Int, (x_center, y_center))
end
if use_common_xy_center
    for idx_IB in 1:n_IB, idx_istp in 1:n_istp
        x_common = round(Int, mean(first.(xy_center[idx_IB, idx_istp, :])))
        y_common = round(Int, mean(last.(xy_center[idx_IB, idx_istp, :])))
        xy_center[idx_IB, idx_istp, :] .= Ref((x_common, y_common))
    end
    println("  [$tag] using common xy_center repeated over reps for each IB, istp")
end

mask_valid_duet = trues(n_IB, n_rep)

count_profile_shot = vec(sum(mask_valid_duet; dims=2))
println("  [$tag] using all duet counts per IB=$(count_profile_shot)")

ids_rep_valid = [findall(@view mask_valid_duet[idx_IB, :]) for idx_IB in 1:n_IB]
dens_core = Array{Vector{Matrix{Float64}}}(undef, n_IB, n_istp)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    dens_core[idx_IB, idx_istp] = [
        crop_center(dens_raw_fmt[idx_IB, idx_istp, idx_rep], xy_center[idx_IB, idx_istp, idx_rep], smwh) |> copy
        for idx_rep in 1:n_rep
        if mask_valid_duet[idx_IB, idx_rep]
    ]
end

idx_strip_center = argmin(abs.(y_dens .- y_strip_offset))
idxs_center = max(1, idx_strip_center - smh_dens_strip):min(length(y_dens), idx_strip_center + smh_dens_strip)
mask_fit = abs.(x_dens .- x_fit_offset) .<= x_max_fit
x_fit_peak = x_dens[mask_fit]

fit_peak = Array{Vector{NamedTuple}}(undef, n_IB, n_istp)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    fit_peak[idx_IB, idx_istp] = map(enumerate(dens_core[idx_IB, idx_istp])) do (idx_rep_valid, dens2d)
        profile = vec(mean(@view(dens2d[idxs_center, :]); dims=1))
        prfl_strip_mean = Float64.(profile[mask_fit])
        p_init_gauss = [amp_gauss_init, sigma_gauss_init, bg_gauss_init]
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
            prfl_modl_mean = Float64.(profile_modl[mask_fit])

            fit_trials = NamedTuple[]
            for phi_modl_seed in phi_modl_init
                p_init_modl = [amp_modl_init, slope_modl_init, quad_modl_init, lambda_modl_init, phi_modl_seed]
                try
                    fit_modl = curve_fit(
                        modl_vary_1d_model,
                        x_fit_peak,
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
            (;
                idx_rep=ids_rep_valid[idx_IB][idx_rep_valid],
                success=true,
                profile,
                profile_modl,
                fit_gauss=(;
                    success=true,
                    params=copy(fit_gauss.param),
                    param_err=err_gauss,
                    fit=fit_gauss_full,
                    resid=copy(fit_gauss.resid),
                    rss=sum(abs2, fit_gauss.resid),
                ),
                fit_modl=(;
                    success=true,
                    params=copy(fit_modl.param),
                    param_err=err_modl,
                    fit=modl_vary_1d_model(x_dens, fit_modl.param),
                    resid=copy(fit_modl.resid),
                    rss=best_trial.rss,
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
                    params=fill(NaN, length(fit_lower_gauss)),
                    param_err=fill(NaN, length(fit_lower_gauss)),
                    fit=fill(NaN, length(x_dens)),
                    resid=fill(NaN, length(x_fit_peak)),
                    rss=NaN,
                ),
                fit_modl=(;
                    success=false,
                    params=fill(NaN, length(fit_lower_modl)),
                    param_err=fill(NaN, length(fit_lower_modl)),
                    fit=fill(NaN, length(x_dens)),
                    resid=fill(NaN, length(x_fit_peak)),
                    rss=NaN,
                    phi_modl_init=NaN,
                ),
            )
        end
    end
end
count_fit = sum(sum(f.success for f in fits) for fits in fit_peak)
count_fit_err = sum(sum(f.success && (any(isnan, f.fit_gauss.param_err) || any(isnan, f.fit_modl.param_err)) for f in fits) for fits in fit_peak)
println("  [$tag] fitted two-step profiles for $count_fit crops; singular error estimates for $count_fit_err crops")

ntfr2d_mean = map(dens_core) do ds
    isempty(ds) && throw(ArgumentError("No valid densities available for a condition."))
    dropdims(mean(stack(ds); dims=3); dims=3)
end

isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
fit_config = (;
    tag,
    path_data,
    path_output,
    smwh,
    mag,
    pixsz,
    bin,
    sigma_center_filter,
    use_common_xy_center,
    x_max_fit,
    x_fit_offset,
    smh_dens_strip,
    y_strip_offset,
    x_center_px0,
    amp_gauss_init,
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
)
JLD2.@save path_fit_jld2 fit_config x_dens y_dens val_IB val_istp num xy_center mask_valid_duet ids_rep_valid ntfr2d_mean fit_peak
println("  [$tag] saved phase distro fit data to $path_fit_jld2")

fig_phase_distro = Figure(fontsize=12)
draw_phase_distro_table!(
    fig_phase_distro,
    x_dens,
    y_dens,
    val_IB,
    val_istp,
    ntfr2d_mean,
    fit_peak,
    xy_center;
    x_max_fit,
    x_fit_offset,
    smidx_mean_profile=smh_dens_strip,
    y_strip_offset,
    x_center_px0,
    lambda_hue_min,
    lambda_hue_max,
    lambda_hue_span,
    polar_lightness,
    polar_chroma,
    markersize_fit,
)
for ext in ("png", "pdf")
    save(joinpath(path_output, "$filename_plot_phase_distro.$ext"), fig_phase_distro; backend=CairoMakie)
end
println("  [$tag] saved phase distro table to $(joinpath(path_output, "$filename_plot_phase_distro.png"))")

fig_live = Figure(fontsize=14)
profile_axes = draw_profile_inspector!(
    fig_live,
    x_dens,
    y_dens,
    dens_core,
    fit_peak,
    xy_center,
    val_istp;
    ib,
    istp,
    idx_rep,
    y_row,
    smidx_mean_profile=smh_dens_strip,
    y_strip_offset,
    ylims_profile,
    x_max_fit,
    x_fit_offset,
    hue_scheme,
    lambda_hue_min,
    lambda_hue_max,
    lambda_hue_span,
    polar_lightness,
    polar_chroma,
    markersize_fit,
    markersize_fit_selected,
    x_center_px0,
)
display(fig_live)
