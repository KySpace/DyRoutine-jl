using CairoMakie
using FFTW
using HDF5
using ImageFiltering
using LinearAlgebra
using LsqFit
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
title_anlz = "10.Ntfr2D.Abrr.LinearWeight.ComparePrflModl"
path_data = joinpath(path_root, "0204_interference", "result", "prfl.h5")
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

tag = "SSNTFR"
val_istp = ["162", "164"]
label_x_dens = "position (μm)"
label_x_modl = "wavenum (μm⁻¹)"
r_tail_min_profile = 20.0
range_r_tail_fit = (17.0, 37.0)
fit_center_bound = 12.0
fit_stride_2d = 3
fit_maxiter_2d = 10_000
fit_threshold_log_2d = 1.5e-1
fit_sigma_wide_min = 15.0
model_center = :gaussian
smwh_reconstruct = (150, 150)
xlims_prfl_reconstruct = (0.0, 0.6)
ylims_prfl_reconstruct = (-0.1, 0.6)

function orient_ntfr2d_axes(
    ntfr2d::AbstractArray{<:Real,4},
    x_dens::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
)
    n_x = length(x_dens)
    n_IB = length(val_IB)
    n_istp = length(val_istp)

    size(ntfr2d) == (n_x, n_x, n_istp, n_IB) && return ntfr2d
    size(ntfr2d) == (n_IB, n_istp, n_x, n_x) && return permutedims(ntfr2d, (3, 4, 2, 1))

    throw(DimensionMismatch(
        "ntfr2d_mean size $(size(ntfr2d)) must be either (x, y, istp, IB) " *
        "$((n_x, n_x, n_istp, n_IB)) or (IB, istp, x, y) $((n_IB, n_istp, n_x, n_x)).",
    ))
end

function orient_prfl_axes(
    prfl::AbstractArray{<:Real,3},
    x_modl::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
)
    n_x = length(x_modl)
    n_IB = length(val_IB)
    n_istp = length(val_istp)

    size(prfl) == (n_x, n_istp, n_IB) && return prfl
    size(prfl) == (n_IB, n_istp, n_x) && return permutedims(prfl, (3, 2, 1))

    throw(DimensionMismatch(
        "profile size $(size(prfl)) must be either (x_modl, istp, IB) " *
        "$((n_x, n_istp, n_IB)) or (IB, istp, x_modl) $((n_IB, n_istp, n_x)).",
    ))
end

function calc_grid_center_profile(dens2d::AbstractMatrix{<:Real})
    idx_center = cld(size(dens2d, 1), 2)
    return vec(@view dens2d[idx_center, :])
end

function calc_symmetric_grid_column_profile(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
)
    profile = calc_grid_center_profile(dens2d)
    idx_center = cld(length(x_dens), 2)
    x_half = abs.(x_dens[idx_center:end])
    profile_half = (profile[idx_center:end] .+ reverse(profile[1:idx_center])) ./ 2
    return x_half, profile_half
end

function gaussian_tail_1d_model(r, params)
    A_wide, σ_wide = params[1:2]
    return @. A_wide * exp(-r^2 / (2 * σ_wide^2))
end

function gaussian_narrow_1d_model(r, params)
    A_narrow, σ_narrow = params
    return @. A_narrow * exp(-r^2 / (2 * σ_narrow^2))
end

function lorentzian_narrow_1d_model(r, params)
    A_narrow, γ_narrow = params
    return @. A_narrow / (1 + (r / γ_narrow)^2)
end

function get_narrow_1d_model(model_center::Symbol)
    model_center == :gaussian && return gaussian_narrow_1d_model
    model_center == :lorentzian && return lorentzian_narrow_1d_model
    throw(ArgumentError("model_center must be :gaussian or :lorentzian, got $model_center."))
end

function fit_two_gaussian_1d_guess(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    r_tail_min::Real,
    model_center::Symbol=:gaussian,
)
    r_profile, profile = calc_symmetric_grid_column_profile(x_dens, dens2d)
    mask_tail = r_profile .> r_tail_min
    any(mask_tail) || throw(ArgumentError("No profile coordinates found above r_tail_min=$r_tail_min."))
    max_r = maximum(r_profile)
    step_r = minimum(diff(r_profile))

    y_tail = Float64.(profile[mask_tail])
    p_init_tail = [max(maximum(y_tail), eps(Float64)), max_r / 2]
    p_lower_tail = [0.0, step_r]
    p_upper_tail = [Inf, max_r * 3]
    fit_tail = curve_fit(
        gaussian_tail_1d_model,
        r_profile[mask_tail],
        y_tail,
        p_init_tail;
        lower=p_lower_tail,
        upper=p_upper_tail,
        maxIter=20_000,
    )
    params_tail = coef(fit_tail)
    wide = gaussian_tail_1d_model(r_profile, params_tail)

    y_narrow = max.(Float64.(profile) .- wide, 0.0)
    mask_narrow = r_profile .<= r_tail_min
    p_init_narrow = [max(maximum(y_narrow[mask_narrow]), eps(Float64)), max(r_tail_min / 4, step_r)]
    p_lower_narrow = [0.0, step_r]
    p_upper_narrow = [Inf, r_tail_min]
    fit_narrow = curve_fit(
        get_narrow_1d_model(model_center),
        r_profile[mask_narrow],
        y_narrow[mask_narrow],
        p_init_narrow;
        lower=p_lower_narrow,
        upper=p_upper_narrow,
        maxIter=20_000,
    )
    params_narrow = coef(fit_narrow)

    return (;
        A_narrow=params_narrow[1],
        σ_narrow=params_narrow[2],
        A_wide=params_tail[1],
        σ_wide=params_tail[2],
        B=0.0,
        r_profile,
        profile,
        wide,
        model_center,
    )
end

function two_gaussian_2d_model(coords, params)
    x0, y0, A_narrow, σ_narrow, A_wide, σ_wide = params[1:6]
    x = @view coords[1, :]
    y = @view coords[2, :]
    r2 = @. (x - x0)^2 + (y - y0)^2
    return @. A_narrow * exp(-r2 / (2 * σ_narrow^2)) +
              A_wide * exp(-r2 / (2 * σ_wide^2))
end

function lorentzian_gaussian_2d_model(coords, params)
    x0, y0, A_narrow, γ_narrow, A_wide, σ_wide = params[1:6]
    x = @view coords[1, :]
    y = @view coords[2, :]
    r2 = @. (x - x0)^2 + (y - y0)^2
    return @. A_narrow / (1 + r2 / γ_narrow^2) +
              A_wide * exp(-r2 / (2 * σ_wide^2))
end

function get_density_2d_model(model_center::Symbol)
    model_center == :gaussian && return two_gaussian_2d_model
    model_center == :lorentzian && return lorentzian_gaussian_2d_model
    throw(ArgumentError("model_center must be :gaussian or :lorentzian, got $model_center."))
end

function double_gaussian_disk_2d_model(coords, params)
    x0, y0, A_narrow, σx_narrow, σy_narrow, A_wide, σ_wide = params
    x = @view coords[1, :]
    y = @view coords[2, :]
    dx2 = @. (x - x0)^2
    dy2 = @. (y - y0)^2
    r2 = @. dx2 + dy2
    return @. A_narrow * exp(-dx2 / (2 * σx_narrow^2) - dy2 / (2 * σy_narrow^2)) +
              A_wide * exp(-r2 / (2 * σ_wide^2))
end

function double_gaussian_disk_2d_model_abrr(coords, params)
    x0, y0, A_narrow, σx_narrow, σy_narrow, A_wide, σ_wide, β = params
    x = @view coords[1, :]
    y = @view coords[2, :]
    dx2 = @. (x - x0)^2
    dy2 = @. (y - y0)^2
    r2 = @. dx2 + dy2
    narrow = @. A_narrow * exp(-dx2 / (2 * σx_narrow^2) - dy2 / (2 * σy_narrow^2))
    tail = @. A_wide * exp(-r2 / (2 * σ_wide^2))
    return @. tail + narrow + β * narrow^2
end

function log_double_gaussian_disk_2d_model(coords, params)
    return log.(max.(double_gaussian_disk_2d_model(coords, params), eps(Float64)))
end

function log_double_gaussian_disk_2d_model_abrr(coords, params)
    return log.(max.(double_gaussian_disk_2d_model_abrr(coords, params), eps(Float64)))
end

function calc_fit_coords_2d(x_dens::AbstractVector{<:Real}, stride::Integer)
    idx = 1:stride:length(x_dens)
    x_fit = Float64.(x_dens[idx])
    coords = Matrix{Float64}(undef, 2, length(x_fit)^2)
    coords[1, :] .= repeat(x_fit; outer=length(x_fit))
    coords[2, :] .= repeat(x_fit; inner=length(x_fit))
    return idx, coords
end

function log_gaussian_rim_model(r2, params)
    log_A, σ = params
    return @. log_A - r2 / (2 * σ^2)
end

function fit_log_gaussian_rim_tail(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    center::Tuple{<:Real,<:Real},
    range_r_tail::Tuple{<:Real,<:Real},
)
    r_tail_lo, r_tail_hi = range_r_tail
    r_tail_lo < r_tail_hi || throw(ArgumentError("range_r_tail lower bound $r_tail_lo must be less than upper bound $r_tail_hi."))
    x0, y0 = center
    x_grid = Float64.(x_dens)
    r2 = [
        (x - x0)^2 + (y - y0)^2
        for y in x_grid
        for x in x_grid
    ]
    z = vec(Float64.(dens2d))
    r = sqrt.(r2)
    mask = @. (r_tail_lo <= r <= r_tail_hi) & isfinite(z) & (z > 0)
    any(mask) || throw(ArgumentError("No positive rim density points found in range_r_tail=$range_r_tail."))

    r2_fit = r2[mask]
    log_z_fit = log.(z[mask])
    step_x = minimum(diff(x_dens))
    max_x = maximum(abs, x_dens)
    p_init = [maximum(log_z_fit), max(max_x / 2, step_x)]
    p_lower = [-100.0, step_x]
    p_upper = [100.0, max_x * 3]
    fit = curve_fit(log_gaussian_rim_model, r2_fit, log_z_fit, p_init; lower=p_lower, upper=p_upper, maxIter=20_000)
    params = coef(fit)
    rss_log_rel = norm(residuals(fit)) / max(norm(log_z_fit), eps(Float64))
    return (; A_wide=exp(params[1]), σ_wide=params[2], params, rss_log_rel)
end

function fit_two_gaussian_2d(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    guess_1d;
    center_bound::Real,
    stride::Integer,
    maxiter::Integer=30_000,
    model_center::Symbol=:gaussian,
)
    idx_fit, coords = calc_fit_coords_2d(x_dens, stride)
    z_fit = vec(Float64.(@view dens2d[idx_fit, idx_fit]))
    mask_fit = @. isfinite(z_fit) & (z_fit > fit_threshold_log_2d)
    any(mask_fit) || throw(ArgumentError("No 2D density points above fit_threshold_log_2d=$fit_threshold_log_2d."))
    coords_fit = @view coords[:, mask_fit]
    z_fit_sel = z_fit[mask_fit]
    step_x = minimum(diff(x_dens))
    max_x = maximum(abs, x_dens)
    p_init = [
        0.0,
        0.0,
        max(guess_1d.A_narrow, eps(Float64)),
        max(guess_1d.σ_narrow, step_x),
        max(guess_1d.σ_narrow, step_x),
        max(guess_1d.A_wide, eps(Float64)),
        max(guess_1d.σ_wide, fit_sigma_wide_min),
        0.01,
    ]
    p_lower = [-center_bound, -center_bound, eps(Float64), step_x, step_x, eps(Float64), fit_sigma_wide_min, 0.0]
    p_upper = [center_bound, center_bound, Inf, r_tail_min_profile, r_tail_min_profile, Inf, max_x * 3, 2.0]
    fit = curve_fit(
        double_gaussian_disk_2d_model_abrr,
        coords_fit,
        z_fit_sel,
        p_init;
        lower=p_lower,
        upper=p_upper,
        maxIter=maxiter,
    )
    params = coef(fit)
    rss_rel = norm(residuals(fit)) / max(norm(z_fit_sel), eps(Float64))
    return (; params, rss_rel, guess_1d, model_center, threshold=fit_threshold_log_2d)
end

function interp1_linear(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::Real)
    (xq < first(x) || xq > last(x)) && return NaN
    idx_hi = searchsortedfirst(x, xq)
    idx_hi == 1 && return Float64(y[1])
    idx_hi > length(x) && return Float64(y[end])
    idx_lo = idx_hi - 1
    t = (xq - x[idx_lo]) / (x[idx_hi] - x[idx_lo])
    return (1 - t) * y[idx_lo] + t * y[idx_hi]
end

function interp2_bilinear(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    xq::Real,
    yq::Real,
)
    (xq < first(x_dens) || xq > last(x_dens) || yq < first(x_dens) || yq > last(x_dens)) && return NaN
    idx_x_hi = searchsortedfirst(x_dens, xq)
    idx_y_hi = searchsortedfirst(x_dens, yq)
    idx_x_hi == 1 && (idx_x_hi = 2)
    idx_y_hi == 1 && (idx_y_hi = 2)
    idx_x_hi > length(x_dens) && (idx_x_hi = length(x_dens))
    idx_y_hi > length(x_dens) && (idx_y_hi = length(x_dens))
    idx_x_lo = idx_x_hi - 1
    idx_y_lo = idx_y_hi - 1
    tx = (xq - x_dens[idx_x_lo]) / (x_dens[idx_x_hi] - x_dens[idx_x_lo])
    ty = (yq - x_dens[idx_y_lo]) / (x_dens[idx_y_hi] - x_dens[idx_y_lo])
    z00 = dens2d[idx_x_lo, idx_y_lo]
    z10 = dens2d[idx_x_hi, idx_y_lo]
    z01 = dens2d[idx_x_lo, idx_y_hi]
    z11 = dens2d[idx_x_hi, idx_y_hi]
    return (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 + (1 - tx) * ty * z01 + tx * ty * z11
end

function calc_centered_cross_profile(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    fit_density;
    axis::Symbol,
)
    x0, y0, A_narrow, σx_narrow, σy_narrow, A_wide, σ_wide = fit_density.params[1:7]
    β = length(fit_density.params) >= 8 ? fit_density.params[8] : 0.0
    s_profile = Float64.(x_dens)
    profile =
        axis == :column ? [interp2_bilinear(x_dens, dens2d, x0, y0 + s) for s in s_profile] :
        axis == :row ? [interp2_bilinear(x_dens, dens2d, x0 + s, y0) for s in s_profile] :
        throw(ArgumentError("axis must be :column or :row, got $axis."))
    σ_narrow = axis == :column ? σy_narrow : σx_narrow
    tail_raw = @. A_wide * exp(-s_profile^2 / (2 * σ_wide^2))
    narrow_raw = @. A_narrow * exp(-s_profile^2 / (2 * σ_narrow^2))
    tail = tail_raw
    narrow = @. narrow_raw + β * narrow_raw^2
    narrow_abrr = narrow
    tailess = profile .- tail
    return (; axis, s_profile, profile, tail, narrow, narrow_abrr, tailess, tail_raw, narrow_raw, β, fit_density)
end

function fit_centered_density_profiles(
    x_dens::AbstractVector{<:Real},
    ntfr2d::AbstractArray{<:Real,4},
    r_tail_min::Real;
    center_bound::Real,
    stride::Integer,
    maxiter::Integer=30_000,
    model_center::Symbol=:gaussian,
)
    get_narrow_1d_model(model_center)
    get_density_2d_model(model_center)
    fit_density = Array{NamedTuple}(undef, size(ntfr2d, 4), size(ntfr2d, 3))
    profile_column = Array{NamedTuple}(undef, size(ntfr2d, 4), size(ntfr2d, 3))
    profile_row = similar(profile_column)

    for idx_IB in axes(ntfr2d, 4), idx_istp in axes(ntfr2d, 3)
        dens2d = @view ntfr2d[:, :, idx_istp, idx_IB]
        guess_1d = fit_two_gaussian_1d_guess(x_dens, dens2d, r_tail_min, model_center)
        fit_2d = fit_two_gaussian_2d(x_dens, dens2d, guess_1d; center_bound, stride, maxiter, model_center)
        fit_density[idx_IB, idx_istp] = fit_2d
        profile_column[idx_IB, idx_istp] = calc_centered_cross_profile(x_dens, dens2d, fit_2d; axis=:column)
        profile_row[idx_IB, idx_istp] = calc_centered_cross_profile(x_dens, dens2d, fit_2d; axis=:row)
    end

    return (; fit_density, profile_fits=(profile_column, profile_row))
end

function draw_folded_branch!(
    ax::Axis,
    s::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    side::Symbol;
    color,
    linewidth::Real,
    linestyle=:solid,
)
    mask_side =
        side == :pos ? s .>= 0 :
        side == :neg ? s .<= 0 :
        throw(ArgumentError("side must be :pos or :neg, got $side."))
    x_branch = abs.(Float64.(s[mask_side]))
    y_branch = Float64.(y[mask_side])
    mask_valid = isfinite.(x_branch) .& isfinite.(y_branch) .& (y_branch .> 0)
    count(mask_valid) >= 2 || return nothing
    order = sortperm(x_branch[mask_valid])
    lines!(
        ax,
        x_branch[mask_valid][order],
        y_branch[mask_valid][order];
        color,
        linewidth,
        linestyle,
    )
    return nothing
end

function draw_density_row!(
    fig::Figure,
    row::Integer,
    x_dens::AbstractVector{<:Real},
    x_modl::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    ntfr2d::AbstractArray{<:Real,4},
    idx_IB::Integer,
    IB::Real;
    colorrange,
    ylims_profile,
    profile_fits,
    fit_density,
    xlims_profile,
    xlims_folded,
    ylims_diag,
    prfl_modl_fit,
    prfl_inco,
    prfl_cohr,
    xlims_prfl,
    ylims_prfl,
    is_bottom_row::Bool=false,
)
    Label(fig[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:right)

    axs_dens = Vector{Axis}(undef, length(val_istp))
    axs_profile = Array{Axis}(undef, length(profile_fits), length(val_istp))
    axs_diag = Vector{Axis}(undef, length(val_istp))
    axs_prfl = Vector{Axis}(undef, length(val_istp))

    for (idx_istp, istp) in enumerate(val_istp)
        ax = Axis(
            fig[row, idx_istp];
            xlabel=is_bottom_row ? label_x_dens : "",
            ylabel=idx_istp == 1 ? label_x_dens : "",
            aspect=DataAspect(),
        )
        axs_dens[idx_istp] = ax
        dens2d = @view ntfr2d[:, :, idx_istp, idx_IB]
        clrmap = gen_clrmap_solo(hue_theme_istp[istp])
        heatmap!(ax, x_dens, x_dens, dens2d; colormap=clrmap, colorrange, rasterize=true)
        x0, y0 = fit_density[idx_IB, idx_istp].params[1:2]
        vlines!(ax, x0; color=(:black, 0.16), linewidth=0.7)
        hlines!(ax, y0; color=(:black, 0.16), linewidth=0.7)
        hidexdecorations!(ax; label=is_bottom_row ? false : true, ticklabels=is_bottom_row ? false : true, ticks=is_bottom_row ? false : true, grid=false)
        hideydecorations!(ax; label=idx_istp == 1 ? false : true, ticklabels=false, ticks=false, grid=false)

        for (idx_fit, profile_fit) in enumerate(profile_fits)
            idx_col = length(val_istp) + (idx_fit - 1) * length(val_istp) + idx_istp
            ax_profile = Axis(
                fig[row, idx_col];
                xlabel=is_bottom_row ? label_x_dens : "",
                ylabel=idx_istp == 1 ? "$(profile_fit[idx_IB, idx_istp].axis)" : "",
                yaxisposition=idx_istp == 1 ? :left : :right,
                xticks=-40:20:40,
            )
            axs_profile[idx_fit, idx_istp] = ax_profile
            profile_data = profile_fit[idx_IB, idx_istp]
            s = profile_data.s_profile
            profile = profile_data.profile
            tail = profile_data.tail
            narrow_raw = profile_data.narrow_raw
            narrow = profile_data.narrow
            tailess = profile_data.tailess
            clr_strong = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.95)
            clr_faint = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.32)
            clr_center = Oklch(0.62, 0.16, 145)
            band!(ax_profile, s, zero.(narrow_raw), narrow_raw; color=(clr_center, 0.30))
            lines!(ax_profile, s, profile; color=clr_faint, linewidth=1.0)
            lines!(ax_profile, s, tail; color=(:gray20, 0.55), linewidth=1.0)
            lines!(ax_profile, s, tailess; color=clr_strong, linewidth=1.8)
            band!(ax_profile, s, zero.(narrow), narrow; color=(clr_center, 0.14))
            vlines!(ax_profile, [-r_tail_min_profile, r_tail_min_profile]; color=(:gray20, 0.28), linewidth=0.7)
            xlims!(ax_profile, xlims_profile)
            ylims!(ax_profile, ylims_profile)
            text!(
                ax_profile,
                xlims_profile[1] + 0.04 * (xlims_profile[2] - xlims_profile[1]),
                ylims_profile[2] - 0.08 * (ylims_profile[2] - ylims_profile[1]);
                text=@sprintf("β=%.3f", profile_data.β),
                color=(clr_center, 0.9),
                fontsize=8,
                align=(:left, :top),
            )
            !is_bottom_row && hidexdecorations!(ax_profile; grid=false)
        end

        idx_col_diag = length(val_istp) + length(profile_fits) * length(val_istp) + idx_istp
        ax_diag = Axis(
            fig[row, idx_col_diag];
            xlabel=is_bottom_row ? label_x_dens : "",
            ylabel=idx_istp == 1 ? "folded log" : "",
            yscale=log10,
            yaxisposition=idx_istp == 1 ? :left : :right,
            xticks=0:20:40,
        )
        axs_diag[idx_istp] = ax_diag
        clr_column = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.38)
        clr_row = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.72)
        clr_tail_fit = RGBAf(Oklch(0.58, 0.17, 145), 0.88)
        profile_column = profile_fits[1][idx_IB, idx_istp]
        profile_row = profile_fits[2][idx_IB, idx_istp]
        draw_folded_branch!(ax_diag, profile_column.s_profile, profile_column.profile, :pos; color=clr_column, linewidth=1.0)
        draw_folded_branch!(ax_diag, profile_column.s_profile, profile_column.profile, :neg; color=clr_column, linewidth=1.0)
        draw_folded_branch!(ax_diag, profile_row.s_profile, profile_row.profile, :pos; color=clr_row, linewidth=1.35)
        draw_folded_branch!(ax_diag, profile_row.s_profile, profile_row.profile, :neg; color=clr_row, linewidth=1.35)
        draw_folded_branch!(ax_diag, profile_row.s_profile, profile_row.tail, :pos; color=clr_tail_fit, linewidth=1.5)
        hlines!(ax_diag, fit_threshold_log_2d; color=(:gray20, 0.45), linewidth=0.8)
        xlims!(ax_diag, xlims_folded)
        ylims!(ax_diag, ylims_diag)
        text!(
            ax_diag,
            xlims_folded[1] + 0.04 * (xlims_folded[2] - xlims_folded[1]),
            ylims_diag[1] * 1.25;
            text=@sprintf("σ=%.1f μm", profile_row.fit_density.params[7]),
            color=clr_tail_fit,
            fontsize=8,
            align=(:left, :bottom),
        )
        !is_bottom_row && hidexdecorations!(ax_diag; grid=false)

        idx_col_prfl = length(val_istp) + length(profile_fits) * length(val_istp) + length(val_istp) + idx_istp
        ax_prfl = Axis(
            fig[row, idx_col_prfl];
            xlabel=is_bottom_row ? label_x_modl : "",
            ylabel=idx_istp == 1 ? "FT profile" : "",
            yaxisposition=idx_istp == 1 ? :left : :right,
            xticks=0:0.2:0.6,
        )
        axs_prfl[idx_istp] = ax_prfl
        clr_theme = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.92)
        clr_theme_faint = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.60)
        clr_fit = RGBAf(Oklch(0.60, 0.17, 145), 0.95)
        lines!(ax_prfl, x_modl, @view(prfl_modl_fit[:, idx_istp, idx_IB]); color=clr_fit, linewidth=1.8)
        lines!(ax_prfl, x_modl, @view(prfl_inco[:, idx_istp, idx_IB]); color=clr_theme_faint, linewidth=1.3, linestyle=:dash)
        lines!(ax_prfl, x_modl, @view(prfl_cohr[:, idx_istp, idx_IB]); color=clr_theme, linewidth=1.5)
        xlims!(ax_prfl, xlims_prfl)
        ylims!(ax_prfl, ylims_prfl)
        !is_bottom_row && hidexdecorations!(ax_prfl; grid=false)
    end

    return axs_dens, axs_profile, axs_diag, axs_prfl
end

function tucky1d(sml; alpha=0.2)
    edge = floor(alpha * sml)
    [abs(x) |> x -> x < edge ? 1.0 : (1 - cos(π * x / sml)) / 2
     for x in -sml:sml]
end

function calc_prfl_modl_1d(dens, smwh; step_modl=1)
    smwh .* 2 .+ 1 == size(dens) || error("expect matching smwh and density")
    smw, smh = smwh
    smw_dens_strip = 20
    smh_modl = 120
    tucky = tucky1d(smh; alpha=0.2)
    idx_strip = (smw_dens_strip |> s -> (-s:1:s) .+ smw .+ 1)
    idx_modl = (smh_modl |> s -> (-s:1:s) .+ smh .+ 1)
    dens_strip = @view dens[:, idx_strip]
    dens_mean = imfilter(dens_strip, Kernel.gaussian(2.5)) |> ds -> vec(mean(ds; dims=2))
    modl_full = abs.(fftshift(fft(dens_mean .* tucky)))
    modl = modl_full[idx_modl]
    return modl / (sum(modl) * step_modl / 2)
end

function reconstruct_density_2d(
    x_dens::AbstractVector{<:Real},
    fit_density::AbstractArray,
)
    coords = Matrix{Float64}(undef, 2, length(x_dens)^2)
    x_grid = Float64.(x_dens)
    coords[1, :] .= repeat(x_grid; outer=length(x_grid))
    coords[2, :] .= repeat(x_grid; inner=length(x_grid))

    dens_fit = Array{Float64}(undef, length(x_dens), length(x_dens), size(fit_density, 2), size(fit_density, 1))
    for idx_IB in axes(fit_density, 1), idx_istp in axes(fit_density, 2)
        dens_fit[:, :, idx_istp, idx_IB] .= reshape(
            double_gaussian_disk_2d_model_abrr(coords, fit_density[idx_IB, idx_istp].params),
            length(x_dens),
            length(x_dens),
        )
    end
    return dens_fit
end

function calc_reconstructed_prfl_modl(
    dens_fit::AbstractArray{<:Real,4},
    smwh::Tuple{<:Integer,<:Integer};
    step_modl::Real,
)
    n_modl = length(calc_prfl_modl_1d(@view(dens_fit[:, :, 1, 1]), smwh; step_modl))
    prfl_modl = Array{Float64}(undef, n_modl, size(dens_fit, 3), size(dens_fit, 4))
    for idx_IB in axes(dens_fit, 4), idx_istp in axes(dens_fit, 3)
        prfl_modl[:, idx_istp, idx_IB] .= calc_prfl_modl_1d(
            @view(dens_fit[:, :, idx_istp, idx_IB]),
            smwh;
            step_modl,
        )
    end
    return prfl_modl
end

x_dens, x_modl, val_IB, ntfr2d_mean, prfl_inco, prfl_cohr = h5open(path_data, "r") do file
    x_dens = read(file["x_dens"])
    x_modl = read(file["x_modl"])
    val_IB = read(file["val_IB"])
    ntfr2d_mean = orient_ntfr2d_axes(read(file["ntfr2d_mean"]), x_dens, val_IB, val_istp)
    prfl_inco = orient_prfl_axes(read(file["prfl_inco"]), x_modl, val_IB, val_istp)
    prfl_cohr = orient_prfl_axes(read(file["prfl_cohr"]), x_modl, val_IB, val_istp)
    return x_dens, x_modl, val_IB, ntfr2d_mean, prfl_inco, prfl_cohr
end
step_modl = median(diff(x_modl))

colorrange_ntfr = (0.0, maximum(ntfr2d_mean))
max_profile = maximum([
    maximum(calc_grid_center_profile(@view ntfr2d_mean[:, :, idx_istp, idx_IB]))
    for idx_istp in axes(ntfr2d_mean, 3), idx_IB in axes(ntfr2d_mean, 4)
])
ylims_profile = (0.0, max_profile * 1.05)
fit_centered = fit_centered_density_profiles(
    x_dens,
    ntfr2d_mean,
    r_tail_min_profile;
    center_bound=fit_center_bound,
    stride=fit_stride_2d,
    maxiter=fit_maxiter_2d,
    model_center,
)
fit_density = fit_centered.fit_density
profile_fits = fit_centered.profile_fits
ntfr2d_fit = reconstruct_density_2d(x_dens, fit_density)
prfl_modl_fit = calc_reconstructed_prfl_modl(ntfr2d_fit, smwh_reconstruct; step_modl)
length(x_modl) == size(prfl_modl_fit, 1) || throw(DimensionMismatch(
    "x_modl length $(length(x_modl)) must match reconstructed profile length $(size(prfl_modl_fit, 1)).",
))
min_tailess_profile = minimum(
    minimum(skipmissing(replace(profile_data.tailess, NaN => missing)))
    for fit in profile_fits
    for profile_data in vec(fit)
)
max_original_profile = maximum(
    maximum(skipmissing(replace(profile_data.profile, NaN => missing)))
    for fit in profile_fits
    for profile_data in vec(fit)
)
xlims_profile = (minimum(x_dens), maximum(x_dens))
val_profile_positive = [
    v
    for fit in profile_fits
    for profile_data in vec(fit)
    for v in profile_data.profile
    if isfinite(v) && v > 0
]
ylims_diag = (1e-2, maximum(val_profile_positive) * 1.1)
xlims_folded = (0.0, maximum(abs, x_dens))
ylims_profile_centered = (
    min(0.0, min_tailess_profile * 1.05),
    max_original_profile * 1.05,
)

fig_ntfr = Figure(fontsize=14)
Label(
    fig_ntfr[0, 1:10];
    text=@sprintf(
        "%s 2D NTFR mean densities, cocenter Gaussian tail + Gaussian peak |> (_ + β_²) fit, mask > %.1g, σ_wide ≥ %.0f μm, common max %.3g, Δk=%.5f",
        tag,
        fit_threshold_log_2d,
        fit_sigma_wide_min,
        colorrange_ntfr[2],
        step_modl,
    ),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_ntfr[1, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
for (idx_fit, profile_fit) in enumerate(profile_fits)
    axis_name = profile_fit[1, 1].axis
    for (idx_istp, istp) in enumerate(val_istp)
        idx_col = 2 + (idx_fit - 1) * length(val_istp) + idx_istp
        Label(
            fig_ntfr[1, idx_col];
            text="$(axis_name) tailess $istp",
            tellwidth=false,
            tellheight=true,
            halign=:center,
            font=:bold,
        )
    end
end
for (idx_istp, istp) in enumerate(val_istp)
    idx_col = length(val_istp) + length(profile_fits) * length(val_istp) + idx_istp
    Label(
        fig_ntfr[1, idx_col];
        text="folded log $istp",
        tellwidth=false,
        tellheight=true,
        halign=:center,
        font=:bold,
    )
end
for (idx_istp, istp) in enumerate(val_istp)
    idx_col = length(val_istp) + length(profile_fits) * length(val_istp) + length(val_istp) + idx_istp
    Label(
        fig_ntfr[1, idx_col];
        text="FT prfl $istp",
        tellwidth=false,
        tellheight=true,
        halign=:center,
        font=:bold,
    )
end

for (idx_IB, IB) in enumerate(val_IB)
    row = idx_IB + 1
    draw_density_row!(
        fig_ntfr,
        row,
        x_dens,
        x_modl,
        val_istp,
        ntfr2d_mean,
        idx_IB,
        IB;
        colorrange=colorrange_ntfr,
        ylims_profile=ylims_profile_centered,
        profile_fits,
        fit_density,
        xlims_profile,
        xlims_folded,
        ylims_diag,
        prfl_modl_fit,
        prfl_inco,
        prfl_cohr,
        xlims_prfl=xlims_prfl_reconstruct,
        ylims_prfl=ylims_prfl_reconstruct,
        is_bottom_row=idx_IB == length(val_IB),
    )
    rowsize!(fig_ntfr.layout, row, Fixed(105))
end

colsize!(fig_ntfr.layout, 1, Fixed(105))
colsize!(fig_ntfr.layout, 2, Fixed(105))
colsize!(fig_ntfr.layout, 3, Fixed(170))
colsize!(fig_ntfr.layout, 4, Fixed(170))
colsize!(fig_ntfr.layout, 5, Fixed(170))
colsize!(fig_ntfr.layout, 6, Fixed(170))
colsize!(fig_ntfr.layout, 7, Fixed(170))
colsize!(fig_ntfr.layout, 8, Fixed(170))
colsize!(fig_ntfr.layout, 9, Fixed(170))
colsize!(fig_ntfr.layout, 10, Fixed(170))
colgap!(fig_ntfr.layout, 1, 8)
colgap!(fig_ntfr.layout, 2, 16)
colgap!(fig_ntfr.layout, 3, 0)
colgap!(fig_ntfr.layout, 4, 14)
colgap!(fig_ntfr.layout, 5, 0)
colgap!(fig_ntfr.layout, 6, 14)
colgap!(fig_ntfr.layout, 7, 0)
colgap!(fig_ntfr.layout, 8, 14)
colgap!(fig_ntfr.layout, 9, 0)
rowgap!(fig_ntfr.layout, 1, 4)

resize_to_layout!(fig_ntfr)
filename_plot_ntfr = "$(tag)_ntfr2d_table"
for ext in ("png", "pdf")
    save(joinpath(path_output, "$filename_plot_ntfr.$ext"), fig_ntfr; backend=CairoMakie)
end
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
println("saved $(joinpath(path_output, "$filename_plot_ntfr.png"))")
