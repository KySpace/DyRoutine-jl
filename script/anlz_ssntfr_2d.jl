using CairoMakie
using HDF5
using LinearAlgebra
using LsqFit
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
title_anlz = "03.Ntfr2D"
path_data = joinpath(path_root, "0204_interference", "result", "prfl.h5")
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

tag = "SSNTFR"
val_istp = ["162", "164"]
label_x_dens = "position (μm)"
r_tail_min_profile = 20.0

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

function calc_center_profile(dens2d::AbstractMatrix{<:Real})
    idx_center = cld(size(dens2d, 2), 2)
    return vec(@view dens2d[:, idx_center])
end

function calc_center_row_profile(dens2d::AbstractMatrix{<:Real})
    idx_center = cld(size(dens2d, 1), 2)
    return vec(@view dens2d[idx_center, :])
end

function calc_symmetric_profile(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real};
    axis::Symbol,
)
    profile =
        axis == :column ? calc_center_profile(dens2d) :
        axis == :row ? calc_center_row_profile(dens2d) :
        throw(ArgumentError("axis must be :column or :row, got $axis."))
    idx_center = cld(length(x_dens), 2)
    x_half = abs.(x_dens[idx_center:end])
    profile_half = (profile[idx_center:end] .+ reverse(profile[1:idx_center])) ./ 2
    return x_half, profile_half
end

function fit_gaussian_tail_profile(r, params)
    A_gauss, σ_gauss = params
    return @. A_gauss * exp(-r^2 / (2 * σ_gauss^2))
end

function fit_symmetric_profile_tails(
    x_dens::AbstractVector{<:Real},
    ntfr2d::AbstractArray{<:Real,4},
    r_tail_min::Real,
    axis::Symbol,
)
    r_profile, _ = calc_symmetric_profile(x_dens, @view(ntfr2d[:, :, 1, 1]); axis)
    mask_tail = r_profile .> r_tail_min
    any(mask_tail) || throw(ArgumentError("No profile coordinates found above r_tail_min=$r_tail_min."))

    fit_tail = Array{NamedTuple}(undef, size(ntfr2d, 4), size(ntfr2d, 3))
    profile_symm = Array{Vector{Float64}}(undef, size(ntfr2d, 4), size(ntfr2d, 3))
    profile_tailess = similar(profile_symm)
    profile_tail = similar(profile_symm)
    max_r = maximum(r_profile)
    step_r = minimum(diff(r_profile))

    for idx_IB in axes(ntfr2d, 4), idx_istp in axes(ntfr2d, 3)
        r, y = calc_symmetric_profile(x_dens, @view(ntfr2d[:, :, idx_istp, idx_IB]); axis)
        profile_symm[idx_IB, idx_istp] = Float64.(y)
        y_tail = Float64.(y[mask_tail])
        p_init = [max(maximum(y_tail), eps(Float64)), max_r / 2]
        p_lower = [0.0, step_r]
        p_upper = [Inf, max_r * 2]
        fit = curve_fit(fit_gaussian_tail_profile, r[mask_tail], y_tail, p_init; lower=p_lower, upper=p_upper, maxIter=20_000)
        params = coef(fit)
        tail = fit_gaussian_tail_profile(r, params)
        rss_rel = norm(residuals(fit)) / max(norm(y_tail), eps(Float64))
        fit_tail[idx_IB, idx_istp] = (; params, rss_rel)
        profile_tail[idx_IB, idx_istp] = tail
        profile_tailess[idx_IB, idx_istp] = Float64.(y) .- tail
    end

    return (; axis, r_profile, r_tail_min, profile_symm, profile_tail, profile_tailess, fit_tail)
end

function draw_density_row!(
    fig::Figure,
    row::Integer,
    x_dens::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    ntfr2d::AbstractArray{<:Real,4},
    idx_IB::Integer,
    IB::Real;
    colorrange,
    ylims_profile,
    profile_tail_fits,
    is_bottom_row::Bool=false,
)
    Label(fig[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:right)

    axs_dens = Vector{Axis}(undef, length(val_istp))
    axs_profile = Array{Axis}(undef, length(profile_tail_fits), length(val_istp))

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
        pos_center = x_dens[cld(length(x_dens), 2)]
        vlines!(ax, pos_center; color=(:black, 0.14), linewidth=0.6)
        hlines!(ax, pos_center; color=(:black, 0.14), linewidth=0.6)
        hidedecorations!(ax; label=is_bottom_row || idx_istp == 1 ? false : true, ticklabels=!is_bottom_row, ticks=!is_bottom_row, grid=false)

        for (idx_fit, profile_tail_fit) in enumerate(profile_tail_fits)
            idx_col = length(val_istp) + (idx_fit - 1) * length(val_istp) + idx_istp
            ax_profile = Axis(
                fig[row, idx_col];
                xlabel=is_bottom_row ? label_x_dens : "",
                ylabel=idx_istp == 1 ? "symm. $(profile_tail_fit.axis)" : "",
                yaxisposition=idx_istp == 1 ? :left : :right,
                xticks=idx_istp == 1 ? (40:-20:20) : (20:20:40),
            )
            axs_profile[idx_fit, idx_istp] = ax_profile
            r = profile_tail_fit.r_profile
            profile = profile_tail_fit.profile_symm[idx_IB, idx_istp]
            tail = profile_tail_fit.profile_tail[idx_IB, idx_istp]
            tailess = profile_tail_fit.profile_tailess[idx_IB, idx_istp]
            clr_strong = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.95)
            clr_faint = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.32)
            band!(ax_profile, r, zero.(tail), tail; color=(:gray45, 0.25))
            lines!(ax_profile, r, profile; color=clr_faint, linewidth=1.1)
            lines!(ax_profile, r, tailess; color=clr_strong, linewidth=1.8)
            lines!(ax_profile, r, tail; color=(:gray20, 0.55), linewidth=1.0)
            vlines!(ax_profile, profile_tail_fit.r_tail_min; color=(:gray20, 0.35), linewidth=0.8)
            idx_istp == 1 ? xlims!(ax_profile, maximum(r), 0) : xlims!(ax_profile, 0, maximum(r))
            ylims!(ax_profile, ylims_profile)
            !is_bottom_row && hidexdecorations!(ax_profile; grid=false)
        end
    end

    return axs_dens, axs_profile
end

x_dens, val_IB, ntfr2d_mean = h5open(path_data, "r") do file
    x_dens = read(file["x_dens"])
    val_IB = read(file["val_IB"])
    ntfr2d_mean = orient_ntfr2d_axes(read(file["ntfr2d_mean"]), x_dens, val_IB, val_istp)
    return x_dens, val_IB, ntfr2d_mean
end

colorrange_ntfr = (0.0, maximum(ntfr2d_mean))
max_profile = maximum([
    maximum(calc_center_profile(@view ntfr2d_mean[:, :, idx_istp, idx_IB]))
    for idx_istp in axes(ntfr2d_mean, 3), idx_IB in axes(ntfr2d_mean, 4)
])
ylims_profile = (0.0, max_profile * 1.05)
profile_tail_fit_column = fit_symmetric_profile_tails(x_dens, ntfr2d_mean, r_tail_min_profile, :column)
profile_tail_fit_row = fit_symmetric_profile_tails(x_dens, ntfr2d_mean, r_tail_min_profile, :row)
profile_tail_fits = (profile_tail_fit_column, profile_tail_fit_row)
min_tailess_profile = minimum(
    minimum(profile)
    for fit in profile_tail_fits
    for profile in vec(fit.profile_tailess)
)
max_original_profile = maximum(
    maximum(profile)
    for fit in profile_tail_fits
    for profile in vec(fit.profile_symm)
)
ylims_profile_symm = (
    min(0.0, min_tailess_profile * 1.05),
    max_original_profile * 1.05,
)

fig_ntfr = Figure(size=(1580, 2200), fontsize=14)
Label(
    fig_ntfr[0, 1:6];
    text=@sprintf("%s 2D NTFR mean densities, common max %.3g", tag, colorrange_ntfr[2]),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_ntfr[1, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
for (idx_fit, profile_tail_fit) in enumerate(profile_tail_fits)
    for (idx_istp, istp) in enumerate(val_istp)
        idx_col = 2 + (idx_fit - 1) * length(val_istp) + idx_istp
        Label(
            fig_ntfr[1, idx_col];
            text="$(profile_tail_fit.axis) tailess $istp",
            tellwidth=false,
            tellheight=true,
            halign=:center,
            font=:bold,
        )
    end
end

for (idx_IB, IB) in enumerate(val_IB)
    row = idx_IB + 1
    draw_density_row!(
        fig_ntfr,
        row,
        x_dens,
        val_istp,
        ntfr2d_mean,
        idx_IB,
        IB;
        colorrange=colorrange_ntfr,
        ylims_profile=ylims_profile_symm,
        profile_tail_fits,
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
colgap!(fig_ntfr.layout, 1, 8)
colgap!(fig_ntfr.layout, 2, 16)
colgap!(fig_ntfr.layout, 3, 0)
colgap!(fig_ntfr.layout, 4, 14)
colgap!(fig_ntfr.layout, 5, 0)
rowgap!(fig_ntfr.layout, 1, 4)

path_plot_ntfr = joinpath(path_output, "$(tag)_ntfr2d_table.png")
save(path_plot_ntfr, fig_ntfr; backend=CairoMakie)
println("saved $path_plot_ntfr")
