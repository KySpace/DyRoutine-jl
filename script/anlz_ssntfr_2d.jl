using CairoMakie
using FFTW
using HDF5
using ImageFiltering
using LinearAlgebra
using LsqFit
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
title_anlz = "20.Ntfr2D.Abrr.Rotated.SkewedYWithTail"
path_data = joinpath(path_root, "0204_interference", "result", "prfl.h5")
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

tag = "SSNTFR"
log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))
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
xlims_prfl_reconstruct = (0.0, 0.5)
ylims_prfl_reconstruct = (-0.02, 1.6)

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

function calc_rotated_ellipse(
    x0::Real,
    y0::Real,
    sigma_x::Real,
    sigma_y::Real,
    θ::Real;
    n::Integer=160,
    scale::Real=1.0,
)
    ϕ = range(0, 2π; length=n + 1)
    cosθ = cos(θ)
    sinθ = sin(θ)
    x_axis = @. scale * sigma_x * cos(ϕ)
    y_axis = @. scale * sigma_y * sin(ϕ)
    x = @. x0 + cosθ * x_axis - sinθ * y_axis
    y = @. y0 + sinθ * x_axis + cosθ * y_axis
    return x, y
end

function calc_rotated_y_axis(
    x0::Real,
    y0::Real,
    sigma_y::Real,
    θ::Real;
    scale_inner::Real=1.08,
    scale_outer::Real=1.55,
)
    cosθ = cos(θ)
    sinθ = sin(θ)
    dx = -sinθ
    dy = cosθ
    r_inner = scale_inner * sigma_y
    r_outer = scale_outer * sigma_y
    x = [
        x0 + dx * r_inner,
        x0 + dx * r_outer,
        NaN,
        x0 - dx * r_inner,
        x0 - dx * r_outer,
    ]
    y = [
        y0 + dy * r_inner,
        y0 + dy * r_outer,
        NaN,
        y0 - dy * r_inner,
        y0 - dy * r_outer,
    ]
    return x, y
end

function draw_density_row!(
    fig::Figure,
    row::Integer,
    x_dens::AbstractVector{<:Real},
    x_modl::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    ntfr2d_fmt::AbstractMatrix{<:AbstractMatrix},
    ntfr2d_diff_fmt::AbstractMatrix{<:AbstractMatrix},
    idx_IB::Integer,
    IB::Real;
    colorrange,
    colorrange_diff_fmt,
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
    axs_diff = Vector{Axis}(undef, length(val_istp))
    axs_profile = Array{Axis}(undef, length(profile_fits), length(val_istp))
    axs_diag = Vector{Axis}(undef, length(val_istp))
    axs_prfl = Vector{Axis}(undef, length(val_istp))
    clrmap_diff = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360)

    for (idx_istp, istp) in enumerate(val_istp)
        clr_strong = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.95)
        clr_faint = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.32)
        clr_center = Oklch(0.62, 0.16, 145)
        ax = Axis(
            fig[row, idx_istp];
            xlabel=is_bottom_row ? label_x_dens : "",
            ylabel=idx_istp == 1 ? label_x_dens : "",
            aspect=DataAspect(),
        )
        axs_dens[idx_istp] = ax
        dens2d = ntfr2d_fmt[idx_IB, idx_istp]
        clrmap = gen_clrmap_solo(hue_theme_istp[istp])
        heatmap!(ax, x_dens, x_dens, dens2d; colormap=clrmap, colorrange, rasterize=true)
        params_density = fit_density[idx_IB, idx_istp].params
        x0, y0 = params_density[1:2]
        x_ellipse, y_ellipse = calc_rotated_ellipse(
            x0,
            y0,
            params_density[4],
            params_density[5],
            params_density[11],
        )
        x_axis_y, y_axis_y = calc_rotated_y_axis(x0, y0, params_density[5], params_density[11])
        vlines!(ax, x0; color=(:black, 0.16), linewidth=0.7)
        hlines!(ax, y0; color=(:black, 0.16), linewidth=0.7)
        lines!(ax, x_ellipse, y_ellipse; color=(:white, 0.95), linewidth=0.7)
        lines!(ax, x_axis_y, y_axis_y; color=(:white, 0.95), linewidth=0.7)
        hidexdecorations!(ax; label=is_bottom_row ? false : true, ticklabels=is_bottom_row ? false : true, ticks=is_bottom_row ? false : true, grid=false)
        hideydecorations!(ax; label=idx_istp == 1 ? false : true, ticklabels=false, ticks=false, grid=false)
        text!(
            ax,
            0.05, 0.95;
            text=@sprintf("x_0=%.3f\ny_0=%.2f", x0, y0),
            space=:relative,
            color=(clr_center, 0.9),
            fontsize=8,
            align=(:left, :top),
        )

        idx_col_diff = length(val_istp) + idx_istp
        ax_diff = Axis(
            fig[row, idx_col_diff];
            xlabel=is_bottom_row ? label_x_dens : "",
            ylabel=idx_istp == 1 ? label_x_dens : "",
            aspect=DataAspect(),
        )
        axs_diff[idx_istp] = ax_diff
        heatmap!(
            ax_diff,
            x_dens,
            x_dens,
            ntfr2d_diff_fmt[idx_IB, idx_istp];
            colormap=clrmap_diff,
            colorrange=colorrange_diff_fmt[idx_IB, idx_istp],
            rasterize=true,
        )
        vlines!(ax_diff, x0; color=(:black, 0.16), linewidth=0.7)
        hlines!(ax_diff, y0; color=(:black, 0.16), linewidth=0.7)
        hidexdecorations!(ax_diff; label=is_bottom_row ? false : true, ticklabels=is_bottom_row ? false : true, ticks=is_bottom_row ? false : true, grid=false)
        hideydecorations!(ax_diff; label=false, ticklabels=false, ticks=false, grid=false)

        for (idx_fit, profile_fit) in enumerate(profile_fits)
            idx_col = 2 * length(val_istp) + (idx_fit - 1) * length(val_istp) + idx_istp
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
            params_fit = profile_data.fit_density.params
            text_fit =
                profile_data.axis == :column ?
                @sprintf("β=%.3f\nσ_y=%.2f α_y=%.2f\nθ=%.2f", profile_data.beta, params_fit[5], params_fit[10], params_fit[11]) :
                @sprintf("β=%.3f\nσ_x=%.2f α_x=%.2f\nθ=%.2f", profile_data.beta, params_fit[4], params_fit[9], params_fit[11])
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
                text=text_fit,
                color=(clr_center, 0.9),
                fontsize=8,
                align=(:left, :top),
            )
            !is_bottom_row && hidexdecorations!(ax_profile; grid=false)
        end

        idx_col_diag = 2 * length(val_istp) + length(profile_fits) * length(val_istp) + idx_istp
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

        idx_col_prfl = 2 * length(val_istp) + length(profile_fits) * length(val_istp) + length(val_istp) + idx_istp
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
        lines!(ax_prfl, x_modl, prfl_modl_fit[idx_IB, idx_istp]; color=clr_fit, linewidth=1.8)
        lines!(ax_prfl, x_modl, @view(prfl_inco[:, idx_istp, idx_IB]); color=clr_theme_faint, linewidth=1.3, linestyle=:dash)
        lines!(ax_prfl, x_modl, @view(prfl_cohr[:, idx_istp, idx_IB]); color=clr_theme, linewidth=1.5)
        xlims!(ax_prfl, xlims_prfl)
        ylims!(ax_prfl, ylims_prfl)
        !is_bottom_row && hidexdecorations!(ax_prfl; grid=false)
    end

    return axs_dens, axs_diff, axs_profile, axs_diag, axs_prfl
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

t_stage = log_step("preparing density/profile ranges")
colorrange_ntfr = (0.0, maximum(maximum, ntfr2d_mean))
max_profile = maximum([
    maximum(calc_grid_center_profile(dens2d))
    for dens2d in ntfr2d_mean
])
ylims_profile = (0.0, max_profile * 1.05)
log_done("prepared density/profile ranges", t_stage)

t_stage = log_step("fitting centered 2D NTFR densities")
fit_centered = fit_centered_density_profiles(
    x_dens,
    ntfr2d_mean,
    r_tail_min_profile;
    center_bound=fit_center_bound,
    stride=fit_stride_2d,
    threshold=fit_threshold_log_2d,
    sigma_wide_min=fit_sigma_wide_min,
    maxiter=fit_maxiter_2d,
    model_center,
    log_tag=tag,
    val_IB,
    val_istp,
)
log_done("fit centered 2D NTFR densities", t_stage)
fit_density = fit_centered.fit_density
profile_fits = fit_centered.profile_fits

t_stage = log_step("reconstructing fitted densities and modulation profiles")
ntfr2d_fit = calc_reconstructed_ntfr2d(x_dens, fit_density)
ntfr2d_diff = map(-, ntfr2d_mean, ntfr2d_fit)
max_abs_diff = maximum(maximum(abs, diff2d) for diff2d in ntfr2d_diff)
colorrange_diff = map(ntfr2d_diff) do diff2d
    c = maximum(abs, diff2d)
    (-c, c)
end
prfl_modl_fit = calc_reconstructed_prfl_modl(ntfr2d_fit, smwh_reconstruct; step_modl)
length(x_modl) == length(prfl_modl_fit[1]) || throw(DimensionMismatch(
    "x_modl length $(length(x_modl)) must match reconstructed profile length $(length(prfl_modl_fit[1])).",
))
log_done("reconstructed fitted densities and modulation profiles", t_stage)

t_stage = log_step("building plot ranges")
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
log_done("built plot ranges", t_stage)

t_stage = log_step("building figure axes and labels")
fig_ntfr = Figure(fontsize=14)
Label(
    fig_ntfr[0, 1:12];
    text=@sprintf(
        "%s 2D NTFR mean densities, fit residuals with per-panel scale, cocenter Gaussian tail + rotated narrow peak |> (_ + βN²) fit, mask > %.1g, σ_wide ≥ %.0f μm, common max %.3g, max residual ±%.3g, Δk=%.5f",
        tag,
        fit_threshold_log_2d,
        fit_sigma_wide_min,
        colorrange_ntfr[2],
        max_abs_diff,
        step_modl,
    ),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_ntfr[1, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
for (idx_istp, istp) in enumerate(val_istp)
    idx_col = length(val_istp) + idx_istp
    Label(fig_ntfr[1, idx_col]; text="diff $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
for (idx_fit, profile_fit) in enumerate(profile_fits)
    axis_name = profile_fit[1, 1].axis
    for (idx_istp, istp) in enumerate(val_istp)
        idx_col = 2 * length(val_istp) + (idx_fit - 1) * length(val_istp) + idx_istp
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
    idx_col = 2 * length(val_istp) + length(profile_fits) * length(val_istp) + idx_istp
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
    idx_col = 2 * length(val_istp) + length(profile_fits) * length(val_istp) + length(val_istp) + idx_istp
    Label(
        fig_ntfr[1, idx_col];
        text="FT prfl $istp",
        tellwidth=false,
        tellheight=true,
        halign=:center,
        font=:bold,
    )
end
log_done("built figure axes and labels", t_stage)

t_stage = log_step("plotting density rows")
for (idx_IB, IB) in enumerate(val_IB)
    println("  [$tag] plotting IB_idx=$idx_IB IB=$IB")
    flush(stdout)
    row = idx_IB + 1
    draw_density_row!(
        fig_ntfr,
        row,
        x_dens,
        x_modl,
        val_istp,
        ntfr2d_mean,
        ntfr2d_diff,
        idx_IB,
        IB;
        colorrange=colorrange_ntfr,
        colorrange_diff_fmt=colorrange_diff,
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
log_done("plotted density rows", t_stage)

colsize!(fig_ntfr.layout, 1, Fixed(105))
colsize!(fig_ntfr.layout, 2, Fixed(105))
colsize!(fig_ntfr.layout, 3, Fixed(105))
colsize!(fig_ntfr.layout, 4, Fixed(105))
colsize!(fig_ntfr.layout, 5, Fixed(170))
colsize!(fig_ntfr.layout, 6, Fixed(170))
colsize!(fig_ntfr.layout, 7, Fixed(170))
colsize!(fig_ntfr.layout, 8, Fixed(170))
colsize!(fig_ntfr.layout, 9, Fixed(170))
colsize!(fig_ntfr.layout, 10, Fixed(170))
colsize!(fig_ntfr.layout, 11, Fixed(170))
colsize!(fig_ntfr.layout, 12, Fixed(170))
colgap!(fig_ntfr.layout, 1, 8)
colgap!(fig_ntfr.layout, 2, 16)
colgap!(fig_ntfr.layout, 3, 0)
colgap!(fig_ntfr.layout, 4, 14)
colgap!(fig_ntfr.layout, 5, 0)
colgap!(fig_ntfr.layout, 6, 14)
colgap!(fig_ntfr.layout, 7, 0)
colgap!(fig_ntfr.layout, 8, 14)
colgap!(fig_ntfr.layout, 9, 0)
colgap!(fig_ntfr.layout, 10, 14)
colgap!(fig_ntfr.layout, 11, 0)
rowgap!(fig_ntfr.layout, 1, 4)

resize_to_layout!(fig_ntfr)
filename_plot_ntfr = "$(tag)_ntfr2d_table"
t_stage = log_step("saving figure outputs")
for ext in ("png", "pdf")
    save(joinpath(path_output, "$filename_plot_ntfr.$ext"), fig_ntfr; backend=CairoMakie)
end
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
cp(joinpath(@__DIR__, "..", "src", "modlntfr.jl"), joinpath(path_output, "modlntfr.jl"); force=true)
log_done("saved figure outputs", t_stage)
println("saved $(joinpath(path_output, "$filename_plot_ntfr.png"))")
