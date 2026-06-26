using CairoMakie
using HDF5
using LinearAlgebra
using LsqFit
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "fitmodels.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
title_anlz = "02.PrflStacked"
path_data = joinpath(path_root, "0204_interference", "result", "prfl.h5")
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

tag = "SSNTFR"

val_istp = ["162", "164"]
label_x_modl = "wavenum (μm⁻¹)"
range_x_plot = (0.0, 1.2)
range_x_colorrange = (0.1, 0.6)
range_x_fit_inco = (0.0, 0.6)
range_x_fit_cohr_tail = (0.1, 0.6)
clr_IB_endpoints = (
    low=(l=0.34, c=0.10, h=255.0),
    high=(l=0.72, c=0.18, h=25.0),
)

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

function calc_prfl_colorrange(
    prfl::AbstractArray{<:Real,3},
    x_modl::AbstractVector{<:Real},
    range_x::Tuple{<:Real,<:Real},
)
    xmin, xmax = range_x
    mask_x = (x_modl .> xmin) .& (x_modl .< xmax)
    any(mask_x) || throw(ArgumentError("No x_modl values found in colorrange window $range_x."))
    val_max = maximum(@view prfl[mask_x, :, :])
    val_max > 0 || throw(ArgumentError("Nonpositive profile maximum $val_max in colorrange window $range_x."))
    return (0.0, val_max)
end

function calc_log_ylims(colorrange::Tuple{<:Real,<:Real})
    _, ymax = colorrange
    ymax > 0 || throw(ArgumentError("Log-scale y upper limit must be positive, got $ymax."))
    return (ymax * 1e-4, ymax)
end

function select_x_range(x_modl::AbstractVector{<:Real}, range_x::Tuple{<:Real,<:Real})
    xmin, xmax = range_x
    mask = (x_modl .>= xmin) .& (x_modl .<= xmax)
    any(mask) || throw(ArgumentError("No x_modl values found in range $range_x."))
    return mask
end

function gen_IB_ticks(val_IB::AbstractVector{<:Real}; n_ticks::Integer=9)
    return round.(range(minimum(val_IB), maximum(val_IB); length=n_ticks); digits=3)
end

function calc_IB_color(val_IB::Real, val_IB_min::Real, val_IB_max::Real, clr_endpoints::NamedTuple)
    t = val_IB_max == val_IB_min ? 0.0 : (val_IB - val_IB_min) / (val_IB_max - val_IB_min)
    l = (1 - t) * clr_endpoints.low.l + t * clr_endpoints.high.l
    c = (1 - t) * clr_endpoints.low.c + t * clr_endpoints.high.c
    h = (1 - t) * clr_endpoints.low.h + t * clr_endpoints.high.h
    return RGBAf(Oklch(l, c, h), 0.88)
end

function gen_clrmap_IB(clr_endpoints::NamedTuple; n::Integer=256, alpha=1.0)
    return [
        begin
            l = (1 - t) * clr_endpoints.low.l + t * clr_endpoints.high.l
            c = (1 - t) * clr_endpoints.low.c + t * clr_endpoints.high.c
            h = (1 - t) * clr_endpoints.low.h + t * clr_endpoints.high.h
            RGBAf(Oklch(l, c, h), alpha)
        end
        for t in range(0, 1; length=n)
    ]
end

function calc_inco_tail_fits(
    x_modl::AbstractVector{<:Real},
    prfl_inco::AbstractArray{<:Real,3};
    range_x_fit::Tuple{<:Real,<:Real},
)
    mask_fit = select_x_range(x_modl, range_x_fit)
    k_fit = Float64.(x_modl[mask_fit])
    n_IB = size(prfl_inco, 3)
    n_istp = size(prfl_inco, 2)
    fit_inco = Array{NamedTuple}(undef, n_IB, n_istp)
    prfl_inco_tailess = similar(prfl_inco, Float64)
    comp_center = similar(prfl_inco, Float64)
    comp_tail = similar(prfl_inco, Float64)
    comp_side = similar(prfl_inco, Float64)

    for idx_IB in 1:n_IB, idx_istp in 1:n_istp
        y_fit = Float64.(@view prfl_inco[mask_fit, idx_istp, idx_IB])
        amp_main = max(maximum(y_fit), eps(Float64))
        idx_side_hint = findmax(y_fit .* (k_fit .> 0.08))[2]
        p_hint = clamp(k_fit[idx_side_hint], 0.12, 0.30)
        p_init = [amp_main, 0.045, max(y_fit[idx_side_hint] / 2, 1e-4), 0.05, p_hint, max(y_fit[end], 1e-4), 0.18]
        p_lower = [0.0, 0.005, 0.0, 0.010, 0.06, 0.0, 0.02]
        p_upper = [Inf, 0.20, Inf, 0.180, 0.45, Inf, 1.00]
        fit = curve_fit(fit_prfl_modl_twinpeak_decay_1d_model, k_fit, y_fit, p_init; lower=p_lower, upper=p_upper, maxIter=20_000)
        params = coef(fit)
        rss_rel = norm(residuals(fit)) / max(norm(y_fit), eps(Float64))

        center = @. params[1] * exp(-x_modl^2 / (2 * params[2]^2))
        side = @. params[3] * exp(-(x_modl - params[5])^2 / (2 * params[4]^2))
        tail = fit_prfl_modl_twinpeak_decay_1d_tail(x_modl, params)
        comp_center[:, idx_istp, idx_IB] .= center
        comp_side[:, idx_istp, idx_IB] .= side
        comp_tail[:, idx_istp, idx_IB] .= tail
        prfl_inco_tailess[:, idx_istp, idx_IB] .= @view(prfl_inco[:, idx_istp, idx_IB]) .- center .- tail
        fit_inco[idx_IB, idx_istp] = (; params, rss_rel)
    end

    return (; fit_inco, prfl_inco_tailess, comp_center, comp_side, comp_tail)
end

function fit_log_prfl_modl_twinpeak_decay_1d_model(k, params)
    return log.(max.(fit_prfl_modl_twinpeak_decay_1d_model(k, params), eps(Float64)))
end

function fit_log_prfl_modl_sidepeak_decay_1d_model(k, params)
    return log.(max.(fit_prfl_modl_sidepeak_decay_1d_model(k, params), eps(Float64)))
end

function fit_common_cohr_tail(
    x_modl::AbstractVector{<:Real},
    prfl_cohr::AbstractArray{<:Real,3};
    range_x_fit::Tuple{<:Real,<:Real},
)
    mask_fit = select_x_range(x_modl, range_x_fit)
    k_fit = Float64.(x_modl[mask_fit])
    prfl_total_avg = vec(mean(prfl_cohr; dims=(2, 3)))
    y_fit = log.(max.(Float64.(prfl_total_avg[mask_fit]), eps(Float64)))
    mask_side_hint = (k_fit .> 0.16) .& (k_fit .< 0.34)
    idx_side_hint = findmax(ifelse.(mask_side_hint, y_fit, -Inf))[2]
    p_hint = clamp(k_fit[idx_side_hint], 0.20, 0.30)
    y_fit_linear = Float64.(prfl_total_avg[mask_fit])
    p_init = [max(y_fit_linear[idx_side_hint], 1e-4), 0.05, p_hint, max(y_fit_linear[end], 1e-5), 0.18]
    p_lower = [0.0, 0.010, 0.16, 0.0, 0.02]
    p_upper = [Inf, 0.180, 0.36, Inf, 1.00]
    fit = curve_fit(fit_log_prfl_modl_sidepeak_decay_1d_model, k_fit, y_fit, p_init; lower=p_lower, upper=p_upper, maxIter=20_000)
    params = coef(fit)
    tail = fit_prfl_modl_sidepeak_decay_1d_tail(x_modl, params)
    side = @. params[1] * exp(-(x_modl - params[3])^2 / (2 * params[2]^2))
    prfl_cohr_tailess = similar(prfl_cohr, Float64)
    for idx_IB in axes(prfl_cohr, 3), idx_istp in axes(prfl_cohr, 2)
        prfl_cohr_tailess[:, idx_istp, idx_IB] .= @view(prfl_cohr[:, idx_istp, idx_IB]) .- tail
    end
    rss_log_rel = norm(residuals(fit)) / max(norm(y_fit), eps(Float64))
    return (; params, rss_log_rel, tail, side, prfl_total_avg, prfl_cohr_tailess)
end

function build_profile_variant_arrays(prfl_inco, prfl_cohr, prfl_inco_tailess, prfl_cohr_tailess)
    val_variant = [:inco, :cohr]
    prfl_original_fmt = Array{Vector{Float64}}(undef, size(prfl_inco, 3), size(prfl_inco, 2), length(val_variant))
    prfl_tailess_fmt = similar(prfl_original_fmt)
    for idx_IB in axes(prfl_inco, 3), idx_istp in axes(prfl_inco, 2)
        prfl_original_fmt[idx_IB, idx_istp, 1] = vec(Float64.(@view prfl_inco[:, idx_istp, idx_IB]))
        prfl_original_fmt[idx_IB, idx_istp, 2] = vec(Float64.(@view prfl_cohr[:, idx_istp, idx_IB]))
        prfl_tailess_fmt[idx_IB, idx_istp, 1] = vec(Float64.(@view prfl_inco_tailess[:, idx_istp, idx_IB]))
        prfl_tailess_fmt[idx_IB, idx_istp, 2] = vec(Float64.(@view prfl_cohr_tailess[:, idx_istp, idx_IB]))
    end
    return (; val_variant, prfl_original_fmt, prfl_tailess_fmt)
end

function draw_stacked_profile_heatmaps!(
    fig::Figure,
    row::Integer,
    prfl::AbstractArray{<:Real,3},
    x_modl::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    label_prfl::AbstractString;
    range_x_plot::Tuple{<:Real,<:Real},
    colorrange::Tuple{<:Real,<:Real},
)
    axs = Vector{Axis}(undef, length(val_istp))
    hm = nothing
    ticks_IB = gen_IB_ticks(val_IB)

    Label(fig[row - 2, 1:length(val_istp)]; text=label_prfl, tellwidth=false, tellheight=true, halign=:center, font=:bold)
    for (idx_istp, istp) in enumerate(val_istp)
        ax = Axis(
            fig[row, idx_istp];
            xlabel=label_x_modl,
            ylabel="IB (A)",
            yaxisposition=idx_istp == 1 ? :left : :right,
            xticks=idx_istp == 1 ? (1.2:-0.2:0.2) : (0.2:0.2:1.2),
            yticks=ticks_IB,
        )
        axs[idx_istp] = ax

        clrmap = gen_clrmap_solo(hue_theme_istp[istp])
        hm = heatmap!(
            ax,
            x_modl,
            val_IB,
            @view(prfl[:, idx_istp, :]);
            colormap=clrmap,
            colorrange,
            rasterize=true,
        )
        idx_istp == 1 ? xlims!(ax, reverse(range_x_plot)) : xlims!(ax, range_x_plot)
        Label(fig[row - 1, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center)
    end

    Colorbar(fig[row, length(val_istp) + 1], hm; label="profile")
    return axs
end

function draw_overlaid_profile_lines!(
    fig::Figure,
    row::Integer,
    prfl::AbstractArray{<:Real,3},
    x_modl::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    label_prfl::AbstractString;
    range_x_plot::Tuple{<:Real,<:Real},
    colorrange::Tuple{<:Real,<:Real},
    clr_endpoints::NamedTuple,
)
    axs = Vector{Axis}(undef, length(val_istp))
    val_IB_min, val_IB_max = extrema(val_IB)

    Label(fig[row - 2, 1:length(val_istp)]; text=label_prfl, tellwidth=false, tellheight=true, halign=:center, font=:bold)
    for (idx_istp, istp) in enumerate(val_istp)
        ax = Axis(
            fig[row, idx_istp];
            xlabel=label_x_modl,
            ylabel="profile",
            yaxisposition=idx_istp == 1 ? :left : :right,
            yscale=log10,
            xticks=idx_istp == 1 ? (1.2:-0.2:0.2) : (0.2:0.2:1.2),
        )
        axs[idx_istp] = ax

        for (idx_IB, IB) in enumerate(val_IB)
            clr = calc_IB_color(IB, val_IB_min, val_IB_max, clr_endpoints)
            lines!(ax, x_modl, @view(prfl[:, idx_istp, idx_IB]); color=clr, linewidth=1.8)
        end

        idx_istp == 1 ? xlims!(ax, reverse(range_x_plot)) : xlims!(ax, range_x_plot)
        ylims!(ax, calc_log_ylims(colorrange))
        Label(fig[row - 1, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center)
    end

    Colorbar(
        fig[row, length(val_istp) + 1];
        colormap=gen_clrmap_IB(clr_endpoints; alpha=0.8),
        limits=extrema(val_IB),
        label="IB (A)",
    )
    return axs
end

function draw_cohr_average_lines!(
    fig::Figure,
    row::Integer,
    prfl_cohr::AbstractArray{<:Real,3},
    tail_cohr::NamedTuple,
    x_modl::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString};
    range_x_plot::Tuple{<:Real,<:Real},
    ylims,
)
    axs = Vector{Axis}(undef, length(val_istp))
    prfl_total_avg = tail_cohr.prfl_total_avg

    for (idx_istp, istp) in enumerate(val_istp)
        ax = Axis(
            fig[row, idx_istp];
            xlabel=label_x_modl,
            ylabel="profile",
            yaxisposition=idx_istp == 1 ? :left : :right,
            yscale=log10,
            xticks=idx_istp == 1 ? (1.2:-0.2:0.2) : (0.2:0.2:1.2),
        )
        axs[idx_istp] = ax
        prfl_istp_avg = vec(mean(@view(prfl_cohr[:, idx_istp, :]); dims=2))
        band!(ax, x_modl, tail_cohr.tail, tail_cohr.tail .+ tail_cohr.side; color=(:mediumseagreen, 0.30))
        lines!(ax, x_modl, tail_cohr.tail; color=(:seagreen4, 0.75), linewidth=2.2)
        lines!(ax, x_modl, prfl_total_avg; color=(:gray30, 0.45), linewidth=3.0)
        lines!(ax, x_modl, prfl_istp_avg; color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.95), linewidth=2.2)
        idx_istp == 1 ? xlims!(ax, reverse(range_x_plot)) : xlims!(ax, range_x_plot)
        ylims!(ax, ylims)
        Label(fig[row - 1, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center)
    end

    return axs, prfl_total_avg
end

function draw_tail_diagnostic_side!(
    ax::Axis,
    x_modl::AbstractVector{<:Real},
    y_original::AbstractVector{<:Real},
    y_tailess::AbstractVector{<:Real},
    idx_istp::Integer;
    range_x_plot::Tuple{<:Real,<:Real},
    ylims,
    color_original,
    color_tailess,
)
    lines!(ax, x_modl, y_original; color=color_original, linewidth=1.6)
    lines!(ax, x_modl, y_tailess; color=color_tailess, linewidth=1.2)
    idx_istp == 1 ? xlims!(ax, reverse(range_x_plot)) : xlims!(ax, range_x_plot)
    ylims!(ax, ylims)
    return ax
end

function draw_inco_tail_diagnostic_side!(
    ax::Axis,
    x_modl::AbstractVector{<:Real},
    y_original::AbstractVector{<:Real},
    center::AbstractVector{<:Real},
    side::AbstractVector{<:Real},
    tail::AbstractVector{<:Real},
    idx_istp::Integer;
    range_x_plot::Tuple{<:Real,<:Real},
    ylims,
    color_original,
)
    base = center .+ tail
    side_only = y_original .- base
    band!(ax, x_modl, zero.(tail), tail; color=(:gray55, 0.22))
    band!(ax, x_modl, tail, base; color=(:gray25, 0.25))
    band!(ax, x_modl, base, base .+ side; color=(:mediumseagreen, 0.28))
    lines!(ax, x_modl, y_original; color=color_original, linewidth=1.6)
    lines!(ax, x_modl, side_only; color=(:seagreen4, 0.95), linewidth=1.0)
    idx_istp == 1 ? xlims!(ax, reverse(range_x_plot)) : xlims!(ax, range_x_plot)
    ylims!(ax, ylims)
    return ax
end

function draw_tail_diagnostic_duet!(
    fig::Figure,
    row::Integer,
    col_start::Integer,
    x_modl::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    label_kind::AbstractString,
    idx_IB::Integer,
    prfl_original::AbstractArray{<:Real,3},
    prfl_tailess::AbstractArray{<:Real,3};
    comp_center=nothing,
    comp_side=nothing,
    comp_tail=nothing,
    range_x_plot::Tuple{<:Real,<:Real},
    ylims,
    is_bottom_row::Bool=false,
)
    axs = Vector{Axis}(undef, length(val_istp))
    for (idx_istp, istp) in enumerate(val_istp)
        ax = Axis(
            fig[row, col_start + idx_istp - 1];
            xlabel=is_bottom_row ? label_x_modl : "",
            ylabel=idx_istp == 1 ? label_kind : "",
            yaxisposition=idx_istp == 1 ? :left : :right,
            xticks=idx_istp == 1 ? (1.2:-0.4:0.4) : (0.4:0.4:1.2),
        )
        axs[idx_istp] = ax
        color_original = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.35)
        color_tailess = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.92)

        if isnothing(comp_tail)
            draw_tail_diagnostic_side!(
                ax,
                x_modl,
                @view(prfl_original[:, idx_istp, idx_IB]),
                @view(prfl_tailess[:, idx_istp, idx_IB]),
                idx_istp;
                range_x_plot,
                ylims,
                color_original,
                color_tailess,
            )
        else
            draw_inco_tail_diagnostic_side!(
                ax,
                x_modl,
                @view(prfl_original[:, idx_istp, idx_IB]),
                @view(comp_center[:, idx_istp, idx_IB]),
                @view(comp_side[:, idx_istp, idx_IB]),
                @view(comp_tail[:, idx_istp, idx_IB]),
                idx_istp;
                range_x_plot,
                ylims,
                color_original,
            )
        end
        !is_bottom_row && hidexdecorations!(ax; grid=false)
    end
    colgap!(fig.layout, col_start, 0)
    return axs
end

x_modl, val_IB, prfl_inco, prfl_cohr = h5open(path_data, "r") do file
    x_modl = read(file["x_modl"])
    val_IB = read(file["val_IB"])
    prfl_inco = orient_prfl_axes(read(file["prfl_inco"]), x_modl, val_IB, val_istp)
    prfl_cohr = orient_prfl_axes(read(file["prfl_cohr"]), x_modl, val_IB, val_istp)
    return x_modl, val_IB, prfl_inco, prfl_cohr
end

colorrange_inco = calc_prfl_colorrange(prfl_inco, x_modl, range_x_colorrange)
colorrange_cohr = calc_prfl_colorrange(prfl_cohr, x_modl, range_x_colorrange)

tail_inco = calc_inco_tail_fits(x_modl, prfl_inco; range_x_fit=range_x_fit_inco)
tail_cohr = fit_common_cohr_tail(x_modl, prfl_cohr; range_x_fit=range_x_fit_cohr_tail)
prfl_tail = build_profile_variant_arrays(
    prfl_inco,
    prfl_cohr,
    tail_inco.prfl_inco_tailess,
    tail_cohr.prfl_cohr_tailess,
)

fig_prfl = Figure(size=(980, 760))
Label(
    fig_prfl[0, 1:3];
    text=@sprintf(
        "%s stacked modulation profiles, colorrange from %.1f < x_modl < %.1f",
        tag,
        range_x_colorrange...
    ),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)

axs_inco = draw_stacked_profile_heatmaps!(
    fig_prfl,
    3,
    prfl_inco,
    x_modl,
    val_IB,
    val_istp,
    "inco";
    range_x_plot,
    colorrange=colorrange_inco,
)
axs_cohr = draw_stacked_profile_heatmaps!(
    fig_prfl,
    6,
    prfl_cohr,
    x_modl,
    val_IB,
    val_istp,
    "cohr";
    range_x_plot,
    colorrange=colorrange_cohr,
)

rowgap!(fig_prfl.layout, 1, 2)
rowgap!(fig_prfl.layout, 2, 4)
rowgap!(fig_prfl.layout, 3, 14)
rowgap!(fig_prfl.layout, 4, 2)
rowgap!(fig_prfl.layout, 5, 4)
colgap!(fig_prfl.layout, 1, 0)
colgap!(fig_prfl.layout, 2, 10)

path_plot_prfl = joinpath(path_output, "$(tag)_prfl_stacked.png")
save(path_plot_prfl, fig_prfl; backend=CairoMakie)
println("saved $path_plot_prfl")

fig_prfl_lines = Figure(size=(980, 760))
Label(
    fig_prfl_lines[0, 1:3];
    text=@sprintf(
        "%s overlaid modulation profiles, IB OKLCH low=(%.2f, %.2f, %.1f), high=(%.2f, %.2f, %.1f)",
        tag,
        clr_IB_endpoints.low.l,
        clr_IB_endpoints.low.c,
        clr_IB_endpoints.low.h,
        clr_IB_endpoints.high.l,
        clr_IB_endpoints.high.c,
        clr_IB_endpoints.high.h,
    ),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)

axs_inco_lines = draw_overlaid_profile_lines!(
    fig_prfl_lines,
    3,
    prfl_inco,
    x_modl,
    val_IB,
    val_istp,
    "inco";
    range_x_plot,
    colorrange=colorrange_inco,
    clr_endpoints=clr_IB_endpoints,
)
axs_cohr_lines = draw_overlaid_profile_lines!(
    fig_prfl_lines,
    6,
    prfl_cohr,
    x_modl,
    val_IB,
    val_istp,
    "cohr";
    range_x_plot,
    colorrange=colorrange_cohr,
    clr_endpoints=clr_IB_endpoints,
)

rowgap!(fig_prfl_lines.layout, 1, 2)
rowgap!(fig_prfl_lines.layout, 2, 4)
rowgap!(fig_prfl_lines.layout, 3, 14)
rowgap!(fig_prfl_lines.layout, 4, 2)
rowgap!(fig_prfl_lines.layout, 5, 4)
colgap!(fig_prfl_lines.layout, 1, 0)
colgap!(fig_prfl_lines.layout, 2, 10)

path_plot_prfl_lines = joinpath(path_output, "$(tag)_prfl_lines.png")
save(path_plot_prfl_lines, fig_prfl_lines; backend=CairoMakie)
println("saved $path_plot_prfl_lines")

mask_x_tail_plot = select_x_range(x_modl, range_x_plot)
ylim_inco_tail = (-0.5, 2.5)
ylim_cohr_tail = (-0.2, 0.8)

fig_tail = Figure(size=(1500, 2250), fontsize=14)
Label(
    fig_tail[0, 1:4];
    text=@sprintf(
        "%s profile tail removal, inco fit range %.1f-%.1f, cohr common tail range %.1f-%.1f",
        tag,
        range_x_fit_inco...,
        range_x_fit_cohr_tail...
    ),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)
Label(fig_tail[1, 1:2]; text="inco", tellwidth=false, tellheight=true, halign=:center, font=:bold)
Label(fig_tail[1, 3:4]; text="cohr", tellwidth=false, tellheight=true, halign=:center, font=:bold)
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_tail[2, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center)
    Label(fig_tail[2, 2 + idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center)
end

for (idx_IB, IB) in enumerate(val_IB)
    row = idx_IB + 2
    Label(fig_tail[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:right)
    draw_tail_diagnostic_duet!(
        fig_tail,
        row,
        1,
        x_modl,
        val_istp,
        "profile",
        idx_IB,
        prfl_inco,
        tail_inco.prfl_inco_tailess;
        comp_center=tail_inco.comp_center,
        comp_side=tail_inco.comp_side,
        comp_tail=tail_inco.comp_tail,
        range_x_plot,
        ylims=ylim_inco_tail,
        is_bottom_row=idx_IB == length(val_IB),
    )
    draw_tail_diagnostic_duet!(
        fig_tail,
        row,
        3,
        x_modl,
        val_istp,
        "profile",
        idx_IB,
        prfl_cohr,
        tail_cohr.prfl_cohr_tailess;
        range_x_plot,
        ylims=ylim_cohr_tail,
        is_bottom_row=idx_IB == length(val_IB),
    )
    rowsize!(fig_tail.layout, row, Fixed(90))
end
colgap!(fig_tail.layout, 1, 0)
colgap!(fig_tail.layout, 2, 18)
colgap!(fig_tail.layout, 3, 0)
rowgap!(fig_tail.layout, 1, 4)
rowgap!(fig_tail.layout, 2, 4)

path_plot_tail = joinpath(path_output, "$(tag)_prfl_tail_diagnostic.png")
save(path_plot_tail, fig_tail; backend=CairoMakie)
println("saved $path_plot_tail")

fig_cohr_avg = Figure(size=(980, 420))
Label(
    fig_cohr_avg[0, 1:2];
    text="$(tag) cohr average profiles",
    tellwidth=false,
    tellheight=true,
    halign=:left,
)
Label(
    fig_cohr_avg[1, 1:2];
    text="per-istp mean over IB, gray: total mean, green: fitted tail + sidepeak",
    tellwidth=false,
    tellheight=true,
    halign=:center,
    font=:bold,
)
axs_cohr_avg, prfl_cohr_total_avg = draw_cohr_average_lines!(
    fig_cohr_avg,
    3,
    prfl_cohr,
    tail_cohr,
    x_modl,
    val_istp;
    range_x_plot,
    ylims=calc_log_ylims(colorrange_cohr),
)
rowgap!(fig_cohr_avg.layout, 1, 4)
rowgap!(fig_cohr_avg.layout, 2, 4)
colgap!(fig_cohr_avg.layout, 1, 0)

path_plot_cohr_avg = joinpath(path_output, "$(tag)_prfl_cohr_average.png")
save(path_plot_cohr_avg, fig_cohr_avg; backend=CairoMakie)
println("saved $path_plot_cohr_avg")
