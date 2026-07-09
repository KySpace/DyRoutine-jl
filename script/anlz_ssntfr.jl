using CairoMakie
using FFTW
using HDF5
using ImageFiltering
using JLD2
using LinearAlgebra
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
use_src_profiles = false
title_anlz = use_src_profiles ? "28.IncoCohrModlNtfr.[WL-migration]" : "30.IncoCohrModlNtfr.[reconstr.29]"
path_data = joinpath(path_root, "0204_interference", "result", "prfl.h5")
path_model_results = joinpath(
    path_root,
    "AnlzRoutine",
    "29.Ntfr2D.Abrr.LinearWeight.Lib",
    "SSNTFR_ntfr2d_model_results.jld2",
)
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

tag = "SSNTFR"

val_istp = ["162", "164"]
label_x_modl = "wavenum (μm⁻¹)"
range_x_plot = (0.0, 1.2)
range_x_colorrange = (0.1, 0.6)
r_tail_min_profile = 20.0
fit_center_bound = 12.0
fit_stride_2d = 3
fit_maxiter_2d = 10_000
fit_threshold_log_2d = 1.5e-1
fit_sigma_wide_min = 15.0
model_center = :gaussian
smwh_reconstruct = (150, 150)
clr_IB_endpoints = (
    low=(l=0.34, c=0.10, h=255.0),
    high=(l=0.72, c=0.18, h=25.0),
)

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

function calc_modl_tail(
    x_modl::AbstractVector{<:Real},
    prfl::AbstractArray{<:Real,3},
    prfl_modl_fit::AbstractArray{<:Real,3},
)
    size(prfl) == size(prfl_modl_fit) || throw(DimensionMismatch(
        "profile size $(size(prfl)) must match prfl_modl_fit size $(size(prfl_modl_fit)).",
    ))
    prfl_tailess = similar(prfl, Float64)
    sel_center = abs.(x_modl) .<= 0.1
    for idx_IB in axes(prfl, 3), idx_istp in axes(prfl, 2)
        ids = (:, idx_istp, idx_IB)
        scale = sum(prfl[ids...][sel_center]) / sum(prfl_modl_fit[ids...][sel_center])
        prfl_tailess[ids...] .=
            @view(prfl[ids...]) .- scale .* @view(prfl_modl_fit[ids...])
    end
    prfl_total_avg = vec(mean(prfl; dims=(2, 3)))
    tail = vec(mean(prfl_modl_fit; dims=(2, 3)))
    tail_istp_avg = dropdims(mean(prfl_modl_fit; dims=3); dims=3)
    return (; prfl_modl_fit, tail, tail_istp_avg, prfl_total_avg, prfl_tailess)
end

function load_prfl_modl_fit_jld2(
    path_results::AbstractString,
    x_modl::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
)
    isfile(path_results) || throw(ArgumentError("Missing reconstructed model results file: $path_results"))
    payload = JLD2.jldopen(path_results, "r") do file
        (;
            x_modl=file["x_modl"],
            val_IB=file["val_IB"],
            val_istp=String.(file["val_istp"]),
            prfl_modl_fit=file["prfl_modl_fit"],
        )
    end

    length(payload.x_modl) == length(x_modl) || throw(DimensionMismatch(
        "Reconstructed model x_modl length $(length(payload.x_modl)) must match current x_modl length $(length(x_modl)).",
    ))
    all(isapprox.(payload.x_modl, x_modl; rtol=0.0, atol=1e-12)) || @warn(DimensionMismatch(
        "Reconstructed model x_modl values do not match current x_modl.",
    ))
    payload.val_IB == val_IB || throw(DimensionMismatch(
        "Reconstructed model val_IB $(payload.val_IB) must match current val_IB $val_IB.",
    ))
    payload.val_istp == String.(val_istp) || throw(DimensionMismatch(
        "Reconstructed model val_istp $(payload.val_istp) must match current val_istp $(String.(val_istp)).",
    ))
    size(payload.prfl_modl_fit) == (length(x_modl), length(val_istp), length(val_IB)) || throw(DimensionMismatch(
        "Reconstructed model prfl_modl_fit size $(size(payload.prfl_modl_fit)) must match " *
        "(x_modl, istp, IB) $((length(x_modl), length(val_istp), length(val_IB))).",
    ))
    return Float64.(payload.prfl_modl_fit)
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
        lines!(ax, x_modl, @view(tail_cohr.tail_istp_avg[:, idx_istp]); color=(:seagreen4, 0.80), linewidth=2.2)
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
        !is_bottom_row && hidexdecorations!(ax; grid=false)
    end
    colgap!(fig.layout, col_start, 0)
    return axs
end

if use_src_profiles
    save_src_profiles = true
    include(joinpath(@__DIR__, "anlz_ssntfr_src.jl"))
else
    x_dens, x_modl, val_IB, ntfr2d_mean, prfl_inco, prfl_cohr = h5open(path_data, "r") do file
        x_dens = read(file["x_dens"])
        x_modl = read(file["x_modl"])
        val_IB = read(file["val_IB"])
        ntfr2d_mean = orient_ntfr2d_axes(read(file["ntfr2d_mean"]), x_dens, val_IB, val_istp)
        prfl_inco = orient_prfl_axes(read(file["prfl_inco"]), x_modl, val_IB, val_istp)
        prfl_cohr = orient_prfl_axes(read(file["prfl_cohr"]), x_modl, val_IB, val_istp)
        return x_dens, x_modl, val_IB, ntfr2d_mean, prfl_inco, prfl_cohr
    end
end
step_modl = median(diff(x_modl))

colorrange_inco = calc_prfl_colorrange(prfl_inco, x_modl, range_x_colorrange)
colorrange_cohr = calc_prfl_colorrange(prfl_cohr, x_modl, range_x_colorrange)

prfl_modl_fit = load_prfl_modl_fit_jld2(path_model_results, x_modl, val_IB, val_istp)
tail_inco = calc_modl_tail(x_modl, prfl_inco, prfl_modl_fit)
tail_cohr = calc_modl_tail(x_modl, prfl_cohr, prfl_modl_fit)
prfl_tail = build_profile_variant_arrays(
    prfl_inco,
    prfl_cohr,
    tail_inco.prfl_tailess,
    tail_cohr.prfl_tailess,
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
        "%s profile tail removal, inco/cohr tails from reconstructed 2D NTFR model, Δk=%.5f",
        tag,
        step_modl,
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
        tail_inco.prfl_tailess;
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
        tail_cohr.prfl_tailess;
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
    text="per-istp mean over IB, gray: total mean, green: reconstructed 2D NTFR model profile",
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
    ylims=(1e-4,1e2)
    # ylims=calc_log_ylims(colorrange_cohr),
)
rowgap!(fig_cohr_avg.layout, 1, 4)
rowgap!(fig_cohr_avg.layout, 2, 4)
colgap!(fig_cohr_avg.layout, 1, 0)

path_plot_cohr_avg = joinpath(path_output, "$(tag)_prfl_cohr_average.png")
save(path_plot_cohr_avg, fig_cohr_avg; backend=CairoMakie)
println("saved $path_plot_cohr_avg")
