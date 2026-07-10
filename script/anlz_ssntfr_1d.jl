using CairoMakie
using FFTW
using HDF5
using ImageFiltering
using JLD2
using LinearAlgebra
using LsqFit
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
title_anlz = "35.Ntfr1D.Abrr.Beta.MovingTail.NoMuOffset.Log"
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
fit_center_bound = 5.0
fit_maxiter_1d = parse(Int, get(ENV, "SSNTFR_1D_MAXITER", "30000"))
fit_threshold_1d = 1.5e-1
fit_sigma_wide_min = 15.0
fit_log_1d = true
smwh_source = (150, 150)
smwh_reconstruct = (150, 0)
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

function pack_fit_profile_results(fit_density::AbstractMatrix, val_IB::AbstractVector, val_istp::AbstractVector)
    size(fit_density) == (length(val_IB), length(val_istp)) || throw(DimensionMismatch(
        "fit_density size $(size(fit_density)) must match (IB, istp) $((length(val_IB), length(val_istp))).",
    ))
    n_params = length(fit_density[1, 1].params)
    n_dens = length(fit_density[1, 1].fit)
    params = Array{Float64}(undef, n_params, length(val_istp), length(val_IB))
    fit_profile = Array{Float64}(undef, n_dens, length(val_istp), length(val_IB))
    fit_tail = similar(fit_profile)
    fit_narrow = similar(fit_profile)
    fit_tailess = similar(fit_profile)
    fit_tailess_beta0 = similar(fit_profile)
    rss_rel = Matrix{Float64}(undef, length(val_istp), length(val_IB))
    maxiter_reached = Matrix{Bool}(undef, length(val_istp), length(val_IB))
    for idx_IB in eachindex(val_IB), idx_istp in eachindex(val_istp)
        fit = fit_density[idx_IB, idx_istp]
        params[:, idx_istp, idx_IB] .= fit.params
        fit_profile[:, idx_istp, idx_IB] .= fit.fit
        fit_tail[:, idx_istp, idx_IB] .= fit.tail
        fit_narrow[:, idx_istp, idx_IB] .= fit.narrow
        fit_tailess[:, idx_istp, idx_IB] .= fit.tailess
        fit_tailess_beta0[:, idx_istp, idx_IB] .= fit.tailess_beta0
        rss_rel[idx_istp, idx_IB] = fit.rss_rel
        maxiter_reached[idx_istp, idx_IB] = fit.maxiter_reached
    end
    return (; params, fit_profile, fit_tail, fit_narrow, fit_tailess, fit_tailess_beta0, rss_rel, maxiter_reached)
end

function build_model_results_payload(;
    x_dens,
    x_modl,
    val_IB,
    val_istp,
    fit_density,
    prfl_modl_fit,
    fit_maxiter_1d,
    fit_log_1d,
    fit_threshold_1d,
    fit_sigma_wide_min,
    fit_center_bound,
    smwh_source,
    smwh_reconstruct,
    step_modl,
)
    fit_results = pack_fit_profile_results(fit_density, val_IB, val_istp)
    prfl_modl_fit_packed = pack_prfl_modl_fit(prfl_modl_fit, x_modl, val_IB, val_istp)
    meta = (;
        fit_param_names=[
            "x0",
            "A_narrow",
            "sigma_narrow",
            "skew_narrow",
            "beta_narrow",
            "A_wide",
            "sigma_wide",
            "skew_wide",
            "beta_wide",
        ],
        dim_fit_params="param,istp,IB",
        dim_fit_profile="x_dens,istp,IB",
        dim_fit_rss_rel="istp,IB",
        dim_prfl_modl_fit="x_modl,istp,IB",
        fit_maxiter_1d,
        fit_log_1d,
        fit_threshold_1d,
        fit_sigma_wide_min,
        fit_center_bound,
        smwh_source,
        smwh_reconstruct,
        step_modl,
    )
    return (;
        x_dens,
        x_modl,
        val_IB,
        val_istp=String.(val_istp),
        fit_params=fit_results.params,
        fit_profile=fit_results.fit_profile,
        fit_tail=fit_results.fit_tail,
        fit_narrow=fit_results.fit_narrow,
        fit_tailess=fit_results.fit_tailess,
        fit_tailess_beta0=fit_results.fit_tailess_beta0,
        fit_rss_rel=fit_results.rss_rel,
        fit_maxiter_reached=fit_results.maxiter_reached,
        prfl_modl_fit=prfl_modl_fit_packed,
        meta,
    )
end

function save_model_results_h5(path_results::AbstractString, payload)
    h5open(path_results, "w") do file
        write(file, "x_dens", payload.x_dens)
        write(file, "x_modl", payload.x_modl)
        write(file, "val_IB", payload.val_IB)
        write(file, "val_istp", payload.val_istp)
        write(file, "fit_params", payload.fit_params)
        write(file, "fit_profile", payload.fit_profile)
        write(file, "fit_tail", payload.fit_tail)
        write(file, "fit_narrow", payload.fit_narrow)
        write(file, "fit_tailess", payload.fit_tailess)
        write(file, "fit_tailess_beta0", payload.fit_tailess_beta0)
        write(file, "fit_rss_rel", payload.fit_rss_rel)
        write(file, "fit_maxiter_reached", payload.fit_maxiter_reached)
        write(file, "prfl_modl_fit", payload.prfl_modl_fit)
        attrs(file)["fit_param_names"] = join(payload.meta.fit_param_names, ",")
        attrs(file)["dim_fit_params"] = payload.meta.dim_fit_params
        attrs(file)["dim_fit_profile"] = payload.meta.dim_fit_profile
        attrs(file)["dim_fit_rss_rel"] = payload.meta.dim_fit_rss_rel
        attrs(file)["dim_prfl_modl_fit"] = payload.meta.dim_prfl_modl_fit
        attrs(file)["fit_maxiter_1d"] = payload.meta.fit_maxiter_1d
        attrs(file)["fit_log_1d"] = payload.meta.fit_log_1d
        attrs(file)["fit_threshold_1d"] = payload.meta.fit_threshold_1d
        attrs(file)["fit_sigma_wide_min"] = payload.meta.fit_sigma_wide_min
        attrs(file)["fit_center_bound"] = payload.meta.fit_center_bound
        attrs(file)["smwh_source"] = collect(payload.meta.smwh_source)
        attrs(file)["smwh_reconstruct"] = collect(payload.meta.smwh_reconstruct)
        attrs(file)["step_modl"] = payload.meta.step_modl
    end
    return path_results
end

function save_model_results_jld2(path_results::AbstractString, payload)
    JLD2.jldopen(path_results, "w") do file
        file["x_dens"] = payload.x_dens
        file["x_modl"] = payload.x_modl
        file["val_IB"] = payload.val_IB
        file["val_istp"] = payload.val_istp
        file["fit_params"] = payload.fit_params
        file["fit_profile"] = payload.fit_profile
        file["fit_tail"] = payload.fit_tail
        file["fit_narrow"] = payload.fit_narrow
        file["fit_tailess"] = payload.fit_tailess
        file["fit_tailess_beta0"] = payload.fit_tailess_beta0
        file["fit_rss_rel"] = payload.fit_rss_rel
        file["fit_maxiter_reached"] = payload.fit_maxiter_reached
        file["prfl_modl_fit"] = payload.prfl_modl_fit
        file["meta"] = payload.meta
    end
    return path_results
end

function draw_density_row!(
    fig::Figure,
    row::Integer,
    x_dens::AbstractVector{<:Real},
    x_modl::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    idx_IB::Integer,
    IB::Real;
    ntfr2d_mean,
    colorrange_ntfr,
    fit_density,
    xlims_profile,
    xlims_folded,
    ylims_diag,
    ylims_profile,
    prfl_modl_fit,
    prfl_inco,
    prfl_cohr,
    xlims_prfl,
    ylims_prfl,
    idx_col_IB_right::Integer,
    is_bottom_row::Bool=false,
)
    Label(fig[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:right)
    Label(fig[row, idx_col_IB_right]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:left)

    for (idx_istp, istp) in enumerate(val_istp)
        clr_strong = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.95)
        clr_faint = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.32)
        clr_theme_shade = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.20)
        clr_theme_shade_dark = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.26)
        clr_center = Oklch(0.62, 0.16, 145)
        fit_info = fit_density[idx_IB, idx_istp]
        params = fit_info.params

        ax_dens = Axis(
            fig[row, idx_istp];
            xlabel=is_bottom_row ? label_x_dens : "",
            ylabel=idx_istp == 1 ? label_x_dens : "",
            aspect=DataAspect(),
        )
        clrmap = gen_clrmap_solo(hue_theme_istp[istp])
        heatmap!(ax_dens, x_dens, x_dens, ntfr2d_mean[idx_IB, idx_istp]'; colormap=clrmap, colorrange=colorrange_ntfr, rasterize=true)
        hidexdecorations!(ax_dens; label=is_bottom_row ? false : true, ticklabels=is_bottom_row ? false : true, ticks=is_bottom_row ? false : true, grid=false)
        hideydecorations!(ax_dens; label=idx_istp == 1 ? false : true, ticklabels=false, ticks=false, grid=false)

        idx_col_profile = length(val_istp) + idx_istp
        ax_profile = Axis(
            fig[row, idx_col_profile];
            xlabel=is_bottom_row ? label_x_dens : "",
            ylabel=idx_istp == 1 ? "row strip" : "",
            yaxisposition=idx_istp == 1 ? :left : :right,
            xticks=-40:20:40,
        )
        band!(ax_profile, x_dens, zero.(fit_info.tail), fit_info.tail; color=(:gray20, 0.18))
        band!(ax_profile, x_dens, zero.(fit_info.narrow_beta0), fit_info.narrow_beta0; color=clr_theme_shade_dark)
        band!(ax_profile, x_dens, zero.(fit_info.narrow), fit_info.narrow; color=clr_theme_shade)
        lines!(ax_profile, x_dens, fit_info.profile; color=clr_faint, linewidth=1.0)
        lines!(ax_profile, x_dens, fit_info.tail; color=(:gray20, 0.55), linewidth=1.0)
        lines!(ax_profile, x_dens, fit_info.tailess_beta0; color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.48), linewidth=1.25)
        lines!(ax_profile, x_dens, fit_info.tailess; color=clr_strong, linewidth=1.8)
        lines!(ax_profile, x_dens, fit_info.fit; color=clr_center, linewidth=1.4)
        vlines!(ax_profile, [-r_tail_min_profile, r_tail_min_profile]; color=(:gray20, 0.28), linewidth=0.7)
        xlims!(ax_profile, xlims_profile)
        ylims!(ax_profile, ylims_profile)
        text!(
            ax_profile,
            0.05, 0.95;
            text=@sprintf("σN=%.2f αN=%.2f βN=%.2f\nσT=%.1f αT=%.2f βT=%.2f", params[3], params[4], params[5], params[7], params[8], params[9]),
            space=:relative,
            color=(clr_center, 0.9),
            fontsize=8,
            align=(:left, :top),
        )
        !is_bottom_row && hidexdecorations!(ax_profile; grid=false)

        idx_col_diag = 2 * length(val_istp) + idx_istp
        ax_diag = Axis(
            fig[row, idx_col_diag];
            xlabel=is_bottom_row ? label_x_dens : "",
            ylabel=idx_istp == 1 ? "folded log" : "",
            yscale=log10,
            yaxisposition=idx_istp == 1 ? :left : :right,
            xticks=0:20:40,
        )
        clr_tail_fit = RGBAf(Oklch(0.58, 0.17, 145), 0.88)
        draw_folded_branch!(ax_diag, x_dens, fit_info.profile, :pos; color=clr_strong, linewidth=1.5)
        draw_folded_branch!(ax_diag, x_dens, fit_info.profile, :neg; color=clr_strong, linewidth=1.2)
        draw_folded_branch!(ax_diag, x_dens, fit_info.tail, :pos; color=clr_tail_fit, linewidth=1.5)
        draw_folded_branch!(ax_diag, x_dens, fit_info.tail, :neg; color=clr_tail_fit, linewidth=1.2)
        hlines!(ax_diag, fit_threshold_1d; color=(:gray20, 0.45), linewidth=0.8)
        xlims!(ax_diag, xlims_folded)
        ylims!(ax_diag, ylims_diag)
        text!(
            ax_diag,
            0.05, 0.05;
            text=@sprintf("σT=%.1f μm", params[7]),
            space=:relative,
            color=clr_tail_fit,
            fontsize=8,
            align=(:left, :bottom),
        )
        !is_bottom_row && hidexdecorations!(ax_diag; grid=false)

        idx_col_prfl = 3 * length(val_istp) + idx_istp
        ax_prfl = Axis(
            fig[row, idx_col_prfl];
            xlabel=is_bottom_row ? label_x_modl : "",
            ylabel=idx_istp == 1 ? "FT profile" : "",
            yaxisposition=idx_istp == 1 ? :left : :right,
            xticks=0:0.2:0.6,
        )
        clr_theme = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.92)
        clr_theme_faint = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.60)
        clr_fit = RGBAf(Oklch(0.60, 0.17, 145), 0.95)
        clr_rss = fit_info.maxiter_reached ? RGBAf(Oklch(0.50, 0.12, 65), 0.95) : clr_fit
        lines!(ax_prfl, x_modl, prfl_modl_fit[idx_IB, idx_istp]; color=clr_fit, linewidth=1.8)
        lines!(ax_prfl, x_modl, @view(prfl_inco[:, idx_istp, idx_IB]); color=clr_theme_faint, linewidth=1.3, linestyle=:dash)
        lines!(ax_prfl, x_modl, @view(prfl_cohr[:, idx_istp, idx_IB]); color=clr_theme, linewidth=1.5)
        xlims!(ax_prfl, xlims_prfl)
        ylims!(ax_prfl, ylims_prfl)
        text!(
            ax_prfl,
            0.95, 0.95;
            text=@sprintf("rss_rel=%.3f", fit_info.rss_rel),
            space=:relative,
            color=clr_rss,
            fontsize=8,
            align=(:right, :top),
        )
        !is_bottom_row && hidexdecorations!(ax_prfl; grid=false)
    end
    return nothing
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
cfg_prfl = get_prfl_modl_1d_config(smwh_source)

t_stage = log_step("fitting centered 1D NTFR row-strip densities")
fit_density = Array{NamedTuple}(undef, size(ntfr2d_mean))
for idx_IB in axes(ntfr2d_mean, 1), idx_istp in axes(ntfr2d_mean, 2)
    println("  [$tag] fitting 1D density IB_idx=$idx_IB IB=$(val_IB[idx_IB]) istp=$(val_istp[idx_istp])")
    flush(stdout)
    profile = calc_center_row_strip_profile(
        x_dens,
        ntfr2d_mean[idx_IB, idx_istp];
        smh_dens_strip=cfg_prfl.smh_dens_strip,
    )
    fit_1d = fit_skew_beta_gauss_1d(
        x_dens,
        profile;
        center_bound=fit_center_bound,
        threshold=fit_threshold_1d,
        sigma_wide_min=fit_sigma_wide_min,
        r_narrow_max=4*8.0, # r_tail_min_profile,
        maxiter=fit_maxiter_1d,
        fit_log=fit_log_1d,
    )
    params = fit_1d.params
    println(
        "  [$tag] fit done IB_idx=$idx_IB istp_idx=$idx_istp " *
        "rss=$(round(fit_1d.rss_rel; digits=4)) " *
        (fit_1d.maxiter_reached ? "maxiter=true " : "") *
        (fit_1d.fit_log ? "log=true " : "") *
        "αN=$(round(params[4]; digits=4)) βN=$(round(params[5]; digits=4)) " *
        "αT=$(round(params[8]; digits=4)) βT=$(round(params[9]; digits=4))",
    )
    flush(stdout)
    fit_density[idx_IB, idx_istp] = fit_1d
end
log_done("fit centered 1D NTFR row-strip densities", t_stage)

t_stage = log_step("reconstructing fitted modulation profiles")
prfl_modl_fit = map(fit_density) do fit
    dens_strip = reshape(fit.fit, 1, :)
    calc_prfl_modl_1d([dens_strip], smwh_reconstruct; step_modl).prfl_inco
end
length(x_modl) == length(prfl_modl_fit[1]) || throw(DimensionMismatch(
    "x_modl length $(length(x_modl)) must match reconstructed profile length $(length(prfl_modl_fit[1])).",
))
log_done("reconstructed fitted modulation profiles", t_stage)

t_stage = log_step("building plot ranges")
colorrange_ntfr = (0.0, maximum(maximum, ntfr2d_mean))
val_profile = [v for fit in fit_density for v in fit.profile if isfinite(v)]
val_tailess = [
    v
    for fit in fit_density
    for profile in (fit.tailess, fit.tailess_beta0)
    for v in profile
    if isfinite(v)
]
val_profile_positive = [v for v in val_profile if v > 0]
xlims_profile = (minimum(x_dens), maximum(x_dens))
xlims_folded = (0.0, maximum(abs, x_dens))
ylims_profile = (min(0.0, minimum(val_tailess) * 1.05), maximum(val_profile) * 1.05)
ylims_diag = (1e-2, maximum(val_profile_positive) * 1.1)
log_done("built plot ranges", t_stage)

t_stage = log_step("building figure axes and labels")
fig_ntfr = Figure(fontsize=14)
idx_col_IB_right = 4 * length(val_istp) + 1
Label(
    fig_ntfr[0, 1:idx_col_IB_right];
    text=@sprintf(
        "%s log-fit 1D NTFR row-strip fits, center %d-row mean |> skew/beta narrow + skew/beta tail, mask > %.1g, σ_wide ≥ %.0f μm, Δk=%.5f",
        tag,
        2 * cfg_prfl.smh_dens_strip + 1,
        fit_threshold_1d,
        fit_sigma_wide_min,
        step_modl,
    ),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_ntfr[1, idx_istp]; text="density $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_ntfr[1, length(val_istp) + idx_istp]; text="row tailess $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_ntfr[1, 2 * length(val_istp) + idx_istp]; text="folded log $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
for (idx_istp, istp) in enumerate(val_istp)
    Label(fig_ntfr[1, 3 * length(val_istp) + idx_istp]; text="FT prfl $istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
end
Label(fig_ntfr[1, idx_col_IB_right]; text="IB", tellwidth=true, tellheight=true, halign=:left, font=:bold)
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
        idx_IB,
        IB;
        ntfr2d_mean,
        colorrange_ntfr,
        fit_density,
        xlims_profile,
        xlims_folded,
        ylims_diag,
        ylims_profile,
        prfl_modl_fit,
        prfl_inco,
        prfl_cohr,
        xlims_prfl=xlims_prfl_reconstruct,
        ylims_prfl=ylims_prfl_reconstruct,
        idx_col_IB_right,
        is_bottom_row=idx_IB == length(val_IB),
    )
    rowsize!(fig_ntfr.layout, row, Fixed(105))
end
log_done("plotted density rows", t_stage)

for col in 1:length(val_istp)
    colsize!(fig_ntfr.layout, col, Fixed(105))
end
for col in length(val_istp)+1:idx_col_IB_right-1
    colsize!(fig_ntfr.layout, col, Fixed(170))
end
colsize!(fig_ntfr.layout, idx_col_IB_right, Fixed(55))
for col in 1:idx_col_IB_right-1
    colgap!(fig_ntfr.layout, col, isodd(col) ? 0 : 14)
end
colgap!(fig_ntfr.layout, length(val_istp), 16)
colgap!(fig_ntfr.layout, idx_col_IB_right - 1, 8)
rowgap!(fig_ntfr.layout, 1, 4)

resize_to_layout!(fig_ntfr)
filename_plot_ntfr = "$(tag)_ntfr1d_table"
t_stage = log_step("saving figure outputs")
for ext in ("png", "pdf")
    save(joinpath(path_output, "$filename_plot_ntfr.$ext"), fig_ntfr; backend=CairoMakie)
end
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
cp(joinpath(@__DIR__, "..", "src", "modlntfr.jl"), joinpath(path_output, "modlntfr.jl"); force=true)
log_done("saved figure outputs", t_stage)
println("saved $(joinpath(path_output, "$filename_plot_ntfr.png"))")

filename_model_results = "$(tag)_ntfr1d_model_results.h5"
filename_model_results_jld2 = "$(tag)_ntfr1d_model_results.jld2"
path_model_results = joinpath(path_output, filename_model_results)
path_model_results_jld2 = joinpath(path_output, filename_model_results_jld2)
t_stage = log_step("saving model outputs")
payload_model_results = build_model_results_payload(;
    x_dens,
    x_modl,
    val_IB,
    val_istp,
    fit_density,
    prfl_modl_fit,
    fit_maxiter_1d,
    fit_log_1d,
    fit_threshold_1d,
    fit_sigma_wide_min,
    fit_center_bound,
    smwh_source,
    smwh_reconstruct,
    step_modl,
)
save_model_results_h5(path_model_results, payload_model_results)
save_model_results_jld2(path_model_results_jld2, payload_model_results)
log_done("saved model outputs", t_stage)
println("saved $path_model_results")
println("saved $path_model_results_jld2")
