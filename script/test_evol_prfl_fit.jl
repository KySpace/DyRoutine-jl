using Interpolations
using Optim
using GLMakie
using JLD2
using Colors
using Pipe: @pipe
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "fitmodels.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, band!, save, text!, rowgap!, colgap!, Button, hspan!, Observable, Label, GridLayout
using LsqFit: curve_fit, coef
using Statistics: mean, std
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
title_extr = raw"[07.24].100.Extr"
title_corr = raw"[07.24].103.PrflProc.NumberMask.[←100]"
cache_extr = JLD2.load(joinpath(path_root, "AnlzRoutine", title_extr, "CFNM_essn_extr.jld2"))
prfl_axial_evol = JLD2.load(joinpath(path_root, "AnlzRoutine", title_corr, "CFNM_corr.jld2"))["prfl_axial_evol_stacked"]
meta_extr = cache_extr["meta_extr"]
val_vars = meta_extr.val_vars
px_in_um = meta_extr.px_in_um
coor = smwh_core[2] |> s -> (-s:s)*px_in_um

function fit_peak_trajectory(
    model,
    evol_prfl::AbstractMatrix,
    r_blur::Real,
    x::AbstractRange,
    t::AbstractVector,
    p0::AbstractVector;
    optimizer = NelderMead(),
    iterations::Integer = 10_000,
)
    nx, nt = size(evol_prfl)
    nx == length(x) || throw(DimensionMismatch("size(evol_prfl, 1) = $nx, but length(x) = $(length(x))"))
    nt == length(t) || throw(DimensionMismatch("size(evol_prfl, 2) = $nt, but length(t) = $(length(t))"))

    r_blur >= 0 || throw(ArgumentError("r_blur must be nonnegative"))

    σ_px = r_blur / abs(step(x))
    kern = KernelFactors.gaussian(σ_px)
    evol_prfl_blur = mapslices(
        profile -> imfilter(profile, kern, Pad(:replicate)),
        float.(evol_prfl);
        dims=1,
    )
    interpolant =
        @pipe interpolate(evol_prfl_blur, (BSpline(Linear()), NoInterp())) |> scale(_, x, axes(evol_prfl_blur, 2))
    xmin, xmax = extrema(x)
    function score(p)
        traj = [model(ti, p) for ti in t]
        all(xi -> xmin <= xi <= xmax, traj) || return -Inf
        sum(i -> interpolant(traj[i], i), eachindex(t))
    end

    result = optimize(
        p -> -score(p),
        collect(p0),
        optimizer,
        Optim.Options(iterations=iterations),
    )

    params = Optim.minimizer(result)
    trajectory = [model(ti, params) for ti in t]

    return (; params, trajectory, evol_prfl_blur, result)
end

ib = findfirst(==(5.316), val_vars.IB) |> something
istp = 2
evol_prfl_sample = prfl_axial_evol[ib, istp]
function oscl_decay(t, p)
    x0, τ_oscl, t0, τ_decay = p
    @. x0 * exp(- t / τ_decay) * cos(2π * (t - t0) / τ_oscl)
end

t_evol = val_vars.t_hold
mask_t = 10 .<= t_evol .<= 120
t_evol_sel = t_evol[mask_t]
evol_prfl_sample_sel = evol_prfl_sample[:,mask_t]
params_traj, traj_recon, evol_prfl_sample_sel_blur, result_fit_traj = fit_peak_trajectory(
    oscl_decay, evol_prfl_sample_sel, 1.0, coor, t_evol_sel,
    [8.0, 40, 10, 200]
)


fig = Figure()
ax_prfl = Axis(fig[1, 1])

clr = make_masked_prfl_clr(evol_prfl_sample_sel_blur, trues(size(t_evol_sel)), hue_theme_istp[val_istp[istp]], 2.5)
heatmap!(ax_prfl, t_evol_sel, coor, clr')
lines!(ax_prfl, t_evol_sel, traj_recon; color=:white)
fig |> resize_to_layout!
fig |> display