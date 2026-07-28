using Interpolations
using Optim
using GLMakie
using JLD2
using Colors
using Pipe: @pipe
using ImageFiltering
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
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, band!, save, text!, rowgap!, colgap!, Button, hspan!, Observable, Label, GridLayout, Slider, IntervalSlider, linkxaxes!, xlims!, rowsize!, resize_to_layout!
using LsqFit: curve_fit, coef
using Statistics: mean, std
using GeometryBasics: Point2f
using Printf: @sprintf
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
title_extr = raw"[07.24].100.Extr"
title_corr = raw"[07.24].103.PrflProc.NumberMask.[←100]"
cache_extr = JLD2.load(joinpath(path_root, "AnlzRoutine", title_extr, "CFNM_essn_extr.jld2"))
prfl_axial_evol = JLD2.load(joinpath(path_root, "AnlzRoutine", title_corr, "CFNM_corr.jld2"))["prfl_axial_evol_stacked"]
meta_extr = cache_extr["meta_extr"]
val_vars = meta_extr.val_vars
px_in_um = meta_extr.px_in_um
coor = meta_extr.smwh_core[2] |> s -> (-s:s)*px_in_um

function fit_peak_trajectory(
    model,
    evol_prfl::AbstractMatrix,
    r_blur::Real,
    x::AbstractRange,
    t::AbstractVector,
    p0::AbstractVector;
    optimizer = ParticleSwarm(),
    iterations::Integer = 10_000,
)
    nx, nt = size(evol_prfl)
    nx == length(x) || throw(DimensionMismatch("size(evol_prfl, 1) = $nx, but length(x) = $(length(x))"))
    nt == length(t) || throw(DimensionMismatch("size(evol_prfl, 2) = $nt, but length(t) = $(length(t))"))

    r_blur >= 0 || throw(ArgumentError("r_blur must be nonnegative"))

    evol_prfl_blur = blur_profile_evolution(evol_prfl, r_blur, x)
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

function blur_profile_evolution(
    evol_prfl::AbstractMatrix,
    r_blur::Real,
    x::AbstractRange,
)
    r_blur >= 0 || throw(ArgumentError("r_blur must be nonnegative"))
    σ_px = r_blur / abs(step(x))
    iszero(σ_px) && return Float64.(evol_prfl)
    kern = KernelFactors.gaussian(σ_px)
    return mapslices(
        profile -> imfilter(profile, kern, Pad(:replicate)),
        float.(evol_prfl);
        dims=1,
    )
end

function oscl_decay(t, p)
    x0, τ_oscl, t0, τ_decay = p
    @. x0 * exp(- t / τ_decay) * cos(2π * (t - t0) / τ_oscl)
end

t_evol = val_vars.t_hold

function shift_profile_evolution(
    evol_prfl::AbstractMatrix,
    trajectory::AbstractVector,
    x::AbstractVector,
)
    size(evol_prfl, 1) == length(x) || throw(DimensionMismatch(
        "profile position dimension $(size(evol_prfl, 1)) does not match x length $(length(x))"
    ))
    size(evol_prfl, 2) == length(trajectory) || throw(DimensionMismatch(
        "profile time dimension $(size(evol_prfl, 2)) does not match trajectory length $(length(trajectory))"
    ))

    shifted = similar(float.(evol_prfl))
    for idx_t in axes(evol_prfl, 2)
        # The displayed coordinate is x - trajectory(t), so sample the
        # original profile at x + trajectory(t) on the fixed plotting grid.
        interpolant = Interpolations.linear_interpolation(
            x,
            Float64.(evol_prfl[:, idx_t]);
            extrapolation_bc=0.0,
        )
        shifted[:, idx_t] .= interpolant.(x .+ trajectory[idx_t])
    end
    return shifted
end

ib0 = findfirst(==(5.316), val_vars.IB) |> something
istp_values = val_vars.istp
fit_initial_range = (10.0, 150.0)
fit_initial_blur = 1.0
fit_initial_t0 = 10.0
fit_colorrange = 2.5

function make_fit(
    ib::Int,
    rng_t_fit::Tuple{<:Real,<:Real},
    r_blur_um::Real,
    t0_initial::Real,
)
    lo, hi = extrema(Float64.(rng_t_fit))
    mask_t = (lo .<= t_evol) .& (t_evol .<= hi)
    any(mask_t) || throw(ArgumentError("fit time range $rng_t_fit selected no t_hold values"))
    t_evol_sel = t_evol[mask_t]
    evol_prfl_sample = mean(prfl_axial_evol[ib, :])
    evol_prfl_blur_full = blur_profile_evolution(evol_prfl_sample, r_blur_um, coor)
    evol_prfl_sample_sel = evol_prfl_sample[:, mask_t]
    fit = fit_peak_trajectory(
        oscl_decay,
        evol_prfl_sample_sel,
        r_blur_um,
        coor,
        t_evol_sel,
        [8.0, 40, t0_initial, 200];
        optimizer=ParticleSwarm(),
    )
    return merge(fit, (;
        t_evol_sel,
        mask_t,
        evol_prfl_blur_full,
        trajectory_full=oscl_decay.(t_evol, Ref(fit.params)),
    ))
end

obs_ib = Observable(ib0)
obs_r_blur_um = Observable(fit_initial_blur)
obs_rng_t_fit = Observable(fit_initial_range)
obs_t0_initial = Observable(fit_initial_t0)
fit_live = Ref(make_fit(obs_ib[], obs_rng_t_fit[], obs_r_blur_um[], obs_t0_initial[]))

fig = Figure(size=(1500, 920))
ax_prfl = Axis(fig[1, 1], ylabel="position (μm)", title="fitted trajectory")
ax_prfl_162 = Axis(fig[2, 1], ylabel="position (μm)", title="162")
ax_prfl_164 = Axis(fig[3, 1], xlabel="t_hold (ms)", ylabel="position (μm)", title="164")

fit_now = fit_live[]
obs_prfl_clr = Observable(make_masked_prfl_clr(
    fit_now.evol_prfl_blur_full,
    fit_now.mask_t,
    132,
    fit_colorrange,
)')
heatmap!(ax_prfl, t_evol, coor, obs_prfl_clr)
traj_plot = lines!(ax_prfl, fit_now.t_evol_sel, fit_now.trajectory; color=:white, linewidth=2)

obs_prfl_162 = Observable(shift_profile_evolution(
    prfl_axial_evol[ib0, 1], fit_now.trajectory_full, coor,
)')
obs_prfl_164 = Observable(shift_profile_evolution(
    prfl_axial_evol[ib0, 2], fit_now.trajectory_full, coor,
)')
obs_colorrange_162 = Observable(calc_prfl_colorrange_auto(
    obs_prfl_162[]', t_evol, coor,
))
obs_colorrange_164 = Observable(calc_prfl_colorrange_auto(
    obs_prfl_164[]', t_evol, coor,
))
heatmap!(ax_prfl_162, t_evol, coor, obs_prfl_162;
    colormap=gen_clrmap_solo(hue_theme_istp[string(istp_values[1])]),
    colorrange=obs_colorrange_162,
)
heatmap!(ax_prfl_164, t_evol, coor, obs_prfl_164;
    colormap=gen_clrmap_solo(hue_theme_istp[string(istp_values[2])]),
    colorrange=obs_colorrange_164,
)
linkxaxes!([ax_prfl, ax_prfl_162, ax_prfl_164])
xlims!(ax_prfl_164, (0, 200))

fit_caption = Observable("")
text!(ax_prfl, fit_caption; position=(0.98, 0.98), space=:relative,
    align=(:right, :top), justification=:right, fontsize=14)

ctrl = GridLayout(fig[4, 1])
function add_cycle!(col, obs, values)
    prev = Button(ctrl[1, col]; label="←")
    Label(ctrl[1, col + 1], lift(i -> "IB=$(values[i]) A\n$(i)/$(length(values))", obs))
    next = Button(ctrl[1, col + 2]; label="→")
    on(prev.clicks) do _
        obs[] = mod1(obs[] - 1, length(values))
    end
    on(next.clicks) do _
        obs[] = mod1(obs[] + 1, length(values))
    end
end
add_cycle!(1, obs_ib, val_vars.IB)

Label(ctrl[1, 5], lift(r -> "blur $(round(r; digits=2)) μm", obs_r_blur_um))
slider_blur = Slider(ctrl[1, 6:9]; range=0.0:0.1:4.0, startvalue=fit_initial_blur, tellwidth=false)
Label(ctrl[1, 10], lift(r -> "fit $(r[1])–$(r[2]) ms", obs_rng_t_fit))
slider_fit_range = IntervalSlider(ctrl[1, 11:16]; range=0.0:1.0:200.0, startvalues=fit_initial_range)
Label(ctrl[1, 17], lift(t0 -> "initial t₀ $(round(t0; digits=1)) ms", obs_t0_initial))
slider_t0_initial = Slider(ctrl[1, 18:21]; range=-15.0:0.5:15.0, startvalue=fit_initial_t0, tellwidth=false)

on(slider_blur.value) do r_blur
    obs_r_blur_um[] = Float64(r_blur)
end
on(slider_fit_range.interval) do rng_t_fit
    obs_rng_t_fit[] = (Float64(rng_t_fit[1]), Float64(rng_t_fit[2]))
end
on(slider_t0_initial.value) do t0_initial
    obs_t0_initial[] = Float64(t0_initial)
end

function update_plot!()
    fit_now = make_fit(obs_ib[], obs_rng_t_fit[], obs_r_blur_um[], obs_t0_initial[])
    fit_live[] = fit_now
    obs_prfl_clr[] = make_masked_prfl_clr(
        fit_now.evol_prfl_blur_full,
        fit_now.mask_t,
        132,
        fit_colorrange,
    )'
    traj_plot[1][] = Point2f.(fit_now.t_evol_sel, fit_now.trajectory)
    obs_prfl_162[] = shift_profile_evolution(
        prfl_axial_evol[obs_ib[], 1], fit_now.trajectory_full, coor,
    )'
    obs_prfl_164[] = shift_profile_evolution(
        prfl_axial_evol[obs_ib[], 2], fit_now.trajectory_full, coor,
    )'
    obs_colorrange_162[] = calc_prfl_colorrange_auto(obs_prfl_162[]', t_evol, coor)
    obs_colorrange_164[] = calc_prfl_colorrange_auto(obs_prfl_164[]', t_evol, coor)
    x0, τ_oscl, t0, τ_decay = fit_now.params
    fit_caption[] = "x₀=$(@sprintf("%.2f", x0)), τₒ=$(@sprintf("%.2f", τ_oscl))\n" *
        "t₀=$(@sprintf("%.2f", t0)), τd=$(@sprintf("%.2f", τ_decay))"
end
onany(obs_ib, obs_r_blur_um, obs_rng_t_fit, obs_t0_initial) do _...
    update_plot!()
end
update_plot!()

rowsize!(fig.layout, 1, 240)
rowsize!(fig.layout, 2, 240)
rowsize!(fig.layout, 3, 240)
fig |> resize_to_layout!
fig |> display
