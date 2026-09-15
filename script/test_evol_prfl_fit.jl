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
cache_corr = JLD2.load(joinpath(path_root, "AnlzRoutine", title_corr, "CFNM_corr.jld2"))
prfl_axial_evol = cache_corr["prfl_axial_evol_stacked"]
prfl_radial_evol = cache_corr["prfl_radial_evol_stacked"]
prfl_axial_evol_reps = cache_corr["prfl_axial_evol"]
prfl_radial_evol_reps = cache_corr["prfl_radial_evol"]
meta_extr = cache_extr["meta_extr"]
val_vars = meta_extr.val_vars
px_in_um = meta_extr.px_in_um
coor_axial = meta_extr.smwh_core[2] |> s -> (-s:s)*px_in_um
coor_radial = meta_extr.smwh_core[1] |> s -> (-s:s)*px_in_um

function fit_peak_trajectory(
    model,
    evol_prfl::AbstractMatrix,
    r_blur::Real,
    x::AbstractRange,
    t::AbstractVector,
    p0::AbstractVector;
    optimizer = ParticleSwarm(),
    iterations::Integer = 10_000,
    p_upper = nothing,
    p_lower = nothing,
)
    nx, nt = size(evol_prfl)
    nx == length(x) || throw(DimensionMismatch("size(evol_prfl, 1) = $nx, but length(x) = $(length(x))"))
    nt == length(t) || throw(DimensionMismatch("size(evol_prfl, 2) = $nt, but length(t) = $(length(t))"))

    r_blur >= 0 || throw(ArgumentError("r_blur must be nonnegative"))

    p_upper = p_upper === nothing ? fill(Inf, size(p0)) : p_upper
    p_lower = p_lower === nothing ? fill(-Inf, size(p0)) : p_lower
    length(p_upper) == length(p0) || throw(DimensionMismatch(
        "length(p_upper) = $(length(p_upper)), but length(p0) = $(length(p0))"
    ))
    length(p_lower) == length(p0) || throw(DimensionMismatch(
        "length(p_lower) = $(length(p_lower)), but length(p0) = $(length(p0))"
    ))
    all(p_lower[i] <= p_upper[i] for i in eachindex(p0)) || throw(ArgumentError(
        "each p_lower value must be less than or equal to the corresponding p_upper value"
    ))

    evol_prfl_blur = blur_profile_evolution(evol_prfl, r_blur, x)
    interpolant =
        @pipe interpolate(evol_prfl_blur, (BSpline(Linear()), NoInterp())) |> scale(_, x, axes(evol_prfl_blur, 2))
    xmin, xmax = extrema(x)
    function score(p)
        all(p_lower[i] <= p[i] <= p_upper[i] for i in eachindex(p0)) || return -Inf
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

    residue = -score(params)
    return (; params, trajectory, evol_prfl_blur, result, residue)
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
    x0, τ_oscl, t0, τ_decay, x_off = p
    @. x0 * exp(- t / τ_decay) * cos(2π * (t - t0) / τ_oscl) + x_off
end

constant_model(t, p) = p[1]

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
fit_initial_range = (10.0, 150.0)
fit_initial_blur = 1.0
fit_initial_t0 = 10.0
profile_modes = (:axial, :radial)
fit_config_axial  = (; p0=[8.0, 33.0, fit_initial_t0, 50.0, 0.0], p_upper=[20.0, 40.0, 20.0, 10_000.0, 5.0], p_lower=[0.0, 25.0, -20.0, 10.0, -5.0])
fit_config_radial = (; p0=[2.0, 13.0, fit_initial_t0, 50.0, 0.0], p_upper=[ 5.0, 15.0, 20.0, 10_000.0, 5.0], p_lower=[0.0, 10.0, -20.0, 10.0, -5.0])

function profile_data(profile_mode::Symbol)
    profile_mode === :axial && return prfl_axial_evol, prfl_axial_evol_reps, coor_axial, fit_config_axial
    profile_mode === :radial && return prfl_radial_evol, prfl_radial_evol_reps, coor_radial, fit_config_radial
    throw(ArgumentError("unknown profile mode: $profile_mode"))
end

function make_fit(
    ib::Int,
    profile_mode::Symbol,
    istp::Int,
    rng_t_fit::Tuple{<:Real,<:Real},
    r_blur_um::Real,
    t0_initial::Real,
)
    lo, hi = extrema(Float64.(rng_t_fit))
    mask_t = (lo .<= t_evol) .& (t_evol .<= hi)
    any(mask_t) || throw(ArgumentError("fit time range $rng_t_fit selected no t_hold values"))
    t_evol_sel = t_evol[mask_t]
    prfl_evol, _, coor, fit_config = profile_data(profile_mode)
    evol_prfl_sample = prfl_evol[ib, istp]
    evol_prfl_blur_full = blur_profile_evolution(evol_prfl_sample, r_blur_um, coor)
    evol_prfl_sample_sel = evol_prfl_sample[:, mask_t]
    fit = fit_peak_trajectory(
        oscl_decay,
        evol_prfl_sample_sel,
        r_blur_um,
        coor,
        t_evol_sel,
        [fit_config.p0[1], fit_config.p0[2], t0_initial, fit_config.p0[4], fit_config.p0[5]];
        optimizer=ParticleSwarm(),
        p_upper=fit_config.p_upper,
        p_lower=fit_config.p_lower,
    )
    return merge(fit, (;
        t_evol_sel,
        mask_t,
        t0_initial,
        evol_prfl_blur_full,
        trajectory_full=oscl_decay.(t_evol, Ref(fit.params)),
    ))
end

function make_fits(
    ib::Int,
    profile_mode::Symbol,
    rng_t_fit::Tuple{<:Real,<:Real},
    r_blur_um::Real,
    t0_initial::Real,
)
    prfl_evol, _, _, _ = profile_data(profile_mode)
    return [
        make_fit(ib, profile_mode, istp, rng_t_fit, r_blur_um, t0_initial)
        for istp in axes(prfl_evol, 2)
    ]
end

function stack_profile_evolution(prfl_evols::AbstractVector{<:AbstractMatrix})
    isempty(prfl_evols) && throw(ArgumentError("cannot stack an empty profile evolution collection"))
    size_first = size(first(prfl_evols))
    all(size(prfl) == size_first for prfl in prfl_evols) || throw(DimensionMismatch(
        "all profile evolutions must have the same size to be stacked"
    ))
    stacked = zeros(Float64, size_first)
    for prfl in prfl_evols
        stacked .+= prfl
    end
    stacked ./= length(prfl_evols)
    return stacked
end

function fit_shifted_profile_stack(
    prfl_evol_reps::AbstractArray{<:Any,3},
    ib::Int,
    istp::Int,
    trajectory::AbstractVector,
    t::AbstractVector,
    coor::AbstractVector,
)
    n_rep = size(prfl_evol_reps, 2)
    shifted_twice = Matrix{Float64}[]
    x_off = Float64[]
    trajectory_offset = zeros(length(trajectory))
    for rep in axes(prfl_evol_reps, 2)
        shifted_once = shift_profile_evolution(
            prfl_evol_reps[ib, rep, istp], trajectory, coor,
        )
        fit_offset = fit_peak_trajectory(
            constant_model,
            shifted_once,
            0.0,
            coor,
            t,
            [0.0];
            optimizer=NelderMead(),
            iterations=1_000,
        )
        offset = fit_offset.params[1]
        push!(x_off, offset)
        trajectory_offset .= -offset
        push!(shifted_twice, shift_profile_evolution(
            shifted_once, trajectory_offset, coor,
        ))
    end
    length(shifted_twice) == n_rep || throw(AssertionError("unexpected repetition count"))
    return (; evolution=stack_profile_evolution(shifted_twice), x_off)
end

obs_ib = Observable(ib0)
obs_profile_mode = Observable(2)
obs_r_blur_um = Observable(fit_initial_blur)
obs_rng_t_fit = Observable(fit_initial_range)
obs_t0_initial = Observable(fit_initial_t0)
fit_live = Ref(make_fits(obs_ib[], profile_modes[obs_profile_mode[]], obs_rng_t_fit[], obs_r_blur_um[], obs_t0_initial[]))

fig = Figure(size=(600, 480))
ax_prfl_162_orig = Axis(fig[1, 1], ylabel="position (μm)", title="162 original")
ax_prfl_162_shift = Axis(fig[2, 1], ylabel="position (μm)", title="162 shifted")
ax_prfl_164_orig = Axis(fig[3, 1], ylabel="position (μm)", title="164 original")
ax_prfl_164_shift = Axis(fig[4, 1], xlabel="t_hold (ms)", ylabel="position (μm)", title="164 shifted")

fit_now = fit_live[]
profile_mode = profile_modes[obs_profile_mode[]]
_, prfl_evol_reps, coor, _ = profile_data(profile_mode)
fit_162, fit_164 = fit_now

function format_fit_caption(fit)
    x0, τ_oscl, t0, τ_decay, x_off = fit.params
    return "x₀=$(@sprintf("%.2f", x0)), τₒ=$(@sprintf("%.2f", τ_oscl))\n" *
        "t₀ᵢ=$(@sprintf("%.2f", fit.t0_initial)) → t₀=$(@sprintf("%.2f", t0))\n" *
        "τd=$(@sprintf("%.2f", τ_decay)), xₒ=$(@sprintf("%.2f", x_off))\n" *
        "res=$(@sprintf("%.6g", fit.residue))"
end

obs_prfl_162_orig = Observable(fit_162.evol_prfl_blur_full')
obs_prfl_164_orig = Observable(fit_164.evol_prfl_blur_full')
obs_colorrange_162_orig = Observable(calc_prfl_colorrange_auto(
    obs_prfl_162_orig[]', t_evol, coor,
))
obs_colorrange_164_orig = Observable(calc_prfl_colorrange_auto(
    obs_prfl_164_orig[]', t_evol, coor,
))
heatmap!(ax_prfl_162_orig, t_evol, coor, obs_prfl_162_orig;
    colormap=gen_clrmap_solo(hue_theme_istp["162"]),
    colorrange=obs_colorrange_162_orig,
)
heatmap!(ax_prfl_164_orig, t_evol, coor, obs_prfl_164_orig;
    colormap=gen_clrmap_solo(hue_theme_istp["164"]),
    colorrange=obs_colorrange_164_orig,
)
traj_plot_162 = lines!(ax_prfl_162_orig, fit_162.t_evol_sel, fit_162.trajectory;
    color=:white, linewidth=2,
)
traj_plot_164 = lines!(ax_prfl_164_orig, fit_164.t_evol_sel, fit_164.trajectory;
    color=:white, linewidth=2,
)
obs_fit_caption_162 = Observable(format_fit_caption(fit_162))
obs_fit_caption_164 = Observable(format_fit_caption(fit_164))
text!(ax_prfl_162_orig, obs_fit_caption_162;
    position=(0.98, 0.98), space=:relative,
    align=(:right, :top), justification=:right, fontsize=11,
)
text!(ax_prfl_164_orig, obs_fit_caption_164;
    position=(0.98, 0.98), space=:relative,
    align=(:right, :top), justification=:right, fontsize=11,
)

fit_shifted_162 = fit_shifted_profile_stack(
    prfl_evol_reps, ib0, 1, fit_162.trajectory_full, t_evol, coor,
)
fit_shifted_164 = fit_shifted_profile_stack(
    prfl_evol_reps, ib0, 2, fit_164.trajectory_full, t_evol, coor,
)
obs_prfl_162_shift = Observable(fit_shifted_162.evolution')
obs_prfl_164_shift = Observable(fit_shifted_164.evolution')
format_x_off(x_off) = join((@sprintf("%.2f", value) for value in x_off), ", ")
obs_x_off_162 = Observable(format_x_off(fit_shifted_162.x_off))
obs_x_off_164 = Observable(format_x_off(fit_shifted_164.x_off))
obs_colorrange_162_shift = Observable(calc_prfl_colorrange_auto(
    obs_prfl_162_shift[]', t_evol, coor,
))
obs_colorrange_164_shift = Observable(calc_prfl_colorrange_auto(
    obs_prfl_164_shift[]', t_evol, coor,
))
heatmap!(ax_prfl_162_shift, t_evol, coor, obs_prfl_162_shift;
    colormap=gen_clrmap_solo(hue_theme_istp["162"]),
    colorrange=obs_colorrange_162_shift,
)
heatmap!(ax_prfl_164_shift, t_evol, coor, obs_prfl_164_shift;
    colormap=gen_clrmap_solo(hue_theme_istp["164"]),
    colorrange=obs_colorrange_164_shift,
)
text!(ax_prfl_162_shift, lift(offsets -> "xₒ: $offsets", obs_x_off_162);
    position=(0.98, 0.98), space=:relative,
    align=(:right, :top), justification=:right, fontsize=14,
)
text!(ax_prfl_164_shift, lift(offsets -> "xₒ: $offsets", obs_x_off_164);
    position=(0.98, 0.98), space=:relative,
    align=(:right, :top), justification=:right, fontsize=14,
)
linkxaxes!([ax_prfl_162_orig, ax_prfl_162_shift, ax_prfl_164_orig, ax_prfl_164_shift])
xlims!(ax_prfl_164_shift, (0, 200))

ctrl = GridLayout(fig[5, 1])
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
prev_profile = Button(ctrl[1, 5]; label="←")
Label(ctrl[1, 6], lift(i -> "$(profile_modes[i])\n$(i)/$(length(profile_modes))", obs_profile_mode))
next_profile = Button(ctrl[1, 7]; label="→")
on(prev_profile.clicks) do _
    obs_profile_mode[] = mod1(obs_profile_mode[] - 1, length(profile_modes))
end
on(next_profile.clicks) do _
    obs_profile_mode[] = mod1(obs_profile_mode[] + 1, length(profile_modes))
end

Label(ctrl[1, 8], lift(r -> "blur $(round(r; digits=2)) μm", obs_r_blur_um))
slider_blur = Slider(ctrl[1, 9:12]; range=0.0:0.2:10.0, startvalue=fit_initial_blur, tellwidth=false)
Label(ctrl[1, 13], lift(r -> "fit $(r[1])–$(r[2]) ms", obs_rng_t_fit))
slider_fit_range = IntervalSlider(ctrl[1, 14:19]; range=0.0:1.0:200.0, startvalues=fit_initial_range)
Label(ctrl[1, 20], lift(t0 -> "initial t₀ $(round(t0; digits=1)) ms", obs_t0_initial))
slider_t0_initial = Slider(ctrl[1, 21:24]; range=-15.0:0.5:15.0, startvalue=fit_initial_t0, tellwidth=false)

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
    profile_mode = profile_modes[obs_profile_mode[]]
    _, prfl_evol_reps, coor, _ = profile_data(profile_mode)
    fit_now = make_fits(obs_ib[], profile_mode, obs_rng_t_fit[], obs_r_blur_um[], obs_t0_initial[])
    fit_live[] = fit_now
    fit_162, fit_164 = fit_now
    obs_prfl_162_orig[] = fit_162.evol_prfl_blur_full'
    obs_prfl_164_orig[] = fit_164.evol_prfl_blur_full'
    obs_colorrange_162_orig[] = calc_prfl_colorrange_auto(
        obs_prfl_162_orig[]', t_evol, coor,
    )
    obs_colorrange_164_orig[] = calc_prfl_colorrange_auto(
        obs_prfl_164_orig[]', t_evol, coor,
    )
    traj_plot_162[1][] = Point2f.(fit_162.t_evol_sel, fit_162.trajectory)
    traj_plot_164[1][] = Point2f.(fit_164.t_evol_sel, fit_164.trajectory)
    obs_fit_caption_162[] = format_fit_caption(fit_162)
    obs_fit_caption_164[] = format_fit_caption(fit_164)
    fit_shifted_162 = fit_shifted_profile_stack(
        prfl_evol_reps, obs_ib[], 1, fit_162.trajectory_full, t_evol, coor,
    )
    fit_shifted_164 = fit_shifted_profile_stack(
        prfl_evol_reps, obs_ib[], 2, fit_164.trajectory_full, t_evol, coor,
    )
    obs_prfl_162_shift[] = fit_shifted_162.evolution'
    obs_prfl_164_shift[] = fit_shifted_164.evolution'
    obs_x_off_162[] = format_x_off(fit_shifted_162.x_off)
    obs_x_off_164[] = format_x_off(fit_shifted_164.x_off)
    obs_colorrange_162_shift[] = calc_prfl_colorrange_auto(
        obs_prfl_162_shift[]', t_evol, coor,
    )
    obs_colorrange_164_shift[] = calc_prfl_colorrange_auto(
        obs_prfl_164_shift[]', t_evol, coor,
    )
end
onany(obs_ib, obs_profile_mode, obs_r_blur_um, obs_rng_t_fit, obs_t0_initial) do _...
    update_plot!()
end
update_plot!()

rowsize!(fig.layout, 1, 96)
rowsize!(fig.layout, 2, 96)
rowsize!(fig.layout, 3, 96)
rowsize!(fig.layout, 4, 96)
fig |> resize_to_layout!
fig |> display
