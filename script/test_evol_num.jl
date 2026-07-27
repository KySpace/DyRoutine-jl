using GLMakie
using JLD2
using Colors
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
using GeometryBasics: Point2f


path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
title_load = "[07.24].100.Extr"
cache_extr = JLD2.load(joinpath(path_root, "AnlzRoutine", title_load, "NTRC_essn_extr.jld2"))
prfl_axial_evol = JLD2.load(joinpath(path_root, "AnlzRoutine", title_load, "NTRC_corr.jld2"))["prfl_axial_evol"]
meta_extr = cache_extr["meta_extr"]
val_vars = meta_extr.val_vars
px_in_um = meta_extr.px_in_um
coor = (-60:1:60)*px_in_um
max_err_rel = 1.6

num_fit = map(cache_extr["extr_fmt"]) do extr
    p = extr.envelope.params_asymm
    2π * prod(p.size) * p.max / (px_in_um^2)
end


function twostep_decay(t::AbstractVector, p::AbstractVector)
    D1, λ1, D2, λ2 = p
    @. D1 * exp(- t / λ1) + D2 * exp(- t / λ2)
end
function fit_num_decay(t::AbstractVector{<:Real}, nums::AbstractVector{<:Real})
    ymax = maximum(nums)
    ymax > 0 || throw(ArgumentError("number values must contain a positive maximum"))
    p0 = [ymax / 2, 50.0, ymax / 2, 200.0]
    lower = [0.0, 1e-6, 0.0, 1e-6]
    upper = [ymax, Inf, ymax, Inf]
    fit = curve_fit(twostep_decay, Float64.(t), Float64.(nums), p0; lower, upper)
    return coef(fit)
end

rng_t_num_decay = (0.0, 150.0)
t_hold = Float64.(val_vars.t_hold)
mask_t = (rng_t_num_decay[1] .<= t_hold .<= rng_t_num_decay[2])

function to_num_heat_clr(value::Real, valid::Bool, hue; max_value=5, thres_alpha=0.1, alpha_base=0.1)
    value_norm = clamp(value, 0, max_value) / max_value
    alpha = thres_alpha <= 0 || value_norm > thres_alpha ? 1.0 :
        clamp(value_norm / thres_alpha * (1 - alpha_base) + alpha_base, 0, 1)
    chroma = 0.24 * value_norm
    return RGBAf(Oklch(1 - 0.8 * value_norm, chroma, hue + (valid ? 0 : 60)), alpha)
end

function make_num_heat_clr(profile::AbstractMatrix, valid::AbstractVector, hue)
    size(profile, 1) == length(valid) || throw(DimensionMismatch("profile time dimension and validity mask differ"))
    return [to_num_heat_clr(profile[i, j], valid[i], hue) for i in axes(profile, 1), j in axes(profile, 2)]
end

function make_num_fit(ib::Int, istp::Int)
    t_fit = t_hold[mask_t]
    nums_full = [Float64.(vec(num_fit[ib, r, :, istp])) for r in 1:3]
    nums_fit = [nums_full[r][mask_t] for r in 1:3]
    params = fit_num_decay(repeat(t_fit, 3), vcat(nums_fit...))
    fitted_full = [twostep_decay(t_hold, params) for _ in 1:3]
    errors_rel_full = [(nums_full[r] .- fitted_full[r]) ./ fitted_full[r] for r in 1:3]
    σ = std(vcat([errors_rel_full[r][mask_t] for r in 1:3]...))
    errors_full = [nums_full[r] .- fitted_full[r] for r in 1:3]
    valid_num_local = reduce(vcat, (permutedims(abs.(errors_rel) .<= max_err_rel * σ) for errors_rel in errors_rel_full))
    return (; t=t_hold, nums=nums_full, params, fitted=fitted_full, errors=errors_full, errors_rel=errors_rel_full, σ, valid_num=valid_num_local)
end

ib0 = findfirst(==(5.316), val_vars.IB)
ib0 = something(ib0, 1)
obs_ib = Observable(ib0)
obs_rep = Observable(2)
obs_istp = Observable(2)

fits_num = [make_num_fit(ib, istp) for ib in axes(num_fit, 1), istp in axes(num_fit, 4)]
valid_num = permutedims(
    reduce((a, b) -> cat(a, b; dims=4),
        [cat([fits_num[ib, istp].valid_num for istp in axes(num_fit, 4)]...; dims=3)
         for ib in axes(num_fit, 1)]),
    (4, 1, 2, 3),
)

fig_live = Figure(size=(1500, 500))
axs_prfl = Axis(fig_live[1, 1], ylabel="position (μm)")
axs_num = Axis(fig_live[2, 1], ylabel="number")
axs_err = Axis(fig_live[3, 1], xlabel="t_hold (ms)", ylabel="error")
hue_rep = [279, 127, 85, 52, 318]
fit_live = fits_num[obs_ib[], obs_istp[]]
profile_initial = prfl_axial_evol[obs_ib[], obs_rep[], obs_istp[]]'
valid_initial = valid_num[obs_ib[], obs_rep[], :, obs_istp[]]
obs_profile_clr = Observable(make_num_heat_clr(profile_initial, valid_initial, hue_theme_istp[string(val_vars.istp[obs_istp[]])]))
hm_profile = heatmap!(axs_prfl, t_hold, coor, obs_profile_clr)
fit_lines = [lines!(axs_num, fit_live.t, fit_live.fitted[r]; color=:black, linewidth=2) for r in 1:3]
scatter_plots = [scatter!(axs_num, fit_live.t, fit_live.nums[r]; color=RGBAf(Oklch(0.6331, 0.0923, hue_rep[r]), r == obs_rep[] ? 1.0 : 0.35)) for r in 1:3]
err_plots = [scatter!(axs_err, fit_live.t, fit_live.errors[r]; color=RGBAf(Oklch(0.6331, 0.0923, hue_rep[r]), r == obs_rep[] ? 1.0 : 0.35)) for r in 1:3]
band_low = Observable(-max_err_rel * fit_live.σ .* fit_live.fitted[obs_rep[]])
band_high = Observable(max_err_rel * fit_live.σ .* fit_live.fitted[obs_rep[]])
band_err = band!(axs_err, fit_live.t, band_low, band_high; color=(:gray, 0.25))
fit_caption = Observable("")
text!(axs_num, fit_caption; position=Point2f(0.98, 0.98), space=:relative, align=(:right, :top),
    justification=:right, fontsize=14)
linkxaxes!([axs_num, axs_prfl])
linkxaxes!([axs_num, axs_err])

for row in 1:3
    rowsize!(fig_live.layout, row, 140)
end

ctrl = GridLayout(fig_live[4, 1])
function add_cycle!(col, label, obs, values)
    prev = Button(ctrl[1, col]; label="←")
    value_label = if label == "IB"
        x -> "$(values[x]) A\n$(x)/$(length(values))"
    elseif label == "istp"
        x -> "$(values[x])\n$(x)/$(length(values))"
    else
        x -> "rep $(values[x])\n$(x)/$(length(values))"
    end
    Label(ctrl[1, col + 1], lift(value_label, obs))
    next = Button(ctrl[1, col + 2]; label="→")
    on(prev.clicks) do _; obs[] = mod1(obs[] - 1, length(values)); end
    on(next.clicks) do _; obs[] = mod1(obs[] + 1, length(values)); end
end
add_cycle!(1, "IB", obs_ib, val_vars.IB)
add_cycle!(4, "rep", obs_rep, 1:3)
add_cycle!(7, "istp", obs_istp, val_vars.istp)

function update_plot!()
    fit_now = fits_num[obs_ib[], obs_istp[]]
    for r in 1:3
        scatter_plots[r].color[] = RGBAf(Oklch(0.6331, 0.0923, hue_rep[r]), r == obs_rep[] ? 1.0 : 0.35)
        err_plots[r].color[] = scatter_plots[r].color[]
        scatter_plots[r][1][] = Point2f.(fit_now.t, fit_now.nums[r])
        err_plots[r][1][] = Point2f.(fit_now.t, fit_now.errors[r])
        fit_lines[r][1][] = Point2f.(fit_now.t, fit_now.fitted[r])
    end
    valid_now = valid_num[obs_ib[], obs_rep[], :, obs_istp[]]
    profile_now = prfl_axial_evol[obs_ib[], obs_rep[], obs_istp[]]'
    obs_profile_clr[] = make_num_heat_clr(profile_now, valid_now, hue_theme_istp[string(val_vars.istp[obs_istp[]])])
    band_low[] = -max_err_rel * fit_now.σ .* fit_now.fitted[obs_rep[]]
    band_high[] = max_err_rel * fit_now.σ .* fit_now.fitted[obs_rep[]]
    D1, λ1, D2, λ2 = fit_now.params
    fit_caption[] = "D₁=$(@sprintf("%.2g", D1)), λ₁=$(@sprintf("%.2g", λ1))\n" *
        "D₂=$(@sprintf("%.2g", D2)), λ₂=$(@sprintf("%.2g", λ2))\n" *
        "σᵣ=$(round(fit_now.σ; digits=4))"
end
onany(obs_ib, obs_rep, obs_istp) do _...; update_plot!(); end
update_plot!()
fig_live |> resize_to_layout!
fig_live |> display
