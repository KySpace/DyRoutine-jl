using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
using ImageFiltering
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))

path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Demo"
path_output = joinpath(path_demo, "12.MultiModes")
isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
step_grid = 0.25 / 10;
smwh_roi = (20, 80) .* 10

x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_grid

# assume all changes surrounds the center
function gen_dens(λ_crys, Δ_sf, W_sf; η=0.6, σx=2, σy=10, A=10, δ=0.1, x0=0, y0=0)
    return (x, y) -> begin
        nvlp = @. A * 2 * 10 / (σx * σy) * exp(-(((x - x0) / σx)^2 + ((y - y0) / σy)^2))
        crystal = @. (1 + η * cos(2 * π * (y - y0) / λ_crys)) / 2
        sf = @. (1 + W_sf / atan(1 / δ) * atan(cos(π * (y - y0) / (λ_crys * Δ_sf)) / δ)) / 2
        @. nvlp * crystal * sf
    end
end

dens_ground = [gen_dens(3, 1, 0)(x, y) for x in x_posi, y in y_posi]
dens_crys = [gen_dens(3, 1, 0.25)(x, y) for x in x_posi, y in y_posi]
dens_crys_wide = [gen_dens(3, 3, 0.25)(x, y) for x in x_posi, y in y_posi]
dens_sf = [gen_dens(3.2, 1, 0)(x, y) for x in x_posi, y in y_posi]
dens_higgs = [gen_dens(3, 1, 0; η=0.8)(x, y) for x in x_posi, y in y_posi]
dens_breath_x = [gen_dens(3, 1, 0; σx=2.5)(x, y) for x in x_posi, y in y_posi]
dens_breath_y = [gen_dens(3, 1, 0; σy=12)(x, y) for x in x_posi, y in y_posi]
dens_quad = [gen_dens(3, 1, 0; σx=1.5, σy=13)(x, y) for x in x_posi, y in y_posi]
dens_dipl_x = [gen_dens(3, 1, 0; x0=0.5)(x, y) for x in x_posi, y in y_posi]
dens_dipl_y = [gen_dens(3, 1, 0; y0=4.0)(x, y) for x in x_posi, y in y_posi]
dens_dipl_x_crys = [gen_dens(3, 1, 0.25; x0=0.5)(x, y) for x in x_posi, y in y_posi]
dens_dipl_x_sf = [gen_dens(3.2, 1, 0; x0=0.5)(x, y) for x in x_posi, y in y_posi]
dens_dipl_x_higgs = [gen_dens(3, 1, 0; η=0.8, x0=0.5)(x, y) for x in x_posi, y in y_posi]
dens_dipl_y_sf = [gen_dens(3, 1, 0; η=0.8, y0=4.0)(x, y) for x in x_posi, y in y_posi]
dens_breath_xy = [gen_dens(3, 1, 0; σx=2.5, σy=12)(x, y) for x in x_posi, y in y_posi]
function set_axis_compare_mode()
    fig = Figure()
    ax_ground_dens = Axis(fig[1, 1]; width=400, height=100)
    ax_excit_dens = Axis(fig[2, 1]; width=400, height=100)
    ax_mode_dens = Axis(fig[3, 1]; width=400, height=100)
    ax_1d_profile = Axis(fig[4, 1]; width=400, height=100)
    axs = Dict("dens2d-ground" => ax_ground_dens, "dens2d-excit" => ax_excit_dens, "dens2d-mode" => ax_mode_dens, "1d-profile" => ax_1d_profile)
    fig, axs
end

function plot_mode(axs, dens_ground, dens_excit)
    clrmap_mono = gen_clrmap_solo(280)
    clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360)
    mode_dens = dens_excit - dens_ground
    max_mode = maximum(abs, mode_dens)
    max_dens = maximum(abs, stack((dens_ground, dens_excit)))
    heatmap!(axs["dens2d-ground"], y_posi, x_posi, dens_ground'; colormap=clrmap_mono, colorrange=(0, max_dens))
    heatmap!(axs["dens2d-excit"], y_posi, x_posi, dens_excit'; colormap=clrmap_mono, colorrange=(0, max_dens))
    heatmap!(axs["dens2d-mode"], y_posi, x_posi, mode_dens'; colormap=clrmap, colorrange=(-max_mode, max_mode))
    for ax in matching_axes(axs, r"dens2d")
        ax.aspect = DataAspect()
    end
    lines!(axs["1d-profile"], y_posi, dens_ground[smwh_roi[1], :] |> vec; color=:black, linestyle=:dash, label="ground")
    lines!(axs["1d-profile"], y_posi, dens_excit[smwh_roi[1], :] |> vec; color=:black, label="excit")
    xlims!(axs["1d-profile"], (y_posi[1], y_posi[end]))
end

for (name, dens) in [
    ("crys", dens_crys),
    ("crys_wide", dens_crys_wide),
    ("sf", dens_sf),
    ("higgs", dens_higgs),
    ("breath_x", dens_breath_x),
    ("breath_y", dens_breath_y),
    ("quad", dens_quad),
    ("dipl_x", dens_dipl_x),
    ("dipl_y", dens_dipl_y),
    ("dipl_x_crys", dens_dipl_x_crys),
    ("dipl_x_sf", dens_dipl_x_sf),
    ("dipl_x_higgs", dens_dipl_x_higgs),
    ("dipl_y_sf", dens_dipl_y_sf),
    ("breath_xy", dens_breath_xy)
    ]
    fig_modes, axs_modes = set_axis_compare_mode()
    plot_mode(axs_modes, dens_ground, dens)
    fig_modes |> resize_to_layout!
    fig_modes |> display
    save(joinpath(path_output, "test_$(name)_modes.png"), fig_modes)
end
