include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))

path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Demo"
path_output = joinpath(path_demo, "11.SFMode.[3→3.2]")
isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
step_grid = 0.25 / 10;
smwh_roi = (20, 80) .* 10

x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_grid

# assume all changes surrounds the center
function gen_dens(λ_crys, Δ_sf, W_sf)
    σx, σy = (2, 10)
    A = 10
    η = 0.6
    δ = 0.1
    return (x, y) -> begin
        nvlp = @. A * exp(-((x / σx)^2 + (y / σy)^2))
        crystal = @. (1 + η * cos(2 * π * y / λ_crys)) / 2
        sf = @. (1 + W_sf / atan(1 / δ) * atan(cos(π * y / (λ_crys * Δ_sf)) / δ)) / 2
        @. nvlp * crystal * sf
    end
end

dens_ground = [gen_dens(3, 1, 0)(x, y) for x in x_posi, y in y_posi]
dens_crys = [gen_dens(3, 1, 0.25)(x, y) for x in x_posi, y in y_posi]
dens_sf = [gen_dens(3.2, 1, 0)(x, y) for x in x_posi, y in y_posi]

fig_modes, axs_modes =
    begin
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

# plot_mode(axs_modes, dens_ground, dens_crys)
plot_mode(axs_modes, dens_ground, dens_sf)
fig_modes |> resize_to_layout!
fig_modes |> display
save(joinpath(path_output, "test_pca_modes.png"), fig_modes)
