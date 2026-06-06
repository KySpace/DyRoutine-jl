include(joinpath(@__DIR__, "..", "src", "vispca.jl"))

step_grid = 0.25;
smwh_roi = (40, 80)

x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi

# assume all changes surrounds the center
function gen_dens(λ_crys, Δ_sf, W_sf)
    σx, σy = (4, 10)
    A = 10
    η = 0.8
    return (x, y) -> begin
        nvlp = @. A * exp(-((x/σx)^2 + (y/σy)^2))
        crystal = @. 1 + η * cos(2 * π * x / λ_crys) / 2
        sf = @. 1 + W_sf * cos(2 * π * x / (λ_crys * Δ_sf)) / 2
        @. nvlp * crystal * sf
    end
end

dens_ground = [gen_dens(3, 1, 0)(x, y) for x in x_vec, y in y_vec]
dens_crys = [gen_dens(3, 4, 0)(x, y) for x in x_vec, y in y_vec]
dens_sf = [gen_dens(3.2, 1, 0)(x, y) for x in x_vec, y in y_vec]

function plot_mode(axs, dens_ground, dens_excit)
    heatmap!(axs["dens2d-ground"], x_vec, y_vec, dens_ground; colormap = :viridis)
    heatmap!(axs["dens2d-excit"], x_vec, y_vec, dens_excit; colormap = :viridis)

end
