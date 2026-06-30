using CairoMakie: Figure, Axis, heatmap!, save
using GLMakie
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Demo"
path_output = joinpath(path_demo, "15.DualSS")
isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
step_grid = 0.25 / 50;
smwh_roi = (20, 80) .* 50

x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_grid

# assume all changes surrounds the center
function gen_dens(; λ_crys=3, σx=3, σy=8, σx_tf=3, σy_tf=8, A_tf=6, A_halo=4, x0=0, y0=0, φ=0, η=0.6, γ=1, soft=0.1)
    return (x, y) -> begin
        x, y = (x .- x0, y .- y0)
        nvlp_tf = @. A_tf * ((1 - soft) * clamp(1 - (x / σx_tf)^2 - (y / σy_tf)^2, 0, 1) + soft * exp(-((x / σx_tf)^2 + (y / σy_tf)^2)))
        nvlp_halo = @. A_halo * exp(-((x / σx)^2 + (y / σy)^2))
        crystal = @. ((1 + η * cos(2 * π * y / λ_crys - φ)) / 2)^γ * sqrt(γ)
        @. nvlp_tf * crystal + nvlp_halo
    end
end

function set_axis_dual!()
    fig = Figure(; backgroundcolor = :transparent)
    ax_1 = Axis(fig[1, 1]; width=400, height=100, backgroundcolor = :transparent)
    ax_2 = Axis(fig[2, 1]; width=400, height=100, backgroundcolor = :transparent)
    axs = [ax_1, ax_2]
    rowgap!(fig.layout, 0)
    fig, axs
end

function plot_duet!(axs, dens; max_dens=10)
    dens_grid = [dens(x, y) for x in x_posi, y in y_posi]

    for i in 1:2
        clrmap_dens = gen_clrmap_solo(hue_theme_istp[i == 1 ? "162" : "164"]; thres_alpha=0.1, alpha_base=-0.1)
        heatmap!(axs[i], y_posi, x_posi, dens_grid'; colorrange=(0, max_dens), colormap=clrmap_dens, rasterize=true)
        axs[i] |> hidedecorations!
        axs[i].aspect = DataAspect()
        axs[i] |> hidespines!
    end
end
fig_duet, axs_duet = set_axis_dual!()
for (desc, dens) in [
    ("[uniform bec     ]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=1, A_halo=4, x0=0, y0=0, φ=0π, η=0.0, γ=1.0, soft=0.3)),
    ("[supersolid      ]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=1, A_halo=4, x0=0, y0=0, φ=0π, η=1.0, γ=1.5, soft=0.3)),
    ("[isolated droplet]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=1, A_halo=4, x0=0, y0=0, φ=0π, η=1.0, γ=2.5, soft=0.3)),
    ("[roton instable  ]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=0.3, A_halo=5, x0=0, y0=0, φ=0π, η=1.0, γ=1.2, soft=0.5)),
]
    axs_duet |> clear_axes!
    plot_duet!(axs_duet, dens; max_dens=8)
    fig_duet |> resize_to_layout!
    fig_duet |> display
    save(joinpath(path_output, "dualss_demo_$desc.png"), fig_duet; px_per_unit = 1, backend=CairoMakie)
    save(joinpath(path_output, "dualss_demo_$desc.svg"), fig_duet; backend=CairoMakie)
    println("$desc displayed and saved")
end