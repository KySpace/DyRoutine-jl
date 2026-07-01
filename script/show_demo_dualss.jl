using CairoMakie: Figure, Axis, heatmap!, save
using GLMakie
using CSV
using DataFrames
using Pipe: @pipe
using Match: @match
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Demo"
path_simu = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Samples\[07.01].Weijing\density profile\a12=78"
# commit #c018bbf9368558cbb09a629dcdd8a39cda93bbeb
path_output = joinpath(path_demo, "17.DualSS.ManySimu")
isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
step_grid = 0.25 / 10;
smwh_roi = (20, 80) .* 10

x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_grid

function read_dens_simu(dir_dens, names_dens)
    df1 = CSV.read(joinpath(dir_dens, names_dens[1]), DataFrame; delim=' ', header=[:x, :y, :dens, :img], skipto=1)
    df2 = CSV.read(joinpath(dir_dens, names_dens[2]), DataFrame; delim=' ', header=[:x, :y, :dens, :img], skipto=1)
    local x_posi = unique(df1.x)
    local y_posi = unique(df1.y)
    hw = length.((y_posi, x_posi))
    dens = [reshape(df1.dens, hw), reshape(df2.dens, hw)]
    unit_in_um = 1.43
    (; x_posi=y_posi .* unit_in_um, y_posi=x_posi .* unit_in_um, dens)
end

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
    fig = Figure(; backgroundcolor=:transparent)
    ax_1 = Axis(fig[1, 1]; width=400, height=100, backgroundcolor=:transparent)
    ax_2 = Axis(fig[2, 1]; width=400, height=100, backgroundcolor=:transparent)
    axs = [ax_1, ax_2]
    rowgap!(fig.layout, 0)
    fig, axs
end

function plot_duet!(axs, xydens; max_dens=10)
    local x_posi, y_posi, dens = xydens

    for i in 1:2
        clrmap_dens = gen_clrmap_solo(hue_theme_istp[i == 1 ? "162" : "164"]; thres_alpha=0.1, alpha_base=-0.1)
        heatmap!(axs[i], y_posi, x_posi, dens[i]'; colorrange=(0, max_dens), colormap=clrmap_dens, rasterize=true)
        axs[i] |> hidedecorations!
        axs[i].aspect = DataAspect()
        axs[i] |> hidespines!
    end
end

fig_duet, axs_duet = set_axis_dual!()
fmt_demo = gen -> (; x_posi, y_posi, dens=([gen(x, y) for x in x_posi, y in y_posi] |> a -> [a, a]))
gen_desc_dens_simu = (a22, proj) -> begin
    prefix = @match proj begin
            "xy" => "density%d"
            "xz" => "density%dz"
        end
    filename_fmt = Printf.Format(prefix * "_78_%g")
    @sprintf("[simu %s a12=78'a22=%.04f]", proj, a22),
    read_dens_simu(joinpath(path_simu, "outdens_$proj"), [Printf.format(filename_fmt, i, a22) for i in 1:2])
    end
a22_sample_xy = [
    88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 98.2, 98.4, 98.6, 98.8, 99, 99.2, 99.4, 99.6,
    99.8000,
    99.8105,
    99.8211,
    99.8316,
    99.8421,
    99.8526,
    99.8632,
    99.8737,
    99.8842,
    99.8947,
    99.9053,
    99.9158,
    99.9263,
    99.9368,
    99.9474,
    99.9579,
    99.9684,
    99.9789,
    99.9895,
    100.0000,
    101.0000,
    102.0000,
    103.0000,
    104.0000,
    105.0000,
    106.0000,
    107.0000,
    108.0000,
]
a22_sample_xz = [

]
for (desc, dens) in [
    ("[demo uniform bec     ]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=1, A_halo=4, x0=0, y0=0, φ=0π, η=0.0, γ=1.0, soft=0.3) |> fmt_demo),
    ("[demo supersolid      ]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=1, A_halo=4, x0=0, y0=0, φ=0π, η=1.0, γ=1.5, soft=0.3) |> fmt_demo),
    ("[demo isolated droplet]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=1, A_halo=4, x0=0, y0=0, φ=0π, η=1.0, γ=2.5, soft=0.3) |> fmt_demo),
    ("[demo roton instable  ]", gen_dens(; λ_crys=2.5, σx=2.5, σy=10, σx_tf=2.5, σy_tf=10, A_tf=0.3, A_halo=5, x0=0, y0=0, φ=0π, η=1.0, γ=1.2, soft=0.5) |> fmt_demo),
    # ("[simu 1]", read_dens_simu(joinpath(path_simu), ["density1_1", "density2_1"])),
    # ("[simu 2]", read_dens_simu(joinpath(path_simu), ["density1_2", "density2_2"])),
    [gen_desc_dens_simu(a22, "xy") for a22 in a22_sample_xy]...,
    [gen_desc_dens_simu(a22, "xz") for a22 in a22_sample_xz]...,
]
    axs_duet |> clear_axes!
    max_dens, xlim, ylim = contains(desc, "simu") ?
                           (0.5, (-10, 10), (-2.5, 2.5)) :
                           (8., (-20, 20), (-5, 5))
    plot_duet!(axs_duet, dens; max_dens)
    fig_duet |> resize_to_layout!
    fig_duet |> display
    limits!(axs_duet[1], xlim, ylim)
    linkaxes!(axs_duet)
    save(joinpath(path_output, "$desc.png"), fig_duet; px_per_unit=1, backend=CairoMakie)
    save(joinpath(path_output, "$desc.svg"), fig_duet; backend=CairoMakie)
    println("$desc displayed and saved")
end
