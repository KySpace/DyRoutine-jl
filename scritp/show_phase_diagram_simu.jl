using CairoMakie: Figure, Axis, heatmap!, save
using GLMakie
using CSV
using DataFrames
using LaTeXStrings
using Pipe: @pipe
using Match: @match
using ScatteredInterpolation
using Statistics

GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path_simu = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Samples\[07.01].Weijing\phase diagram"
# commit #c018bbf9368558cbb09a629dcdd8a39cda93bbeb
path_output = joinpath(path_simu, "19.DualSS.PhaseDiagram.CW")
isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)

df_contrast_fine = CSV.read(joinpath(path_simu, "W_new_C.txt"), DataFrame; delim='\t', header=[:a12, :a22, :contrast_162, :contrast_164], skipto=2)
df_weightsp_fine = CSV.read(joinpath(path_simu, "W_new2.txt"), DataFrame; delim='\t', header=[:a12, :a22, :weightsp_162, :weightsp_164], skipto=2)
df_contrast_coarse = CSV.read(joinpath(path_simu, "W_coarse_C.txt"), DataFrame; delim='\t', header=[:a12, :a22, :contrast_162, :contrast_164], skipto=2)
df_weightsp_coarse = CSV.read(joinpath(path_simu, "W_coarse2.txt"), DataFrame; delim='\t', header=[:a12, :a22, :weightsp_162, :weightsp_164], skipto=2)
a12_ctrs_cat = vcat(df_contrast_fine.a12, df_contrast_coarse.a12)
a22_ctrs_cat = vcat(df_contrast_fine.a22, df_contrast_coarse.a22)
a12_wght_cat = vcat(df_weightsp_fine.a12, df_weightsp_coarse.a12)
a22_wght_cat = vcat(df_weightsp_fine.a22, df_weightsp_coarse.a22)
df_contrast_162_cat = vcat(df_contrast_fine.contrast_162, df_contrast_coarse.contrast_162)
df_contrast_164_cat = vcat(df_contrast_fine.contrast_164, df_contrast_coarse.contrast_164)
df_weightsp_162_cat = vcat(df_weightsp_fine.weightsp_162, df_weightsp_coarse.weightsp_162)
df_weightsp_164_cat = vcat(df_weightsp_fine.weightsp_164, df_weightsp_coarse.weightsp_164)

using Statistics

function average_duplicate_points(x, y, z1, z2) 
    d = Dict{Tuple{Float64, Float64}, Vector{Tuple{Float64, Float64}}}()
    for (xi, yi, z1i, z2i) in zip(vec(x), vec(y), vec(z1), vec(z2))
        if isfinite(xi) && isfinite(yi) && isfinite(z1i) && isfinite(z2i)
            push!(
                get!(d, (Float64(xi), Float64(yi)), Tuple{Float64, Float64}[]),
                (Float64(z1i), Float64(z2i)),
            )
        end
    end
    xs, ys, z1s, z2s = Float64[], Float64[], Float64[], Float64[]
    for ((xi, yi), vals) in d
        push!(xs, xi)
        push!(ys, yi)
        push!(z1s, mean(first.(vals)))
        push!(z2s, mean(last.(vals)))
    end
    return xs, ys, z1s, z2s
end

a22_ctrs, a12_ctrs, df_contrast_162, df_contrast_164 = average_duplicate_points(a22_ctrs_cat, a12_ctrs_cat, df_contrast_162_cat, df_contrast_164_cat)
a22_wght, a12_wght, df_weightsp_162, df_weightsp_164 = average_duplicate_points(a22_wght_cat, a12_wght_cat, df_weightsp_162_cat, df_weightsp_164_cat)

coor_ctrs = hcat(a22_ctrs, a12_ctrs)'
coor_wght = hcat(a22_wght, a12_wght)'

ntpl_contrast_162 = interpolate(Shepard(), coor_ctrs, df_contrast_162)
ntpl_contrast_164 = interpolate(Shepard(), coor_ctrs, df_contrast_164)
ntpl_weightsp_162 = interpolate(Shepard(), coor_wght, df_weightsp_162)
ntpl_weightsp_164 = interpolate(Shepard(), coor_wght, df_weightsp_164)

a22_q = range(88, 108; length=200)
a12_q = range(70,  96; length=200)
coor_q = reduce(hcat, ([xi, yi] for xi in a22_q, yi in a12_q))

contrast_162_q = @pipe evaluate(ntpl_contrast_162, coor_q) |> reshape(_, length(a22_q), length(a12_q))
contrast_164_q = @pipe evaluate(ntpl_contrast_164, coor_q) |> reshape(_, length(a22_q), length(a12_q))
weightsp_162_q = @pipe evaluate(ntpl_weightsp_162, coor_q) |> reshape(_, length(a22_q), length(a12_q))
weightsp_164_q = @pipe evaluate(ntpl_weightsp_164, coor_q) |> reshape(_, length(a22_q), length(a12_q))

fig = Figure();
Label(fig[0, 1]; text="contrast", valign=:center, halign=:center, font=:bold)
Label(fig[0, 2]; text="side peak weight", valign=:center, halign=:center, font=:bold)
ax_contrast_162 = Axis(fig[1, 1]; ylabel=L"a_{12} (a_0)", xlabel=L"a_{22} (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());
ax_weightsp_162 = Axis(fig[1, 2]; ylabel=L"a_{12} (a_0)", xlabel=L"a_{22} (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());
ax_contrast_164 = Axis(fig[2, 1]; ylabel=L"a_{12} (a_0)", xlabel=L"a_{22} (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());
ax_weightsp_164 = Axis(fig[2, 2]; ylabel=L"a_{12} (a_0)", xlabel=L"a_{22} (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());

clrmp_162 = gen_clrmap_solo(hue_theme_istp["162"])
clrmp_164 = gen_clrmap_solo(hue_theme_istp["164"])

function gen_clrfn(istp; thres_alpha=0.0, alpha_base=1.0)
    hue = hue_theme_istp[istp]
    clrfn= t ->
        begin
            alpha = thres_alpha <= 0 || abs(t) > thres_alpha ?
                    1.0 :
                    clamp(abs(t) / thres_alpha * (1 - alpha_base) + alpha_base, 0, 1)
            Oklch(1 - 0.8 * t, 0.24 * t, hue) |> c -> RGBAf(c, alpha)
        end
    return clrfn
end

clrrng_c = extrema(vcat(vec(contrast_162_q), vec(contrast_164_q)))
clrrng_w = extrema(vcat(vec(weightsp_162_q), vec(weightsp_164_q)))

hm_c1 = heatmap!(ax_contrast_162, a22_q, a12_q, contrast_162_q; colormap=:BrBg, colorrange=clrrng_c)
hm_w1 = heatmap!(ax_weightsp_162, a22_q, a12_q, weightsp_162_q; colormap=:BrBg, colorrange=clrrng_w)
hm_c2 = heatmap!(ax_contrast_164, a22_q, a12_q, contrast_164_q; colormap=:BrBg, colorrange=clrrng_c)
hm_w2 = heatmap!(ax_weightsp_164, a22_q, a12_q, weightsp_164_q; colormap=:BrBg, colorrange=clrrng_w)

Colorbar(fig[3, 1], hm_c1; vertical=false);
Colorbar(fig[4, 1], hm_c2; vertical=false);
Colorbar(fig[3, 2], hm_w1; vertical=false);
Colorbar(fig[4, 2], hm_w2; vertical=false);

colsize!(fig.layout, 1, 500)
colsize!(fig.layout, 2, 500)
rowsize!(fig.layout, 1, 300)
rowsize!(fig.layout, 2, 300)

fig |> resize_to_layout!
fig |> display
# fig |> f -> save(joinpath(path_output, "phase_diagram_cw.png"), f; px_per_unit=2.0, backend=CairoMakie)
# fig |> f -> save(joinpath(path_output, "phase_diagram_cw.svg"), f; px_per_unit=2.0, backend=CairoMakie)
# fig |> f -> save(joinpath(path_output, "phase_diagram_cw.pdf"), f; px_per_unit=2.0, backend=CairoMakie)
