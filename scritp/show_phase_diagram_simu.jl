using CairoMakie: Figure, Axis, heatmap!, save
using GLMakie
using CSV
using DataFrames
using LaTeXStrings
using Pipe: @pipe
using Match: @match
using NaturalNeighbours
using Statistics

GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

path_simu = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Samples\[07.01].Weijing\phase diagram"
path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Demo"
# commit #c018bbf9368558cbb09a629dcdd8a39cda93bbeb
path_output = joinpath(path_demo, "20.DualSS.PhaseDiagram.CW")
isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)

df_contrast_fine = CSV.read(joinpath(path_simu, "W_new_C.txt"), DataFrame; delim='\t', header=[:a12, :a22, :contrast_162, :contrast_164], skipto=2)
df_weightsp_fine = CSV.read(joinpath(path_simu, "W_new2.txt"), DataFrame; delim='\t', header=[:a12, :a22, :weightsp_162, :weightsp_164], skipto=2)
df_contrast_coarse = CSV.read(joinpath(path_simu, "W_coarse_C.txt"), DataFrame; delim='\t', header=[:a12, :a22, :contrast_162, :contrast_164], skipto=2)
df_weightsp_coarse = CSV.read(joinpath(path_simu, "W_coarse2.txt"), DataFrame; delim='\t', header=[:a12, :a22, :weightsp_162, :weightsp_164], skipto=2)
df_contrast_a12_78 = CSV.read(joinpath(path_simu, "contrast_combined.txt"), DataFrame; delim='\t', header=[:a12, :a22, :contrast_162, :contrast_164], skipto=2)
a12_ctrs_cat = vcat(df_contrast_fine.a12, df_contrast_coarse.a12)
a22_ctrs_cat = vcat(df_contrast_fine.a22, df_contrast_coarse.a22)
a12_wght_cat = vcat(df_weightsp_fine.a12, df_weightsp_coarse.a12)
a22_wght_cat = vcat(df_weightsp_fine.a22, df_weightsp_coarse.a22)
df_contrast_162_cat = vcat(df_contrast_fine.contrast_162, df_contrast_coarse.contrast_162)
df_contrast_164_cat = vcat(df_contrast_fine.contrast_164, df_contrast_coarse.contrast_164)
df_weightsp_162_cat = vcat(df_weightsp_fine.weightsp_162, df_weightsp_coarse.weightsp_162)
df_weightsp_164_cat = vcat(df_weightsp_fine.weightsp_164, df_weightsp_coarse.weightsp_164)


function average_duplicate_points(x, y, z1, z2)
    d = Dict{Tuple{Float64,Float64},Vector{Tuple{Float64,Float64}}}()
    for (xi, yi, z1i, z2i) in zip(vec(x), vec(y), vec(z1), vec(z2))
        if isfinite(xi) && isfinite(yi) && isfinite(z1i) && isfinite(z2i)
            push!(
                get!(d, (Float64(xi), Float64(yi)), Tuple{Float64,Float64}[]),
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

ntpl_contrast_162 = interpolate(a22_ctrs, a12_ctrs, df_contrast_162; derivatives=true)
ntpl_contrast_164 = interpolate(a22_ctrs, a12_ctrs, df_contrast_164; derivatives=true)
ntpl_weightsp_162 = interpolate(a22_wght, a12_wght, df_weightsp_162; derivatives=true)
ntpl_weightsp_164 = interpolate(a22_wght, a12_wght, df_weightsp_164; derivatives=true)

a22_g = range(88, 108; length=200)
a12_g = range(70, 96; length=200)
a22_q = vec([xi for xi in a22_g, yi in a12_g])
a12_q = vec([yi for xi in a22_g, yi in a12_g])

contrast_162_q = @pipe ntpl_contrast_162(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))
contrast_164_q = @pipe ntpl_contrast_164(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))
weightsp_162_q = @pipe ntpl_weightsp_162(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))
weightsp_164_q = @pipe ntpl_weightsp_164(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))

sample_contrast = [(96, 78, :utriangle), (99.8421, 78, :diamond), (104, 78, :circle)]
##
fig_full = Figure();
Label(fig_full[1, 0]; text=L"^{162}\text{Dy}", valign=:center, halign=:center, fontsize=16)
Label(fig_full[2, 0]; text=L"^{164}\text{Dy}", valign=:center, halign=:center, fontsize=16)
Label(fig_full[0, 1]; text="contrast", valign=:center, halign=:center, font=:bold)
Label(fig_full[0, 2]; text="side peak weight", valign=:center, halign=:center, font=:bold)
ax_contrast_162 = Axis(fig_full[1, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());
ax_weightsp_162 = Axis(fig_full[1, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());
ax_contrast_164 = Axis(fig_full[2, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());
ax_weightsp_164 = Axis(fig_full[2, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());

fig_ctrs_164 = Figure()
Label(fig_ctrs_164[1, 0]; text=L"^{164}\text{Dy}", valign=:center, halign=:center, fontsize=16)
ax_contrast_sample = Axis(fig_ctrs_164[1, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", xlabelsize=16, ylabelsize=16, aspect=DataAspect());

fig_a1278 = Figure()
Label(fig_a1278[0, 1]; text=L"a_{12} = 78 a_0", valign=:center, halign=:center, fontsize=16)
ax_a1278 = Axis(fig_a1278[1, 1]; ylabel=L"C", xlabel=L"a_{22} \; (a_0)", xlabelsize=16, ylabelsize=16, width=600, height=200);
Box(fig_a1278[1, 1]; color=:white, width=160, height=100, halign=0.14, valign=0.20)
Box(fig_a1278[1, 1]; color=(Oklch(0.90, 0.005, 192), 0.2), width=160, height=100, halign=0.14, valign=0.20)
ax_a1278_zoom = Axis(fig_a1278[1, 1]; backgroundcolor=:white, width=160, height=100, halign=0.14, valign=0.20, xticklabelsize=10, yticklabelsize=10, xgridvisible=false, ygridvisible=false);


function gen_clrmap_parabola(hue, light_maxchroma, chroma_max, light_min; thres_alpha=0.0, alpha_base=1.0, light_max=1.0, chroma_lightmax=0, hue_range=(0, 0))
    clrmap = [
        begin
            l = (1 - t) * (light_max - light_min) + light_min
            c = (chroma_lightmax - chroma_max) / (light_max - light_maxchroma)^2 * (l - light_maxchroma)^2 + chroma_max
            h = hue + (t - 0.5) * (hue_range[2] - hue_range[1])
            alpha = thres_alpha <= 0 || abs(t) > thres_alpha ?
                    1.0 :
                    clamp(abs(t) / thres_alpha * (1 - alpha_base) + alpha_base, 0, 1)
            Oklch(clamp(l, 0, 1), clamp(c, 0, 1), clamp(h, 0, 360)) |> c -> RGBAf(c, alpha)
        end
        for t in range(0, 1; length=256)]
    return clrmap
end

clrmp_162 = gen_clrmap_solo(hue_theme_istp["162"])
clrmp_164 = gen_clrmap_solo(hue_theme_istp["164"])
clrmp_turqoise = gen_clrmap_parabola(192, 0.58, 0.06, 0.55; hue_range=(0, 0), light_max=0.97, chroma_lightmax=0.008, thres_alpha=0.01, alpha_base=0.7)

function gen_clrfn(istp; thres_alpha=0.0, alpha_base=1.0)
    hue = hue_theme_istp[istp]
    clrfn = t ->
        begin
            alpha = thres_alpha <= 0 || abs(t) > thres_alpha ?
                    1.0 :
                    clamp(abs(t) / thres_alpha * (1 - alpha_base) + alpha_base, 0, 1)
            Oklch(1 - 0.8 * t, 0.24 * t, hue) |> c -> RGBAf(c, alpha)
        end
    return clrfn
end




clrrng_c = (0, 1)
clrrng_w = (0, 0.25) # extrema(vcat(vec(weightsp_162_q), vec(weightsp_164_q)))

hm_c1 = heatmap!(ax_contrast_162, a22_g, a12_g, contrast_162_q; colormap=clrmp_162, colorrange=clrrng_c)
hm_w1 = heatmap!(ax_weightsp_162, a22_g, a12_g, weightsp_162_q; colormap=clrmp_162, colorrange=clrrng_w)
hm_c2 = heatmap!(ax_contrast_164, a22_g, a12_g, contrast_164_q; colormap=clrmp_164, colorrange=clrrng_c)
hm_w2 = heatmap!(ax_weightsp_164, a22_g, a12_g, weightsp_164_q; colormap=clrmp_164, colorrange=clrrng_w)
hm_cs = heatmap!(ax_contrast_sample, a22_g, a12_g, contrast_164_q; colormap=clrmp_turqoise, colorrange=clrrng_c)
for (a22, a12, marker) in sample_contrast
    scatter!(ax_contrast_sample, [a22], [a12]; color=:mediumpurple4, marker=marker, markersize=8)
end

Colorbar(fig_full[3, 1], hm_c1; vertical=false, label=L"C");
Colorbar(fig_full[4, 1], hm_c2; vertical=false, label=L"C");
Colorbar(fig_full[3, 2], hm_w1; vertical=false, label=L"W");
Colorbar(fig_full[4, 2], hm_w2; vertical=false, label=L"W");
Colorbar(fig_ctrs_164[1, 2], hm_cs; vertical=true, label=L"\text{contrast} \; ^{164}\text{Dy}");
limits!(ax_contrast_sample, (95, 108), (70, 96))

colsize!(fig_full.layout, 0, 20)
colsize!(fig_full.layout, 1, 300)
colsize!(fig_full.layout, 2, 300)
colgap!(fig_full.layout, 10)
rowsize!(fig_full.layout, 1, 360)
rowsize!(fig_full.layout, 2, 360)

colsize!(fig_ctrs_164.layout, 1, 360)
rowsize!(fig_ctrs_164.layout, 1, 400)
xlims!(ax_contrast_sample, (90, 106))
ylims!(ax_contrast_sample, (70, 90))
ax_contrast_sample.xticks = 90:2:110
ax_contrast_sample.yticks = 70:2:106

xlims!(ax_a1278, (96, 101))
ylims!(ax_a1278, (-0.05, 0.85))
xlims!(ax_a1278_zoom, (99.8, 99.9))
ylims!(ax_a1278_zoom, (-0.02, 0.32))
vspan!(ax_a1278, 99.8, 99.9; color=(Oklch(0.90, 0.005, 192), 0.5))
cut_c_162 = scatterlines!(ax_a1278, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_162; color=Oklch(0.4, 0.14, hue_theme_istp["162"]), label=L"^{162}\text{Dy}")
cut_c_164 = scatterlines!(ax_a1278, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_164; color=Oklch(0.4, 0.14, hue_theme_istp["164"]), label=L"^{162}\text{Dy}")
scatterlines!(ax_a1278_zoom, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_162; color=Oklch(0.4, 0.14, hue_theme_istp["162"]))
scatterlines!(ax_a1278_zoom, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_164; color=Oklch(0.4, 0.14, hue_theme_istp["164"]))
axislegend(ax_a1278; position=:rt, framevisible=false, labelsize=14)
ax_a1278.xticks = 90:2:110
ax_a1278.yticks = 0:0.2:1

fig_full |> resize_to_layout!
fig_full |> display
fig_ctrs_164 |> resize_to_layout!
fig_ctrs_164 |> display
fig_a1278 |> resize_to_layout!
fig_a1278 |> display

##
fig_full |> f -> save(joinpath(path_output, "phase_diagram_W_sample.png"), f; px_per_unit=2.0, backend=CairoMakie)
fig_full |> f -> save(joinpath(path_output, "phase_diagram_W_sample.svg"), f; px_per_unit=2.0, backend=CairoMakie)
fig_a1278 |> f -> save(joinpath(path_output, "contrast_[a12=78a0].png"), f; px_per_unit=2.0, backend=CairoMakie)
fig_a1278 |> f -> save(joinpath(path_output, "contrast_[a12=78a0].svg"), f; px_per_unit=2.0, backend=CairoMakie)
fig_ctrs_164 |> f -> save(joinpath(path_output, "phase_diagram_C_sample.png"), f; px_per_unit=2.0, backend=CairoMakie)
fig_ctrs_164 |> f -> save(joinpath(path_output, "phase_diagram_C_sample.svg"), f; px_per_unit=2.0, backend=CairoMakie)
fig_ctrs_164 |> f -> save(joinpath(path_output, "phase_diagram_C_sample.pdf"), f; px_per_unit=2.0, backend=CairoMakie)
