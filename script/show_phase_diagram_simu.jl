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

path_simu_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Samples\[07.01].Weijing"
path_simu_w = joinpath(path_simu_root, "[07.06] weight recalculation")
path_simu_c = joinpath(path_simu_root, "[07.05] contrast without blur")
path_simu_z = joinpath(path_simu_root, "[07.02] density profiles")
path_simu_sample_xy = joinpath(path_simu_root, "[07.02] density profiles")
path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Demo"
# commit #c018bbf9368558cbb09a629dcdd8a39cda93bbeb
path_output = joinpath(path_demo, "25.DualSS.PhaseDiagram.CWZ")
isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)

function average_duplicate_points(df::AbstractDataFrame, header_dupl::AbstractVector{Symbol})
    header = Symbol.(names(df))
    header_missing = setdiff(header_dupl, header)
    isempty(header_missing) ||
        throw(ArgumentError("header_dupl contains columns not present in df: $(header_missing)"))

    ids_valid = [
        all(col -> begin
                val = row[col]
                !ismissing(val) && (!(val isa Number) || isfinite(val))
            end, header)
        for row in eachrow(df)
    ]
    df_valid = df[ids_valid, :]
    header_avg = setdiff(header, header_dupl)
    df_nondupl = combine(
        groupby(df_valid, header_dupl),
        header_avg .=> (col -> mean(Float64.(col))) .=> header_avg,
    )
    return select(df_nondupl, header)
end

function make_dataframe_nondupl(path, filenames, header; header_dupl=[:a12, :a22], skipto=2)
    df = @pipe [
                         CSV.read(joinpath(path, fn), DataFrame; delim='\t', header, skipto)
                         for fn in filenames
                     ] |> vcat(_...) 
    df_nondupl = average_duplicate_points(df, header_dupl)
    col_nondupl = @pipe df_nondupl |> Tables.columntable(_) |> Tuple(_)
    (;df=df_nondupl, col_nondupl)
end

(df_ctrs, (a12_ctrs, a22_ctrs, df_contrast_162, df_contrast_164)) = make_dataframe_nondupl(path_simu_c, ["C_coarse.txt", "C_coarse2.txt", "C_fine.txt", "C_precise.txt"], [:a12, :a22, :contrast_162, :contrast_164]; header_dupl=[:a12, :a22])
(df_wght, (a12_wght, a22_wght, df_weightsp_162, df_weightsp_164)) = make_dataframe_nondupl(path_simu_w, ["W_coarse.txt", "W_coarse2.txt", "W_fine.txt", "W_precise.txt"], [:a12, :a22, :weightsp_162, :weightsp_164]; header_dupl=[:a12, :a22])
df_contrast_a12_78 = df_ctrs[df_ctrs.a12 .== 78, :]
df_weightsp_a12_78 = df_wght[df_wght.a12 .== 78, :]

##
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

a22_roton_instab = 99.8632
sample_contrast = [
    (96, 78, :utriangle, colorant"rgb(107, 93, 147)", colorant"rgb(179, 162, 209)"),
    (a22_roton_instab, 78, :dtriangle, colorant"rgb(107, 107, 107)", colorant"rgb(217, 217, 217)"),
    (104, 78, :diamond, colorant"rgb(144, 113, 45)", colorant"rgb(217, 195, 131)"),
]
##
fig_full = Figure();
kwargs_axis_common = (; xlabelsize=16, ylabelsize=16, xlabelfont="Helvetica World", ylabelfont="Helvetica World", xticklabelsize=14, yticklabelsize=14, xtickalign=1, ytickalign=1, xminortickalign=1, yminortickalign=1, xgridvisible=false, ygridvisible=false)
Label(fig_full[1, 0]; text=L"^{162}\text{Dy}", valign=:center, halign=:center, fontsize=16)
Label(fig_full[2, 0]; text=L"^{164}\text{Dy}", valign=:center, halign=:center, fontsize=16)
Label(fig_full[0, 1]; text="contrast", valign=:center, halign=:center, font=:bold)
Label(fig_full[0, 2]; text="side peak weight", valign=:center, halign=:center, font=:bold)
ax_contrast_162 = Axis(fig_full[1, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_weightsp_162 = Axis(fig_full[1, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_contrast_164 = Axis(fig_full[2, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_weightsp_164 = Axis(fig_full[2, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());

clr_lines = [
    colorant"rgb(157, 76, 76)",
    colorant"rgb(72, 93, 144)"
]
clr_marker_face = [
    colorant"rgb(214, 163, 164)",
    colorant"rgb(164, 181, 217)",
]


fig_ctrs_164 = Figure()
ax_contrast_sample = Axis(fig_ctrs_164[1, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., width=280, height=280);

fig_a1278 = Figure()
ax_a1278 = Axis(fig_a1278[1, 1]; ylabel="contrast", xlabel=L"a_{22} \; (a_0)", width=400, height=150, kwargs_axis_common...);
kwargs_a1278_zoom = (; width=140, height=80, halign=0.13, valign=0.35)
Box(fig_a1278[1, 1]; color=:white, kwargs_a1278_zoom..., strokewidth=0)
Box(fig_a1278[1, 1]; color=(Oklch(0.90, 0.005, 192), 0.2), kwargs_a1278_zoom..., strokewidth=0)
ax_a1278_zoom = Axis(fig_a1278[1, 1]; backgroundcolor=:white, kwargs_a1278_zoom..., kwargs_axis_common..., xticklabelsize=13, yticklabelsize=13);
ax_a1278_zoom.xticks = [99.8, 99.90]
ax_a1278_zoom.xminorticks = 99.80:0.02:99.90
ax_a1278_zoom.xminorticksvisible = true
ax_a1278.xticks = 90:1:110
ax_a1278.xminorticks = IntervalsBetween(2)
ax_a1278.xminorticksvisible = true
ax_a1278.yticks = 0:0.2:1

function gen_clrmap_parabola(hue, light_maxchroma, chroma_max, light_min; thres_alpha=0.0, alpha_base=1.0, light_max=1.0, chroma_lightmax=0, hue_range=(0, 0), prescale=(t -> t))
    clrmap = [
        begin
            t = prescale(t)
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
clrmp_turqoise = gen_clrmap_parabola(196, 0.58, 0.06, 0.55; hue_range=(0, 0), light_max=0.97, chroma_lightmax=0.008, thres_alpha=0, alpha_base=1.0, prescale=t -> t^5)

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
hm_cs = heatmap!(ax_contrast_sample, a22_g, a12_g, contrast_164_q; colormap=clrmp_turqoise, colorrange=clrrng_c, rasterize=true)
for (a22, a12, marker, clr_stroke, clr_face) in sample_contrast
    scatter!(ax_contrast_sample, [a22], [a12];
        color=clr_face, strokecolor=clr_stroke, strokewidth=1.5, marker=marker, markersize=12)
end

Colorbar(fig_full[3, 1], hm_c1; vertical=false, label=L"C");
Colorbar(fig_full[4, 1], hm_c2; vertical=false, label=L"C");
Colorbar(fig_full[3, 2], hm_w1; vertical=false, label=L"W");
Colorbar(fig_full[4, 2], hm_w2; vertical=false, label=L"W");
Colorbar(fig_ctrs_164[1, 2], hm_cs; vertical=true, label="contrast", labelrotation=-π / 2);
limits!(ax_contrast_sample, (95, 108), (70, 96))

colsize!(fig_full.layout, 0, 20)
colsize!(fig_full.layout, 1, 300)
colsize!(fig_full.layout, 2, 300)
colgap!(fig_full.layout, 10)
rowsize!(fig_full.layout, 1, 360)
rowsize!(fig_full.layout, 2, 360)

# colsize!(fig_ctrs_164.layout, 1, 320)
# rowsize!(fig_ctrs_164.layout, 1, 280)
xlims!(ax_contrast_sample, (90, 106))
ylims!(ax_contrast_sample, (70, 90))
ax_contrast_sample.xticks = 90:2:110
ax_contrast_sample.yticks = 70:2:106

xlims!(ax_a1278, (95.8, 101.2))
ylims!(ax_a1278, (-0.05, 1.05))
xlims!(ax_a1278_zoom, (99.8, 99.9))
ylims!(ax_a1278_zoom, (-0.02, 0.52))
vspan!(ax_a1278, 99.8, 99.9; color=(Oklch(0.90, 0.005, 192), 0.5))
# vlines!(ax_a1278, a22_roton_instab; color=:mediumpurple4, linewidth=0.8)
vlines!(ax_a1278_zoom, a22_roton_instab; color=sample_contrast[2][4], linewidth=0.8, linestyle=:dash)
kwargs_lines = i -> (; linewidth=1, color=clr_lines[i], strokecolor=clr_lines[i], strokewidth=1, markersize=8, marker=(i == 1 ? :rect : :circle), markercolor=clr_marker_face[i])
cut_c_162 = scatterlines!(ax_a1278, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_162; kwargs_lines(1)..., label=L"^{162}\text{Dy}")
cut_c_164 = scatterlines!(ax_a1278, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_164; kwargs_lines(2)..., label=L"^{164}\text{Dy}")
scatterlines!(ax_a1278_zoom, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_162; kwargs_lines(1)...)
scatterlines!(ax_a1278_zoom, df_contrast_a12_78.a22, df_contrast_a12_78.contrast_164; kwargs_lines(2)...)
let (a22, _, marker, clr_stroke, clr_face) = sample_contrast[2]
    scatter!(ax_a1278_zoom, [a22], [0.47];
        color=clr_face, strokecolor=clr_stroke, strokewidth=1, marker=marker, markersize=8)
end
axislegend(ax_a1278; position=:rt, framevisible=false, labelsize=14)


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
