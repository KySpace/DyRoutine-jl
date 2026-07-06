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
path_simu_w = joinpath(path_simu_root, "[07.06] weight recalculation wolfram")
path_simu_c = joinpath(path_simu_root, "[07.05] contrast without blur")
path_simu_z = joinpath(path_simu_root, "[07.01] vert sepr")
path_simu_sample_xy = joinpath(path_simu_root, "[07.02] density profiles")
path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Demo"
# commit #93504925a6f1b5e790e838a1a66c2e5f653afdf5
path_output = joinpath(path_demo, "29.DualSS.PhaseDiagram.CWZ.Wolfram.WithZ")
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

function make_dataframe_nondupl(path, filenames, header; delim='\t', header_dupl=[:a12, :a22], skipto=2)
    df = @pipe [
        CSV.read(joinpath(path, fn), DataFrame; delim, header, skipto)
        for fn in filenames
    ] |> vcat(_...)
    df_nondupl = average_duplicate_points(df, header_dupl)
    col_nondupl = @pipe df_nondupl |> Tables.columntable(_) |> Tuple(_)
    (; df=df_nondupl, col_nondupl)
end

(df_ctrs, (a12_ctrs, a22_ctrs, vec_ctrs_162, vec_ctrs_164)) = make_dataframe_nondupl(path_simu_c, ["C_coarse.txt", "C_coarse2.txt", "C_fine.txt", "C_precise.txt"], [:a12, :a22, :contrast_162, :contrast_164]; header_dupl=[:a12, :a22])
(df_sepr, (a12_sepr, a22_sepr, vec_sepr)) = make_dataframe_nondupl(path_simu_c, ["dz_coarse.txt"], [:a12, :a22, :vert_sepr]; header_dupl=[:a12, :a22])
# (df_wght, (a12_wght, a22_wght, vec_wght_162, vec_wght_164)) = make_dataframe_nondupl(path_simu_w, ["W_coarse.txt", "W_coarse2.txt", "W_fine.txt", "W_precise.txt"], [:a12, :a22, :weightsp_162, :weightsp_164]; header_dupl=[:a12, :a22])
(df_wght, (a12_wght, a22_wght, vec_wght_162, vec_wght_164, _, _)) = make_dataframe_nondupl(path_simu_w, ["pdspw1a0.csv", "pdspw2a0.csv", "pdspwline.csv", "pdspwfine.csv"], [:a12, :a22, :weightsp_162, :weightsp_164, :err_weightsp_162, :err_weightsp_164]; delim=',', header_dupl=[:a12, :a22], skipto=1)
df_ctrs_a12_78 = @pipe df_ctrs[df_ctrs.a12.==78, :] |> sort!(_, :a22)
df_wght_a12_78 = @pipe df_wght[df_wght.a12.==78, :] |> sort!(_, :a22)

##
clrrng_c = (0, 1)
clrrng_w = (0, 0.45) # extrema(vcat(vec(wght_162_q), vec(wght_164_q)))
coor_ctrs = hcat(a22_ctrs, a12_ctrs)'
coor_wght = hcat(a22_wght, a12_wght)'

ntpl_ctrs_162 = interpolate(a22_ctrs, a12_ctrs, vec_ctrs_162; derivatives=true)
ntpl_ctrs_164 = interpolate(a22_ctrs, a12_ctrs, vec_ctrs_164; derivatives=true)
ntpl_wght_162 = interpolate(a22_wght, a12_wght, vec_wght_162; derivatives=true)
ntpl_wght_164 = interpolate(a22_wght, a12_wght, vec_wght_164; derivatives=true)
ntpl_sepr = interpolate(a22_sepr, a12_sepr, vec_sepr; derivatives=true)

a22_g = range(88, 108; length=200)
a12_g = range(70, 96; length=200)
a22_q = vec([xi for xi in a22_g, yi in a12_g])
a12_q = vec([yi for xi in a22_g, yi in a12_g])

ctrs_162_q = @pipe ntpl_ctrs_162(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))
ctrs_164_q = @pipe ntpl_ctrs_164(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))
wght_162_q = @pipe ntpl_wght_162(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))
wght_164_q = @pipe ntpl_wght_164(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))
sepr_q = @pipe ntpl_sepr(a22_q, a12_q; method=Sibson()) |> reshape(_, length(a22_g), length(a12_g))

a22_roton_instab = 99.8632
sample_ctrs = [
    (96, 78, :utriangle, colorant"rgb(107, 93, 147)", colorant"rgb(179, 162, 209)"),
    (a22_roton_instab, 78, :dtriangle, colorant"rgb(107, 107, 107)", colorant"rgb(217, 217, 217)"),
    (104, 78, :diamond, colorant"rgb(144, 113, 45)", colorant"rgb(217, 195, 131)"),
]

## Visualization: Heatmap on contrast and wght
fig_full = Figure();
kwargs_axis_common = (; xlabelsize=16, ylabelsize=16, xlabelfont="Helvetica World", ylabelfont="Helvetica World", xticklabelsize=14, yticklabelsize=14, xtickalign=1, ytickalign=1, xminortickalign=1, yminortickalign=1, xgridvisible=false, ygridvisible=false)
Label(fig_full[1, 0]; text=L"^{162}\text{Dy}", valign=:center, halign=:center, fontsize=16)
Label(fig_full[2, 0]; text=L"^{164}\text{Dy}", valign=:center, halign=:center, fontsize=16)
Label(fig_full[0, 1]; text="contrast", valign=:center, halign=:center, font=:bold)
Label(fig_full[0, 2]; text="side peak weight", valign=:center, halign=:center, font=:bold)
ax_ctrs_162 = Axis(fig_full[1, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_wght_162 = Axis(fig_full[1, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_ctrs_164 = Axis(fig_full[2, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_wght_164 = Axis(fig_full[2, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_sepr = Axis(fig_full[3, 1:2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
clrmp_162 = gen_clrmap_solo(hue_theme_istp["162"])
clrmp_164 = gen_clrmap_solo(hue_theme_istp["164"])

hm_c1 = heatmap!(ax_ctrs_162, a22_g, a12_g, ctrs_162_q; colormap=clrmp_162, colorrange=clrrng_c, rasterize=true)
hm_w1 = heatmap!(ax_wght_162, a22_g, a12_g, wght_162_q; colormap=clrmp_162, colorrange=clrrng_w, rasterize=true)
hm_c2 = heatmap!(ax_ctrs_164, a22_g, a12_g, ctrs_164_q; colormap=clrmp_164, colorrange=clrrng_c, rasterize=true)
hm_w2 = heatmap!(ax_wght_164, a22_g, a12_g, wght_164_q; colormap=clrmp_164, colorrange=clrrng_w, rasterize=true)

Colorbar(fig_full[3, 1], hm_c1; vertical=false, label=L"C");
Colorbar(fig_full[4, 1], hm_c2; vertical=false, label=L"C");
Colorbar(fig_full[3, 2], hm_w1; vertical=false, label=L"W");
Colorbar(fig_full[4, 2], hm_w2; vertical=false, label=L"W");
colsize!(fig_full.layout, 0, 20)
colsize!(fig_full.layout, 1, 300)
colsize!(fig_full.layout, 2, 300)
colgap!(fig_full.layout, 10)
rowsize!(fig_full.layout, 1, 360)
rowsize!(fig_full.layout, 2, 360)
fig_full |> resize_to_layout!
fig_full |> display
for format in ["png", "svg"]
    fig_full |> f -> save(joinpath(path_output, "phase_diagram_W_sample.$format"), f; px_per_unit=2.0, backend=CairoMakie)
end

## Visualization: Heatmap and Linecut for one property
clr_lines = [
    colorant"rgb(157, 76, 76)",
    colorant"rgb(72, 93, 144)"
]
clr_marker_face = [
    colorant"rgb(214, 163, 164)",
    colorant"rgb(164, 181, 217)",
]

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

for (name_prop, abbr_prop, clrrng_map, df_prop, df_prop_line, prop_variant_line, prop_variant_map, lim_prop, lim_prop_zoom, ticks_prop, ticks_prop_zoom, kwargs_a1278_zoom, info_sample, clr_hm) in [
    (
        name_prop="contrast",
        abbr_prop="C",
        clrrng_map=(0, 1.0),
        df_prop=df_ctrs,
        df_prop_line=df_ctrs_a12_78,
        prop_variant_line=[:contrast_162, :contrast_164],
        prop_variant_map=[ctrs_162_q, ctrs_164_q],
        lim_prop=(-0.05, 1.05),
        lim_prop_zoom=(-0.02, 0.52),
        ticks_prop=0:0.2:1,
        ticks_prop_zoom=0:0.2:1,
        kwargs_a1278_zoom=(; width=140, height=80, halign=0.13, valign=0.35),
        info_sample=sample_ctrs,
        clr_hm=(; hue=196, prescale=(t -> t^5))
    ),
    (
        name_prop="sidepeak weight",
        abbr_prop="W",
        clrrng_map=(0, 0.5),
        df_prop=df_wght,
        df_prop_line=df_wght_a12_78,
        prop_variant_line=[:weightsp_162, :weightsp_164],
        prop_variant_map=[wght_162_q, wght_164_q],
        lim_prop=(-0.05, 0.55),
        lim_prop_zoom=(-0.02, 0.32),
        ticks_prop=0:0.1:0.6,
        ticks_prop_zoom=0:0.1:0.6,
        kwargs_a1278_zoom=(; width=140, height=80, halign=0.13, valign=0.35),
        info_sample=sample_ctrs,
        clr_hm=(; hue=84, prescale=(t -> t))
    ),
    (
        name_prop="sidepeak weight",
        abbr_prop="W",
        clrrng_map=(0, 0.5),
        df_prop=df_wght,
        df_prop_line=df_wght_a12_78,
        prop_variant_line=[:weightsp_162, :weightsp_164],
        prop_variant_map=[wght_162_q, wght_164_q],
        lim_prop=(-0.05, 0.55),
        lim_prop_zoom=(-0.02, 0.32),
        ticks_prop=0:0.1:0.6,
        ticks_prop_zoom=0:0.1:0.6,
        kwargs_a1278_zoom=(; width=140, height=80, halign=0.13, valign=0.35),
        info_sample=sample_ctrs,
        clr_hm=(; hue=84, prescale=(t -> t))
    ),
]
    fig_prop = Figure()
    ax_prop = Axis(fig_prop[1, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., width=280, height=280)

    fig_a1278 = Figure()
    ax_a1278 = Axis(fig_a1278[1, 1]; ylabel=name_prop, xlabel=L"a_{22} \; (a_0)", width=400, height=150, kwargs_axis_common...)

    clrmp_map = gen_clrmap_parabola(clr_hm.hue, 0.58, 0.06, 0.55; light_max=0.97, chroma_lightmax=0.008, prescale=clr_hm.prescale)
    Box(fig_a1278[1, 1]; color=:white, kwargs_a1278_zoom..., strokewidth=0)
    Box(fig_a1278[1, 1]; color=(Oklch(0.90, 0.005, 192), 0.2), kwargs_a1278_zoom..., strokewidth=0)
    ax_a1278_zoom = Axis(fig_a1278[1, 1]; backgroundcolor=:white, kwargs_a1278_zoom..., kwargs_axis_common..., xticklabelsize=13, yticklabelsize=13)
    lim_x_zoom = (99.8, 99.9)
    ax_a1278_zoom.xticks = lim_x_zoom |> collect
    ax_a1278_zoom.xminorticks = 99.80:0.02:99.90
    ax_a1278_zoom.xminorticksvisible = true
    ax_a1278.xticks = 90:1:110
    ax_a1278.xminorticks = IntervalsBetween(2)
    ax_a1278.xminorticksvisible = true
    ax_a1278.yticks = ticks_prop

    for (i, prop) in enumerate(prop_variant_map)
        ax_prop |> empty!
        istp = i == 1 ? "162" : "164"
        hm = heatmap!(ax_prop, a22_g, a12_g, prop; colormap=clrmp_map, colorrange=clrrng_map, rasterize=true)
        for (a22, a12, marker, clr_stroke, clr_face) in info_sample
            scatter!(ax_prop, [a22], [a12];
                color=clr_face, strokecolor=clr_stroke, strokewidth=1.5, marker=marker, markersize=12)
        end
        Colorbar(fig_prop[1, 2], hm; vertical=true, label=name_prop, labelrotation=-π / 2)
        limits!(ax_prop, (90, 106), (70, 90))
        ax_prop.xticks = 90:2:110
        ax_prop.yticks = 70:2:106
        fig_prop |> resize_to_layout!
        fig_prop |> display
        for format in ["png", "svg", "pdf"]
            fig_prop |> f -> save(joinpath(path_output, "phase_diagram_$(abbr_prop)_$(istp)_sample.$format"), f; px_per_unit=2.0, backend=CairoMakie)
        end
    end

    limits!(ax_a1278, (95.8, 101.2), lim_prop)
    limits!(ax_a1278_zoom, lim_x_zoom, lim_prop_zoom)
    vspan!(ax_a1278, 99.8, 99.9; color=(Oklch(0.90, 0.005, 192), 0.5))
    # vlines!(ax_a1278, a22_roton_instab; color=:mediumpurple4, linewidth=0.8)
    vlines!(ax_a1278_zoom, a22_roton_instab; color=info_sample[2][4], linewidth=0.8, linestyle=:dash)
    kwargs_lines = i -> (; linewidth=1, color=clr_lines[i], strokecolor=clr_lines[i], strokewidth=1, markersize=8, marker=(i == 1 ? :rect : :circle), markercolor=clr_marker_face[i])
    for (i, col) in enumerate(prop_variant_line)
        istp = i == 1 ? "162" : "164"
        scatterlines!(ax_a1278, df_prop_line.a22, df_prop_line[!, col]; kwargs_lines(i)..., label=L"^{%$istp}\text{Dy}")
        scatterlines!(ax_a1278_zoom, df_prop_line.a22, df_prop_line[!, col]; kwargs_lines(i)...)
    end
    y_marker_zoom = @pipe df_prop_line[lim_x_zoom[1].<=df_prop_line.a22.<=lim_x_zoom[2], :] |>
                          map(symb -> _[!, symb], prop_variant_line) |> vcat(_...) |> maximum |> _ * 0.95
    let (a22, _, marker, clr_stroke, clr_face) = info_sample[2]
        scatter!(ax_a1278_zoom, [a22], [y_marker_zoom];
            color=clr_face, strokecolor=clr_stroke, strokewidth=1, marker=marker, markersize=8)
    end
    axislegend(ax_a1278; position=:rt, framevisible=false, labelsize=14)
    fig_a1278 |> resize_to_layout!
    fig_a1278 |> display
    for format in ["png", "svg", "pdf"]
        fig_a1278 |> f -> save(joinpath(path_output, "linecut_$(abbr_prop)_[a12=78a0].$format"), f; px_per_unit=2.0, backend=CairoMakie)
    end

end
