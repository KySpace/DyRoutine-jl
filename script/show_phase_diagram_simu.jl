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

set_theme!(fonts = (; regular = raw"C:/Windows/Fonts/arial.ttf", italic = raw"C:/Windows/Fonts/ariali.ttf"))

# commit #706657ff45a8cf83882d67b6c273d63425d9966d
path_demo = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Demo"
path_output = joinpath(path_demo, "33.DualSS.PhaseDiagram.CWZ.Size")

path_simu_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS\Samples\[07.01].Weijing"
path_simu_w = joinpath(path_simu_root, "[07.06] weight recalculation wolfram")
path_simu_c = joinpath(path_simu_root, "[07.05] contrast without blur")
path_simu_z = joinpath(path_simu_root, "[07.01] vert sepr")
path_simu_sample_xy = joinpath(path_simu_root, "[07.02] density profiles")

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
    df_nondupl = @pipe average_duplicate_points(df, header_dupl) |> Float64.(_)
    col_nondupl = @pipe df_nondupl |> Tables.columntable(_) |> Tuple(_)
    (; df=df_nondupl, col_nondupl)
end

(df_ctrs, (a12_ctrs, a22_ctrs, vec_ctrs_162, vec_ctrs_164)) = make_dataframe_nondupl(path_simu_c, ["C_coarse.txt", "C_coarse2.txt", "C_fine.txt", "C_precise.txt"], [:a12, :a22, :contrast_162, :contrast_164]; header_dupl=[:a12, :a22])
(df_sepr, (a12_sepr, a22_sepr, vec_sepr)) = make_dataframe_nondupl(path_simu_z, ["dz_coarse.txt"], [:a12, :a22, :vert_sepr]; header_dupl=[:a12, :a22])
# (df_wght, (a12_wght, a22_wght, vec_wght_162, vec_wght_164)) = make_dataframe_nondupl(path_simu_w, ["W_coarse.txt", "W_coarse2.txt", "W_fine.txt", "W_precise.txt"], [:a12, :a22, :weightsp_162, :weightsp_164]; header_dupl=[:a12, :a22])
(df_wght, (a12_wght, a22_wght, vec_wght_162, vec_wght_164, _, _)) = make_dataframe_nondupl(path_simu_w, ["pdspw1a0.csv", "pdspw2a0.csv", "pdspwline.csv", "pdspwfine.csv"], [:a12, :a22, :weightsp_162, :weightsp_164, :err_weightsp_162, :err_weightsp_164]; delim=',', header_dupl=[:a12, :a22], skipto=1)
df_ctrs_a12_78 = @pipe df_ctrs[df_ctrs.a12.==78, :] |> sort!(_, :a22)
df_wght_a12_78 = @pipe df_wght[df_wght.a12.==78, :] |> sort!(_, :a22)
df_sepr_a12_78 = @pipe df_sepr[df_sepr.a12.==78, :] |> sort!(_, :a22)

##
clrrng_c = (0, 1)
clrrng_w = (0, 0.45) # extrema(vcat(vec(wght_162_q), vec(wght_164_q)))
clrrng_z = (0, 2.2) #

ntpl_ctrs_162 = interpolate(a22_ctrs, a12_ctrs, vec_ctrs_162; derivatives=true)
ntpl_ctrs_164 = interpolate(a22_ctrs, a12_ctrs, vec_ctrs_164; derivatives=true)
ntpl_wght_162 = interpolate(a22_wght, a12_wght, vec_wght_162; derivatives=true)
ntpl_wght_164 = interpolate(a22_wght, a12_wght, vec_wght_164; derivatives=true)
ntpl_sepr = interpolate(a22_sepr, a12_sepr, vec_sepr; derivatives=false)

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
sample_tripts = [
    (96, 78, BezierPath("m 0.0,0.0 l 0.93134,1.651 h -1.905 z"; fit=true, flipy=true), colorant"#6f5b98", colorant"#b8a2d6"),
    (a22_roton_instab, 78, BezierPath("m 0.0,0.0 l -0.93134,1.650997 l -0.97102,-1.650997 z"; fit=true, flipy=true), colorant"#555555", colorant"#cfcfcf"),
    (104, 78, BezierPath("m 0.0,0.0 l 1.01335,-1.016008 l 1.016,1.016008 l -1.016,1.01599 z"; fit=true, flipy=true), colorant"#9b7118", colorant"#e3c27a"),
]
sample_sepr = [
    (96, 78, BezierPath("m 0.0,0.0 l 0.93134,1.651 h -1.905 z"; fit=true, flipy=true), colorant"#6f5b98", colorant"#b8a2d6"),
    (104, 78, BezierPath("m 0.0,0.0 l 1.01335,-1.016008 l 1.016,1.016008 l -1.016,1.01599 z"; fit=true, flipy=true), colorant"#9b7118", colorant"#e3c27a"),
]

## Visualization: Heatmap on contrast and wght
fig_full = Figure();
kwargs_style_common = (; spinewidth=0.212/0.132*1.0, 
                         leftspinecolor=colorant"#4a4a4a", rightspinecolor=colorant"#4a4a4a", topspinecolor=colorant"#4a4a4a", bottomspinecolor=colorant"#4a4a4a",
                         )
kwargs_axis_common = (; xlabelsize=18.67, 
                        ylabelsize=18.67, 
                        xlabelfont="Arial", 
                        ylabelfont="Arial", 
                        xticklabelsize=16, 
                        yticklabelsize=16, 
                        xtickalign=1, ytickalign=1, 
                        xminortickalign=1, yminortickalign=1, 
                        xgridvisible=false, ygridvisible=false,
                        xtickcolor=colorant"#4a4a4a", ytickcolor=colorant"#4a4a4a",
                        xlabelcolor=colorant"#222222", ylabelcolor=colorant"#222222",
                        xticklabelcolor=colorant"#333333", yticklabelcolor=colorant"#333333",
                        xminortickcolor=colorant"#4a4a4a", yminortickcolor=colorant"#4a4a4a",
                        xtickwidth=0.212/0.132*1.0, ytickwidth=0.212/0.132*1.0,
                        xminortickwidth=0.169/0.132*1.0, yminortickwidth=0.169/0.132*1.0,
                        kwargs_style_common...)
kwargs_clrbar_common = (; labelsize=18.67,
                          ticklabelsize=16,
                          tickcolor=colorant"#4a4a4a",
                          labelcolor=colorant"#222222",
                          ticklabelcolor=colorant"#333333",
                          tickwidth=0.212/0.132*1.0,
                          kwargs_style_common...)
kwargs_marker_common = (; strokewidth=2*0.296/0.265, markersize=16*1.902/2.117)
Label(fig_full[1, 0]; text=L"^{162}\text{Dy}", valign=:center, halign=:center, fontsize=18.67)
Label(fig_full[2, 0]; text=L"^{164}\text{Dy}", valign=:center, halign=:center, fontsize=18.67)
Label(fig_full[0, 1]; text="contrast", valign=:center, halign=:center, font=:bold)
Label(fig_full[0, 2]; text="side peak weight", valign=:center, halign=:center, font=:bold)
Label(fig_full[0, 3]; text="vertical separation (μm)", valign=:center, halign=:center, font=:bold)
ax_ctrs_162 = Axis(fig_full[1, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_wght_162 = Axis(fig_full[1, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_ctrs_164 = Axis(fig_full[2, 1]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_wght_164 = Axis(fig_full[2, 2]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
ax_sepr = Axis(fig_full[1:2, 3]; ylabel=L"a_{12} \; (a_0)", xlabel=L"a_{22} \; (a_0)", kwargs_axis_common..., aspect=DataAspect());
clrmp_162 = gen_clrmap_solo(hue_theme_istp["162"])
clrmp_164 = gen_clrmap_solo(hue_theme_istp["164"])
clrmp_sepr = gen_clrmap_solo(293)
clrrng_sepr = (0, 2.0)

hm_c1 = heatmap!(ax_ctrs_162, a22_g, a12_g, ctrs_162_q; colormap=clrmp_162, colorrange=clrrng_c, rasterize=true)
hm_w1 = heatmap!(ax_wght_162, a22_g, a12_g, wght_162_q; colormap=clrmp_162, colorrange=clrrng_w, rasterize=true)
hm_c2 = heatmap!(ax_ctrs_164, a22_g, a12_g, ctrs_164_q; colormap=clrmp_164, colorrange=clrrng_c, rasterize=true)
hm_w2 = heatmap!(ax_wght_164, a22_g, a12_g, wght_164_q; colormap=clrmp_164, colorrange=clrrng_w, rasterize=true)
hm_sepr = heatmap!(ax_sepr, a22_g, a12_g, sepr_q; colormap=clrmp_sepr, colorrange=clrrng_sepr, rasterize=true)

Colorbar(fig_full[3, 1], hm_c1; vertical=false, label=L"C");
Colorbar(fig_full[4, 1], hm_c2; vertical=false, label=L"C");
Colorbar(fig_full[3, 2], hm_w1; vertical=false, label=L"W");
Colorbar(fig_full[4, 2], hm_w2; vertical=false, label=L"W");
Colorbar(fig_full[3, 3], hm_sepr; vertical=false, label=L"\Delta z");
colsize!(fig_full.layout, 0, 20)
colsize!(fig_full.layout, 1, 300)
colsize!(fig_full.layout, 2, 300)
colsize!(fig_full.layout, 3, 300)
colgap!(fig_full.layout, 10)
rowsize!(fig_full.layout, 1, 360)
rowsize!(fig_full.layout, 2, 360)
fig_full |> resize_to_layout!
fig_full |> display
for format in ["png", "svg"]
    fig_full |> f -> save(joinpath(path_output, "phase_diagram_W_sample.$format"), f; px_per_unit=2.0, backend=CairoMakie)
end

## Visualization: Heatmap and Linecut for one property
clr_lines_istp = [
    colorant"#c44e52",
    colorant"#4c72b0"
]
clr_stroke_istp = [
    colorant"#b14549",
    colorant"#3e6198"
]
clr_marker_face_istp = [
    colorant"#e6a0a3",
    colorant"#9bb5dd",
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

for (
    name_prop,
    abbr_prop,
    info_sample,
    label_axis_hm,
    label_axis_line,
    clrrng_map,
    clr_hm,
    labels_hm,
    kwargs_axis_prop,
    clrs_line,
    df_prop_line,
    prop_variant_line,
    labels_line,
    kwargs_axis_line,
    prop_variant_map,
    lim_prop,
    ticks_prop,
    zoom
) in [
    (
        name_prop="Contrast",
        abbr_prop="C",
        info_sample=sample_tripts,
        label_axis_hm="Contrast",
        label_axis_line="Contrast",
        clrrng_map=(0, 1.0),
        clr_hm=(; hue=84, prescale=(t -> t^5)),
        labels_hm=["162", "164"],
        kwargs_axis_prop=(;),
        clrs_line=(; line=clr_lines_istp, stroke=clr_stroke_istp, markerface=clr_marker_face_istp),
        df_prop_line=df_ctrs_a12_78,
        prop_variant_line=[:contrast_162, :contrast_164],
        labels_line=["162", "164"],
        kwargs_axis_line=(;),
        prop_variant_map=[ctrs_162_q, ctrs_164_q],
        lim_prop=(-0.05, 1.05),
        ticks_prop=0:0.2:1,
        zoom=(;
            lim_prop=(-0.02, 0.52),
            ticks_prop=0:0.2:1,
            kwargs_a1278=(; width=140, height=80, halign=0.13, valign=0.35),
            lim_x=(99.8, 99.9),
        ),
    ),
    (
        name_prop="Side-peak weight",
        abbr_prop="W",
        info_sample=sample_tripts,
        label_axis_hm=rich("Side-peak weight ", rich("W", font = :italic), subscript("162"),),
        label_axis_line=rich(rich("W", font = :italic), subscript("162,164")),
        clrrng_map=(0, 0.5),
        clr_hm=(; hue=196, prescale=(t -> t^3)),
        labels_hm=["162", "164"],
        kwargs_axis_prop=(;),
        clrs_line=(; line=clr_lines_istp, stroke=clr_stroke_istp, markerface=clr_marker_face_istp),
        df_prop_line=df_wght_a12_78,
        prop_variant_line=[:weightsp_162, :weightsp_164],
        labels_line=["162", "164"],
        kwargs_axis_line=(; width=320),
        prop_variant_map=[wght_162_q, wght_164_q],
        lim_prop=(-0.05, 0.55),
        ticks_prop=0:0.1:0.6,
        zoom=nothing,
    ),
    (
        name_prop="Vertical separation",
        abbr_prop="Z",
        info_sample=sample_sepr,
        label_axis_hm=rich("Vertical separation ", rich("δz", font=:italic), " (μm)"),
        label_axis_line=rich(rich("δz", font=:italic), " (μm)"),
        clrrng_map=(0.3, 2.0),
        clr_hm=(; hue=293, prescale=(t -> t^2)),
        labels_hm=[""],
        kwargs_axis_prop=(; width=500, height=500),
        clrs_line=(; line=[Oklch(0.56, 0.08, 293)], stroke=[Oklch(0.56, 0.08, 293)], markerface=[Oklch(0.80, 0.066, 293)]),
        df_prop_line=df_sepr_a12_78,
        prop_variant_line=[:vert_sepr],
        labels_line=[nothing],
        kwargs_axis_line=(;),
        prop_variant_map=[sepr_q],
        lim_prop=(-0.2, 2.2),
        ticks_prop=0:0.5:2,
        zoom=nothing,
    ),
]
    clrmp_map = gen_clrmap_parabola(clr_hm.hue, 0.58, 0.06, 0.44; light_max=0.99, chroma_lightmax=0.008, prescale=clr_hm.prescale)
    for (i, prop) in enumerate(prop_variant_map)
        fig_prop = Figure(;backgroundcolor = :transparent)
        ax_prop = Axis(fig_prop[1, 1]; ylabel=rich(rich("a", font = :italic), subscript("12"), "\u2005(", rich("a", font = :italic), subscript("0"), ")",), 
                                        xlabel=rich(rich("a", font = :italic), subscript("22"), "\u2005(", rich("a", font = :italic), subscript("0"), ")",), kwargs_axis_common..., width=280, height=280, kwargs_axis_prop...)
        ax_prop |> empty!
        hm = heatmap!(ax_prop, a22_g, a12_g, prop; colormap=clrmp_map, colorrange=clrrng_map, rasterize=true)
        for (a22, a12, marker, clr_stroke, clr_face) in info_sample
            scatter!(ax_prop, [a22], [a12];
                color=clr_face, strokecolor=clr_stroke, marker=marker, kwargs_marker_common...)
        end
        Colorbar(fig_prop[1, 2], hm; vertical=true, label=label_axis_hm, labelrotation=-π / 2, kwargs_clrbar_common...)
        limits!(ax_prop, (94, 106), (70, 86))
        ax_prop.xticks = 90:2:110
        ax_prop.yticks = 70:2:106
        fig_prop |> resize_to_layout!
        fig_prop |> display
        for format in ["png", "svg", "pdf"]
            fig_prop |> f -> save(joinpath(path_output, "phase_diagram_$(abbr_prop)_$(labels_hm[i])_sample.$format"), f; px_per_unit=2.0, backend=CairoMakie)
        end
    end

    fig_a1278 = Figure(;backgroundcolor = :transparent)
    ax_a1278 = Axis(fig_a1278[1, 1]; xticksmirrored = true, yticksmirrored = true, ylabel=label_axis_line, xlabel=rich(rich("a", font = :italic), subscript("22"), "\u2005(", rich("a", font = :italic), subscript("0"), ")",), width=400, height=150, kwargs_axis_common..., kwargs_axis_line...)

    ax_a1278.xticks = 90:1:110
    ax_a1278.xminorticks = IntervalsBetween(2)
    ax_a1278.xminorticksvisible = true
    ax_a1278.yticks = ticks_prop

    lim_axis_line = (95.8, 101.2)
    limits!(ax_a1278, lim_axis_line, lim_prop)
    kwargs_lines = i -> (; linewidth=1, color=clrs_line.line[i], strokecolor=clrs_line.stroke[i], strokewidth=1.2, markersize=8, marker=(i == 1 ? :rect : :circle), markercolor=clrs_line.markerface[i])
    if !isnothing(zoom)
        Box(fig_a1278[1, 1]; color=:white, zoom.kwargs_a1278..., strokewidth=0)
        Box(fig_a1278[1, 1]; color=(Oklch(0.90, 0.005, 192), 0.2), zoom.kwargs_a1278..., strokewidth=0)
        ax_a1278_zoom = Axis(fig_a1278[1, 1]; backgroundcolor=:white, zoom.kwargs_a1278..., kwargs_axis_common..., xticklabelsize=13, yticklabelsize=13)
        ax_a1278_zoom.xticks = zoom.lim_x |> collect
        ax_a1278_zoom.xminorticks = 99.80:0.02:99.90
        ax_a1278_zoom.xminorticksvisible = true
        ax_a1278_zoom.yticks = zoom.ticks_prop
        limits!(ax_a1278_zoom, zoom.lim_x, zoom.lim_prop)
        lines!(ax_a1278_zoom, [a22_roton_instab, a22_roton_instab], [zoom.lim_prop[1], zoom.lim_prop[2]*0.90]; color=info_sample[2][4], linewidth=0.8, linestyle=:dash)
        vspan!(ax_a1278, zoom.lim_x[1], zoom.lim_x[2]; color=(Oklch(0.90, 0.005, 192), 0.5))
        for (i, col) in enumerate(prop_variant_line)
            scatterlines!(ax_a1278_zoom, df_prop_line.a22, df_prop_line[!, col]; kwargs_lines(i)...)
        end
        y_marker = zoom.lim_prop[2] * 0.90
        let (a22, _, marker, clr_stroke, clr_face) = info_sample[2]
            scatter!(ax_a1278_zoom, [a22], [y_marker];
                color=clr_face, strokecolor=clr_stroke, marker=marker, kwargs_marker_common...)
        end
    else
        lines!(ax_a1278, [a22_roton_instab, a22_roton_instab], [lim_prop[1], lim_prop[2]*0.84]; color=colorant"#d0d0d0", linewidth=0.8)
        y_marker = lim_prop[2] * 0.90
        let (a22, _, marker, clr_stroke, clr_face) = info_sample[2]
            scatter!(ax_a1278, [a22], [y_marker];
                color=clr_face, strokecolor=clr_stroke, marker=marker, kwargs_marker_common...)
        end
    end
    for (i, col) in enumerate(prop_variant_line)
        global sel_a22 = lim_axis_line[1] .<= df_prop_line.a22 .<= lim_axis_line[2]
        scatterlines!(ax_a1278, df_prop_line.a22[sel_a22], df_prop_line[sel_a22, col]; kwargs_lines(i)..., label=labels_line[i])
    end
    axislegend(ax_a1278; position=:rt, framevisible=false, labelsize=16)
    fig_a1278 |> resize_to_layout!
    fig_a1278 |> display
    for format in ["png", "svg", "pdf"]
        fig_a1278 |> f -> save(joinpath(path_output, "linecut_$(abbr_prop)_[a12=78a0].$format"), f; px_per_unit=2.0, backend=CairoMakie)
    end

end
