include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_data = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data"
path_result = joinpath(path_data, "Result")

val_istp = (Symbol("162"), Symbol("164"))
loadcfgs = (:SCS, :DCS, :DIS, :DDM)
font_labels = "Noto Sans Math"
fontsize_labels = 22 / 3 # 5.5 pt with CairoMakie's default pt_per_unit = 0.75

# CairoMakie exports 92 units as 69 points (24.3 mm) in PDF.
fig = Figure(size=(106, 38), fontsize=fontsize_labels, figure_padding=2)
border_color = RGBAf(0.25, 0.25, 0.25, 0.45)
ax = Axis(fig[1, 1];
    limits=(0.0, 6.0, 0.0, 2.4),
    xticksvisible=false,
    xticklabelsvisible=false,
    yticksvisible=false,
    yticklabelsvisible=false,
    xgridvisible=false,
    ygridvisible=false,
    bottomspinecolor=border_color,
    leftspinecolor=border_color,
    topspinecolor=border_color,
    rightspinecolor=border_color,
    spinewidth=0.75,
    backgroundcolor=RGBAf(1, 1, 1, 0.86),
)

pos_x = range(1.75, 5.25; length=length(loadcfgs))
pos_y = (1.20, 0.45)
for (x, loadcfg) in zip(pos_x, loadcfgs)
    text!(ax, x, 1.90;
        text=string(loadcfg), align=(:center, :center), font=font_labels)
end
for (y, istp) in zip(pos_y, val_istp)
    text!(ax, 0.16, y;
        text=string(istp), align=(:left, :center), font=font_labels)
    for (x, loadcfg) in zip(pos_x, loadcfgs)
        style = dualmot_curve_style((; loadcfg, istp))
        scatter!(ax, [x], [y];
            color=style.markercolor,
            marker=style.marker,
            markersize=style.markersize,
            strokecolor=style.strokecolor,
            strokewidth=style.strokewidth)
    end
end

mkpath(path_result)
for format in ("svg", "png", "pdf")
    save(joinpath(path_result, "marker_sample.$format"), fig)
end
println("Saved marker sample to $path_result")
