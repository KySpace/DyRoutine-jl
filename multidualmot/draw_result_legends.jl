"""Export compact load-configuration/isotope marker keys for Result figures."""

path_result_legends = joinpath(path_root, "Result")
mkpath(path_result_legends)
val_istp_legend = (Symbol("162"), Symbol("164"))
legend_versions = (
    (name="DCS-DIS-DDM", loadcfgs=(:DCS, :DIS, :DDM)),
    (name="SCS-DCS", loadcfgs=(:SCS, :DCS)),
    (name="DIS-DDM", loadcfgs=(:DIS, :DDM)),
)

# CairoMakie uses 0.75 points per fontsize unit; 8 units gives 6 pt.
legend_fontsize = 6 / 0.75
pos_y_legend = (1.20, 0.45)
pos_x_start_legend = 1.75
pos_x_step_legend = (5.25 - 1.75) / 3
fig_height_legend = 38
axis_limits_y_legend = (0.0, 2.4)

function draw_result_legend(loadcfgs::Tuple, name::AbstractString)
    pos_x = [pos_x_start_legend + (idx - 1) * pos_x_step_legend
        for idx in eachindex(loadcfgs)]
    x_max = last(pos_x) + 0.75
    # Preserve the same data-to-output scale and gaps as the four-column key.
    fig_width = round(Int, (106 - 4) * x_max / 6 + 4)
    fig = Figure(size=(fig_width, fig_height_legend),
        fontsize=legend_fontsize, figure_padding=2)
    border_color = RGBAf(0.25, 0.25, 0.25, 0.45)
    ax = Axis(fig[1, 1];
        width=fig_width - 4,
        height=fig_height_legend - 4,
        halign=0,
        valign=0,
        tellwidth=false,
        tellheight=false,
        limits=(0.0, x_max, axis_limits_y_legend...),
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
    for (x, loadcfg) in zip(pos_x, loadcfgs)
        text!(ax, x, 1.90; text=string(loadcfg),
            align=(:center, :center), font=DUALMOT_FONT,
            fontsize=legend_fontsize)
    end
    for (y, istp) in zip(pos_y_legend, val_istp_legend)
        text!(ax, 0.16, y; text=string(istp),
            align=(:left, :center), font=DUALMOT_FONT,
            fontsize=legend_fontsize)
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
    save(joinpath(path_result_legends,
        "result_legend_162-164_$(name).pdf"), fig)
end
for version in legend_versions
    draw_result_legend(version.loadcfgs, version.name)
end

function draw_result_isotope_loadcfg_legend()
    isotopes = Symbol.(string.(160:164))
    loadcfgs = (:DIS, :DDM)
    pos_x = [pos_x_start_legend + (idx - 1) * pos_x_step_legend
        for idx in eachindex(isotopes)]
    x_max = last(pos_x) + 0.75
    fig_width = round(Int, (106 - 4) * x_max / 6 + 4)
    fig = Figure(size=(fig_width, fig_height_legend),
        fontsize=legend_fontsize, figure_padding=2)
    ax = Axis(fig[1, 1];
        width=fig_width - 4,
        height=fig_height_legend - 4,
        halign=0,
        valign=0,
        tellwidth=false,
        tellheight=false,
        limits=(0.0, x_max, axis_limits_y_legend...),
        xticksvisible=false,
        xticklabelsvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
        xgridvisible=false,
        ygridvisible=false,
        bottomspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        leftspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        topspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        rightspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        spinewidth=0.75,
        backgroundcolor=RGBAf(1, 1, 1, 0.86),
    )
    for (x, istp) in zip(pos_x, isotopes)
        text!(ax, x, 1.90; text=string(istp),
            align=(:center, :center), font=DUALMOT_FONT,
            fontsize=legend_fontsize)
    end
    for (y, loadcfg) in zip(pos_y_legend, loadcfgs)
        text!(ax, 0.16, y; text=string(loadcfg),
            align=(:left, :center), font=DUALMOT_FONT,
            fontsize=legend_fontsize)
        for (x, istp) in zip(pos_x, isotopes)
            style = dualmot_curve_style((; loadcfg, istp))
            scatter!(ax, [x], [y];
                color=style.markercolor,
                marker=style.marker,
                markersize=style.markersize,
                strokecolor=style.strokecolor,
                strokewidth=style.strokewidth)
        end
    end
    save(joinpath(path_result_legends,
        "result_legend_160-164_DIS-DDM.pdf"), fig)
end

function draw_result_isotope_ratio_legend()
    isotopes = Symbol.(string.(160:164))
    pos_x = [pos_x_start_legend + (idx - 1) * pos_x_step_legend
        for idx in eachindex(isotopes)]
    x_max = last(pos_x) + 0.75
    fig_width = round(Int, (106 - 4) * x_max / 6 + 4)
    fig = Figure(size=(fig_width, fig_height_legend),
        fontsize=legend_fontsize, figure_padding=2)
    ax = Axis(fig[1, 1];
        width=fig_width - 4,
        height=fig_height_legend - 4,
        halign=0,
        valign=0,
        tellwidth=false,
        tellheight=false,
        limits=(0.0, x_max, axis_limits_y_legend...),
        xticksvisible=false,
        xticklabelsvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
        xgridvisible=false,
        ygridvisible=false,
        bottomspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        leftspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        topspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        rightspinecolor=RGBAf(0.25, 0.25, 0.25, 0.45),
        spinewidth=0.75,
        backgroundcolor=RGBAf(1, 1, 1, 0.86),
    )
    y_marker = 1.20
    for (x, istp) in zip(pos_x, isotopes)
        text!(ax, x, 1.90; text=string(istp),
            align=(:center, :center), font=DUALMOT_FONT,
            fontsize=legend_fontsize)
        style = dualmot_curve_style((; loadcfg=:DDM, istp))
        scatter!(ax, [x], [y_marker];
            color=style.markercolor,
            marker=:hexagon,
            markersize=style.markersize,
            strokecolor=style.strokecolor,
            strokewidth=style.strokewidth)
    end
    save(joinpath(path_result_legends,
        "result_legend_160-164_hexagon.pdf"), fig)
end

draw_result_isotope_loadcfg_legend()
draw_result_isotope_ratio_legend()
println("Saved five Result legend PDFs to $path_result_legends")
