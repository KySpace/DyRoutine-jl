"""Export the standalone legend boxes for added isotope-pair curve figures."""

path_result_curve_legends = joinpath(path_root, "Result")
mkpath(path_result_curve_legends)
legend_curve_fontsize = 6 / 0.75
legend_curve_font = "NewComputerModern Math"
legend_curve_x_start = 1.25
legend_curve_x_step = 0.95
legend_curve_y_step = 0.45

function draw_selected_curve_legend(isotopes::Tuple, loadcfgs::Tuple,
    filename::AbstractString)
    pos_x = [legend_curve_x_start + (idx - 1) * legend_curve_x_step
        for idx in eachindex(loadcfgs)]
    pos_y = [1.30 - (idx - 1) * legend_curve_y_step
        for idx in eachindex(isotopes)]
    x_max = last(pos_x) + 0.7
    y_min = 0.25
    y_max = 2.05
    legend_box_width = round(Int, max(74, 32 + x_max * 13))
    fig_height = round(Int, max(32, 14 + (y_max - y_min) * 16))
    fig = Figure(size=(fig_width, fig_height), fontsize=legend_curve_fontsize,
        figure_padding=2)
    border_color = RGBAf(0.25, 0.25, 0.25, 0.45)
    ax = Axis(fig[1, 1]; width=legend_box_width - 4, height=fig_height - 4,
        halign=:center, valign=:center, tellwidth=false, tellheight=false,
        limits=(0.0, x_max, y_min, y_max),
        xticksvisible=false, xticklabelsvisible=false,
        yticksvisible=false, yticklabelsvisible=false,
        xgridvisible=false, ygridvisible=false,
        bottomspinecolor=border_color, leftspinecolor=border_color,
        topspinecolor=border_color, rightspinecolor=border_color,
        spinewidth=0.75, backgroundcolor=RGBAf(1, 1, 1, 0.86))
    for (x, loadcfg) in zip(pos_x, loadcfgs)
        text!(ax, x, 1.75; text=string(loadcfg),
            align=(:center, :center), font=legend_curve_font,
            fontsize=legend_curve_fontsize)
    end
    for (y, isotope) in zip(pos_y, isotopes)
        text!(ax, 0.16, y; text=string(isotope),
            align=(:left, :center), font=legend_curve_font,
            fontsize=legend_curve_fontsize)
        for (x, loadcfg) in zip(pos_x, loadcfgs)
            style = dualmot_curve_style((; loadcfg, istp=isotope))
            scatter!(ax, [x], [y]; color=style.markercolor,
                marker=style.marker, markersize=style.markersize,
                strokecolor=style.strokecolor, strokewidth=style.strokewidth)
        end
    end
    save(joinpath(path_result_curve_legends, filename), fig)
    fig
end

draw_selected_curve_legend((Symbol("162"),), (:SIS, :SDS),
    "result_legend_162_SIS-SDS.pdf")
draw_selected_curve_legend((Symbol("161"), Symbol("163")), (:SIS, :SDS),
    "result_legend_161-163_SIS-SDS.pdf")
draw_selected_curve_legend((Symbol("162"), Symbol("163")), (:SCS, :DCS),
    "result_legend_162-163_SCS-DCS.pdf")
draw_selected_curve_legend((Symbol("161"), Symbol("163")), (:SCS, :DCS),
    "result_legend_161-163_SCS-DCS.pdf")
draw_selected_curve_legend((Symbol("163"), Symbol("164")), (:DCS, :SCS),
    "result_legend_163-164_DCS-SCS.pdf")
draw_selected_curve_legend((Symbol("163"),), (:SIS, :SDS),
    "result_legend_163_SIS-SDS.pdf")
draw_selected_curve_legend((Symbol("162"), Symbol("163")), (:DIS, :DDM),
    "result_legend_162-163_DIS-DDM.pdf")
draw_selected_curve_legend((Symbol("163"), Symbol("164")), (:DIS, :DDM),
    "result_legend_163-164_DIS-DDM.pdf")
draw_selected_curve_legend((Symbol("161"), Symbol("163")), (:DIS, :DDM),
    "result_legend_161-163_DIS-DDM.pdf")

println("Saved nine selected-curve Result legends to $path_result_curve_legends")
