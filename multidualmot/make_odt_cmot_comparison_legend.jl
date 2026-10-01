include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_output = joinpath(path_root, "Isotope pair comparison")
mkpath(path_output)
val_istp_legend = Symbol.(string.(160:164))

function comparison_legend_axis!(fig::Figure, slot; width::Real, height::Real, limits)
    border_color = RGBAf(0.25, 0.25, 0.25, 0.45)
    Axis(slot;
        width, height,
        halign=0.05, valign=0.10,
        tellwidth=false, tellheight=false,
        limits,
        xticksvisible=false, xticklabelsvisible=false,
        yticksvisible=false, yticklabelsvisible=false,
        xgridvisible=false, ygridvisible=false,
        bottomspinecolor=border_color, leftspinecolor=border_color,
        topspinecolor=border_color, rightspinecolor=border_color,
        spinewidth=0.75, backgroundcolor=RGBAf(1, 1, 1, 0.86),
    )
end

function save_comparison_legend(fig::Figure, filename::AbstractString)
    for format in ("svg", "png")
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_output, "$filename.$format"), fig; save_options...)
    end
end

function draw_number_legend()
    fig = Figure(size=(90, 62), fontsize=6, figure_padding=1)
    ax = comparison_legend_axis!(fig, fig[1,1];
        width=88, height=60, limits=(0.15, 7.2, 0.5, 7.2))
    x_positions = (1.9, 3.2, 5.0, 6.3)
    text!(ax, 2.55, 6.65; text="DDM", align=(:center,:center), fontsize=5.5)
    text!(ax, 5.65, 6.65; text="DIS", align=(:center,:center), fontsize=5.5)
    for (idx, x) in enumerate(x_positions)
        text!(ax, x, 5.95; text=idx in (1, 3) ? "ODT" : "CMOT",
            align=(:center,:center), fontsize=4)
    end
    for (idx_istp, istp) in enumerate(val_istp_legend)
        y = 5 - idx_istp + 1
        text!(ax, 0.36, y; text=string(istp), align=(:left,:center), fontsize=6)
        for (idx_loadcfg, loadcfg) in enumerate((:DDM, :DIS)),
            (idx_source, open_marker) in enumerate((false, true))
            x = x_positions[2(idx_loadcfg - 1) + idx_source]
            style = dualmot_curve_style((; loadcfg, istp))
            scatter!(ax, [x], [y];
                color=open_marker ? :transparent : style.markercolor,
                marker=style.marker, markersize=5,
                strokecolor=style.strokecolor, strokewidth=style.strokewidth)
        end
    end
    fig
end

function draw_ratio_legend()
    fig = Figure(size=(48, 60), fontsize=7, figure_padding=1)
    ax = comparison_legend_axis!(fig, fig[1,1];
        width=46, height=58, limits=(0.5, 3.5, 0.5, 6.2))
    x_positions = (1.5, 2.7)
    text!(ax, x_positions[1], 5.8; text="DDM", align=(:center,:center), fontsize=7)
    text!(ax, x_positions[2], 5.8; text="DIS", align=(:center,:center), fontsize=7)
    for (idx_istp, istp) in enumerate(val_istp_legend),
        (idx_loadcfg, loadcfg) in enumerate((:DDM, :DIS))
        y = 5 - idx_istp + 1
        if idx_loadcfg == 1
            text!(ax, 0.62, y; text=string(istp), align=(:left,:center), fontsize=6)
        end
        style = dualmot_curve_style((; loadcfg, istp))
        scatter!(ax, [x_positions[idx_loadcfg]], [y];
            color=style.markercolor, marker=style.marker, markersize=5,
            strokecolor=style.strokecolor, strokewidth=style.strokewidth)
    end
    fig
end

fig_numbers_legend = draw_number_legend()
fig_ratios_legend = draw_ratio_legend()
save_comparison_legend(fig_numbers_legend, "[ODT.CMOT.comparison].[legend].[numbers]")
save_comparison_legend(fig_ratios_legend, "[ODT.CMOT.comparison].[legend].[ratios]")
println("Saved number and ratio legend figures to $path_output")
