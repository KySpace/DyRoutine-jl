function pair_istp_position(idx_pair::Integer, pair::AbstractString, istp::Symbol)
    isotopes = Symbol.(split(pair, "-"))
    length(isotopes) == 2 || throw(ArgumentError("invalid isotope pair: $pair"))
    direction = istp == isotopes[1] ? -1 : istp == isotopes[2] ? 1 :
        throw(ArgumentError("$istp is not part of isotope pair $pair"))
    idx_pair + direction * pair_istp_offset
end

marker_style(style::NamedTuple; markersize::Real=style.markersize) = (
    color=style.markercolor,
    marker=style.marker,
    markersize,
    strokecolor=style.strokecolor,
    strokewidth=style.strokewidth,
)

function draw_pair_spans!(ax::Axis)
    for idx_pair in eachindex(val_pair)
        vspan!(ax, idx_pair - 0.5, idx_pair + 0.5; color=:white) |>
            span -> translate!(span, 0, 0, -100)
    end
    nothing
end

function style_key_axis(fig::Figure; width::Real, limits)
    border_color = RGBAf(0.25, 0.25, 0.25, 0.45)
    Axis(fig[1, 1];
        width,
        height=58,
        halign=0.05,
        valign=0.10,
        tellwidth=false,
        tellheight=false,
        limits,
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
end

function draw_number_style_key!(fig::Figure, loadcfgs::Tuple)
    val_istp_key = Symbol.(string.(160:164))
    pos_y = collect(5:-1:1)
    pos_x = 1.10 .+ 0.34 .* (0:length(loadcfgs)-1)
    ax_key = style_key_axis(fig;
        width=length(loadcfgs) > 2 ? 76 : 56,
        limits=(0.65, max(1.70, last(pos_x) + 0.26), 0.5, 6.2),
    )
    text!(ax_key, mean(pos_x), 5.8;
        text=join(string.(loadcfgs), "  "),
        align=(:center, :center), fontsize=22 / 3)
    for (idx_istp, istp) in enumerate(val_istp_key)
        text!(ax_key, 0.72, pos_y[idx_istp]; text=string(istp),
            align=(:left, :center), fontsize=22 / 3)
    end
    for (idx_loadcfg, loadcfg) in enumerate(loadcfgs),
        (idx_istp, istp) in enumerate(val_istp_key)
        style = dualmot_curve_style((; loadcfg, istp))
        scatter!(ax_key, [pos_x[idx_loadcfg]], [pos_y[idx_istp]];
            marker_style(style; markersize=6)...)
    end
    ax_key
end

function draw_ratio_style_key!(fig::Figure)
    val_istp_key = Symbol.(string.(160:164))
    pos_y = collect(5:-1:1)
    ax_key = style_key_axis(fig;
        width=28,
        limits=(0.65, 1.25, 0.5, 5.5),
    )
    for (idx_istp, istp) in enumerate(val_istp_key)
        text!(ax_key, 0.7, pos_y[idx_istp]; text=string(istp),
            align=(:left, :center), fontsize=22 / 3)
        style = dualmot_curve_style((; loadcfg=:DDM, istp))
        scatter!(ax_key, [1.07], [pos_y[idx_istp]];
            merge(marker_style(style; markersize=6), (; marker=:hexagon))...)
    end
    ax_key
end
