function write_number_plot(
    path_plot::AbstractString,
    val_t_hold,
    val_number::AbstractMatrix,
    err_number::AbstractMatrix,
    val_istp,
)
    Base.find_package("CairoMakie") === nothing && throw(
        ArgumentError("CairoMakie is not installed in the current Julia environment."),
    )

    @eval using CairoMakie

    size(val_number) == size(err_number) || throw(
        ArgumentError("val_number and err_number must have the same size."),
    )
    size(val_number, 1) == length(val_t_hold) || throw(
        ArgumentError("t_hold axis length does not match number statistics."),
    )
    size(val_number, 2) == length(val_istp) || throw(
        ArgumentError("istp axis length does not match number statistics."),
    )

    pos_istp_blue = findfirst(==(5), val_istp)
    pos_istp_red = findfirst(==(0), val_istp)
    !isnothing(pos_istp_blue) || throw(ArgumentError("Expected istp value 5 in val_istp."))
    !isnothing(pos_istp_red) || throw(ArgumentError("Expected istp value 0 in val_istp."))

    color_blue_line = color_from_oklch(0.62, 0.18, 255.0, 0.95)
    color_blue_fill = color_from_oklch(0.62, 0.18, 255.0, 0.22)
    color_red_line = color_from_oklch(0.62, 0.18, 25.0, 0.95)
    color_red_fill = color_from_oklch(0.62, 0.18, 25.0, 0.22)

    fig = CairoMakie.Figure(size=(920, 620), backgroundcolor=RGBf(0.98, 0.977, 0.965))
    ax = CairoMakie.Axis(
        fig[1, 1];
        title="Mean Number vs t_hold",
        xlabel="t_hold",
        ylabel="number",
        xgridvisible=true,
        ygridvisible=true,
    )

    plot_number_series!(
        ax,
        Float64.(collect(val_t_hold)),
        vec(val_number[:, pos_istp_blue]),
        vec(err_number[:, pos_istp_blue]),
        color_blue_line,
        color_blue_fill,
        "istp = 5",
    )
    plot_number_series!(
        ax,
        Float64.(collect(val_t_hold)),
        vec(val_number[:, pos_istp_red]),
        vec(err_number[:, pos_istp_red]),
        color_red_line,
        color_red_fill,
        "istp = 0",
    )

    CairoMakie.axislegend(ax; position=:rt, framevisible=true)
    CairoMakie.save(path_plot, fig)
    return path_plot
end

function plot_number_series!(ax, val_t, val_mean, val_err, color_line, color_fill, label)
    CairoMakie.band!(ax, val_t, val_mean .- val_err, val_mean .+ val_err; color=color_fill)
    CairoMakie.lines!(ax, val_t, val_mean; color=color_line, linewidth=2.6, label=label)
    CairoMakie.scatter!(ax, val_t, val_mean; color=color_line, markersize=9)
end

function color_from_oklch(val_l, val_c, val_h_deg, val_alpha)
    val_r, val_g, val_b = oklch_to_srgb(val_l, val_c, val_h_deg)
    return RGBAf(val_r, val_g, val_b, val_alpha)
end

function oklch_to_srgb(val_l, val_c, val_h_deg)
    val_h_rad = deg2rad(val_h_deg)
    val_a = val_c * cos(val_h_rad)
    val_b = val_c * sin(val_h_rad)

    val_l_ = val_l + 0.3963377774 * val_a + 0.2158037573 * val_b
    val_m_ = val_l - 0.1055613458 * val_a - 0.0638541728 * val_b
    val_s_ = val_l - 0.0894841775 * val_a - 1.2914855480 * val_b

    val_l_lin = val_l_^3
    val_m_lin = val_m_^3
    val_s_lin = val_s_^3

    val_r_lin = 4.0767416621 * val_l_lin - 3.3077115913 * val_m_lin + 0.2309699292 * val_s_lin
    val_g_lin = -1.2684380046 * val_l_lin + 2.6097574011 * val_m_lin - 0.3413193965 * val_s_lin
    val_b_lin = -0.0041960863 * val_l_lin - 0.7034186147 * val_m_lin + 1.7076147010 * val_s_lin

    return (
        srgb_from_linear(val_r_lin),
        srgb_from_linear(val_g_lin),
        srgb_from_linear(val_b_lin),
    )
end

function srgb_from_linear(val_rgb)
    val_clamped = clamp(val_rgb, 0.0, 1.0)
    return val_clamped <= 0.0031308 ? 12.92 * val_clamped : 1.055 * val_clamped^(1 / 2.4) - 0.055
end
