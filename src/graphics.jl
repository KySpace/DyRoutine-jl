using CairoMakie: CairoMakie, Figure, Axis, axislegend, save, band!, lines!, scatter!, RGBf, RGBAf
using Colors: Oklch

function write_number_plot(
    path_plot::AbstractString,
    val_t_hold,
    val_number::AbstractMatrix,
    err_number::AbstractMatrix,
    val_istp,
)
    size(val_number) == size(err_number) || throw(
        ArgumentError("val_number and err_number must have the same size."),
    )
    size(val_number, 1) == length(val_t_hold) || throw(
        ArgumentError("t_hold axis length does not match number statistics."),
    )
    size(val_number, 2) == length(val_istp) || throw(
        ArgumentError("istp axis length does not match number statistics."),
    )

    pos_istp_blue = findfirst(==("162"), string.(val_istp))
    pos_istp_red = findfirst(==("164"), string.(val_istp))
    !isnothing(pos_istp_blue) || throw(ArgumentError("Expected istp value 162 in val_istp."))
    !isnothing(pos_istp_red) || throw(ArgumentError("Expected istp value 164 in val_istp."))

    color_blue_line = RGBAf(Oklch(0.62, 0.18, 255.0), 0.95)
    color_blue_fill = RGBAf(Oklch(0.62, 0.18, 255.0), 0.22)
    color_red_line = RGBAf(Oklch(0.62, 0.18, 25.0), 0.95)
    color_red_fill = RGBAf(Oklch(0.62, 0.18, 25.0), 0.22)

    fig = Figure(
        size=(920, 620),
        backgroundcolor=RGBf(0.98, 0.977, 0.965),
    )
    ax = Axis(
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
        "istp = 162",
    )
    plot_number_series!(
        ax,
        Float64.(collect(val_t_hold)),
        vec(val_number[:, pos_istp_red]),
        vec(err_number[:, pos_istp_red]),
        color_red_line,
        color_red_fill,
        "istp = 164",
    )

    axislegend(ax; position=:rt, framevisible=true)
    save(path_plot, fig)
    return path_plot
end

function plot_number_series!(ax, val_t, val_mean, val_err, color_line, color_fill, label)
    band!(ax, val_t, val_mean .- val_err, val_mean .+ val_err; color=color_fill)
    lines!(ax, val_t, val_mean; color=color_line, linewidth=2.6, label=label)
    scatter!(ax, val_t, val_mean; color=color_line, markersize=9)
end
