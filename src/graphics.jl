using CairoMakie, GLMakie
using Colors: Oklch

hue_theme_istp = Dict("162" => 11.0, "164" => 250.0)
GLMakie.activate!()

function set_axis!(title::String)
    fig = Figure(
        size=(920, 620),
    )
    ax = Axis(
        fig[1, 1];
        title=title,
        xgridvisible=true,
        ygridvisible=true,
    )
    return fig, ax
end

function plot_num_stat_evo!(
    ax::Axis,
    val_t::AbstractVector,
    stat_number::AbstractVector,
    val_istp,
)
    hue = hue_theme_istp[val_istp]
    val_mean = map(s -> s[1], stat_number)
    val_err = map(s -> s[2], stat_number)
    clr_line = Oklch(0.62, 0.18, hue) |> c -> RGBAf(c, 0.95)
    clr_bar  = Oklch(0.62, 0.18, hue) |> c -> RGBAf(c, 0.40)

    errorbars!(ax, val_t, val_mean, val_err; color=clr_bar)
    lines!(ax, val_t, val_mean; color=clr_bar, linewidth=2.6, label=val_istp)
    scatter!(ax, val_t, val_mean; color=clr_bar, markersize=9)
    # save(path_plot, fig; backend=CairoMakie)
    # return path_plot
end
