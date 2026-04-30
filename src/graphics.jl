using CairoMakie, GLMakie
using Colors: Oklch

hue_theme_istp = Dict("162" => 11.0, "164" => 250.0)
GLMakie.activate!()

function gen_clrmap_solo(hue)
    return [Oklch(1 - t, 0.24 * t, hue) |> c -> RGBAf(c) for t in range(0, 1; length=256)]
end

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
    clr_bar = Oklch(0.62, 0.18, hue) |> c -> RGBAf(c, 0.40)

    errorbars!(ax, val_t, val_mean, val_err; color=clr_bar)
    lines!(ax, val_t, val_mean; color=clr_bar, linewidth=2.6, label=val_istp)
    scatter!(ax, val_t, val_mean; color=clr_bar, markersize=9)
    # save(path_plot, fig; backend=CairoMakie)
    # return path_plot
end

function set_axis_full(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer})
    CairoMakie.activate!()
    n_row = n_dim_vars[2]
    n_col = n_dim_vars[1] * n_dim_vars[3]
    CairoMakie.activate!()
    fig = Figure()
    axs_solo = Array{Dict}(undef, n_dim_vars)
    for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
        print("\rbuilding axis for rep $i, $t")
        gl = GridLayout()
        # fig[1, 1][t, (r-1)*n_dim_vars[3]+i] = gl
        fig[1, 1][t, r+(i-1)*n_dim_vars[1]] = gl
        axs_solo[r, t, i] = set_panel_solo_essn_2d!(gl)
    end
    return fig, axs_solo
end

function set_panel_solo_essn_2d!(gl::GridLayout)
    for obj in contents(gl)
        obj isa Axis && delete!(obj)
    end
    trim!(gl)
    ax_dens = Axis(gl[1:2, 1])
    ax_modl = Axis(gl[1, 2])
    ax_prfl_ft = Axis(gl[2, 2])
    colsize!(gl, 1, Fixed(200))
    colsize!(gl, 2, Fixed(300))
    rowsize!(gl, 1, Fixed(200))
    rowsize!(gl, 2, Fixed(100))
    return Dict("dens" => ax_dens, "modl" => ax_modl, "prfl_ft" => ax_prfl_ft)
end
