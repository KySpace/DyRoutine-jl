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

function set_axis_full(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer}, panel_setter::Function)
    CairoMakie.activate!()
    CairoMakie.activate!()
    fig = Figure()
    length(n_dim_vars) == 3 || throw(ArgumentError("n_dim_vars must be a 3-tuple"))
    axs_solo = Array{Dict}(undef, n_dim_vars)
    axs_stacked = Array{Dict}(undef, n_dim_vars[2:end])
    for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
        print("\rbuilding solo axis for rep $r, $t")
        gl = GridLayout()
        # fig[1, 1][t, (r-1)*n_dim_vars[3]+i] = gl
        fig[1, 3*(i-1)+1][t, r] = gl
        axs_solo[r, t, i] = panel_setter(gl)
    end
    for t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
        print("\rbuilding stack axis for $t")
        gl = GridLayout()
        fig[1, 3*(i-1)+2][t, 1] = gl
        axs_stacked[t, i] = panel_setter(gl)
    end
    colsize!(fig.layout, 3, Fixed(2))
    return fig, axs_solo, axs_stacked
end

function clean_gridlayout!(gl::GridLayout)
    for obj in contents(gl)
        obj isa Axis && delete!(obj)
    end
    trim!(gl)
end

function set_panel_solo_essn_2d!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_dens = Axis(gl[1:2, 1])
    ax_modl = Axis(gl[1, 2])
    ax_prfl_ft = Axis(gl[2, 2])
    colsize!(gl, 1, Fixed(200))
    colsize!(gl, 2, Fixed(300))
    rowsize!(gl, 1, Fixed(200))
    rowsize!(gl, 2, Fixed(100))
    return Dict("dens" => ax_dens, "modl" => ax_modl, "prfl_ft" => ax_prfl_ft)
end

function set_panel_solo_modl!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_modl = Axis(gl[1, 2])
    ax_prfl_ft_upright = Axis(gl[1, 3])
    ax_prfl_ft_sideway = Axis(gl[1, 1])
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(160))
    colsize!(gl, 3, Fixed(160))
    rowsize!(gl, 1, Fixed(80))
    return Dict("modl" => ax_modl, "upright" => ax_prfl_ft_upright, "sideway" => ax_prfl_ft_sideway)
end

function set_axis_pca_4x4!()
    fig = Figure()
    axs_mode = Array{Dict}(undef, 16)
    for r in 1:4, c in 1:4
        gl = GridLayout()
        fig[r, c] = gl
        axs_mode[(r-1)*4+c] = set_panel_pca_duet!(gl)
    end
    return fig, axs_mode
end

function set_axis_pca_dual_4x2!()
    fig = Figure()
    axs_mode = Array{Dict}(undef, (2, 8))
    for r in 1:2, c in 1:4
        gl = GridLayout()
        fig[1, 1][r, c] = gl
        axs_mode[1, (r-1)*4+c] = set_panel_pca_solo!(gl)
        Box(fig[2, 1], color=:black)
        gl = GridLayout()
        fig[3, 1][r, c] = gl
        axs_mode[2, (r-1)*4+c] = set_panel_pca_solo!(gl)
    end
    rowsize!(fig.layout, 2, Fixed(2))
    return fig, axs_mode
end

function set_panel_pca_duet!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_l = Axis(gl[1:2, 1])
    ax_r = Axis(gl[1:2, 2])
    ax_evol = Axis(gl[1, 3])
    ax_freq = Axis(gl[2, 3])
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(100))
    colsize!(gl, 3, Fixed(200))
    rowsize!(gl, 1, Fixed(150))
    rowsize!(gl, 2, Fixed(150))
    return Dict("l" => ax_l, "r" => ax_r, "evol" => ax_evol, "freq" => ax_freq)
end

function set_panel_pca_solo!(gl::GridLayout)
    for obj in contents(gl)
        obj isa Axis && delete!(obj)
    end
    trim!(gl)
    ax_mode = Axis(gl[1:2, 1])
    ax_evol = Axis(gl[1, 2])
    ax_freq = Axis(gl[2, 2])
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(200))
    rowsize!(gl, 1, Fixed(120))
    rowsize!(gl, 2, Fixed(120))
    return Dict("mode" => ax_mode, "evol" => ax_evol, "freq" => ax_freq)
end
