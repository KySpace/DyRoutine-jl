using CairoMakie, GLMakie
using Colors: Oklch
using LaTeXStrings

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

function set_axis_sidepeak!(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer}, panel_setter::Function)
    length(n_dim_vars) == 3 || throw(ArgumentError("n_dim_vars must be a 3-tuple"))
    fig = Figure()
    axs_repeats = Array{Dict}(undef, n_dim_vars[1])
    for r in 1:n_dim_vars[1]
        fig[1, r] = Label(fig, text="repeat $r"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
        print("\rbuilding axis for side peak trend for repeat $r")
        gl = GridLayout()
        fig[2, r] = gl
        axs_repeats[r] = panel_setter(gl, r)
    end
    fig[2, n_dim_vars[1]+1] |> Box
    gl = GridLayout()
    fig[1, n_dim_vars[1]+2] = Label(fig, text="Processed after stacked"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+2] = gl
    axs_stacked = panel_setter(gl, 1)
    gl = GridLayout()
    fig[1, n_dim_vars[1]+3] = Label(fig, text="Reps overlayed"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+3] = gl
    axs_all = panel_setter(gl, 1)
    return fig, Dict("repeats" => axs_repeats, "stacked" => axs_stacked, "all" => axs_all)
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
    ax_modl = Axis(gl[1, 3])
    ax_dens = Axis(gl[1, 2])
    ax_prfl_ft_upright = Axis(gl[1, 4])
    ax_prfl_ft_sideway = Axis(gl[1, 1])
    colsize!(gl, 1, Fixed(200))
    colsize!(gl, 2, Fixed(80))
    colsize!(gl, 3, Fixed(240))
    colsize!(gl, 4, Fixed(240))
    colgap!(gl, 5)
    rowsize!(gl, 1, Fixed(160))
    return Dict("dens" => ax_dens, "modl" => ax_modl, "upright" => ax_prfl_ft_upright, "sideway" => ax_prfl_ft_sideway)
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
    gl |> clean_gridlayout!
    ax_mode = Axis(gl[1:2, 1])
    ax_evol = Axis(gl[1, 2])
    ax_freq = Axis(gl[2, 2])
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(200))
    rowsize!(gl, 1, Fixed(120))
    rowsize!(gl, 2, Fixed(120))
    return Dict("mode" => ax_mode, "evol" => ax_evol, "freq" => ax_freq)
end

function set_panel_trend_sidepeak_nvlp!(gl::GridLayout, col::Int)
    gl |> clean_gridlayout!
    w, h = (400, 200)
    ax_evol_weight = Axis(gl[1, 1]; width=w, height=h, ylabel="side peak \nweight")
    ax_evol_height = Axis(gl[2, 1]; width=w, height=h, ylabel="side peak \nheight")
    ax_evol_width = Axis(gl[3, 1]; width=w, height=h, ylabel="side peak \nwidth (μm⁻¹)")
    ax_evol_wavenum = Axis(gl[4, 1]; width=w, height=h, ylabel="side peak \nwavenum (μm⁻¹)")
    ax_evol_sizes = Axis(gl[5, 1]; width=w, height=h, ylabel="envelope size (μm)")
    ax_freq_weight = Axis(gl[1, 2]; width=w, height=h)
    ax_freq_height = Axis(gl[2, 2]; width=w, height=h)
    ax_freq_width = Axis(gl[3, 2]; width=w, height=h)
    ax_freq_wavenum = Axis(gl[4, 2]; width=w, height=h)
    ax_freq_sizes = Axis(gl[5, 2]; width=w, height=h)
    rowgap!(gl, 4)
    rowgap!(gl, 4)
    dict_axs = Dict(
        "evol-weight" => ax_evol_weight,
        "evol-height" => ax_evol_height,
        "evol-width" => ax_evol_width,
        "evol-wavenum" => ax_evol_wavenum,
        "evol-sizes" => ax_evol_sizes,
        "freq-weight" => ax_freq_weight,
        "freq-height" => ax_freq_height,
        "freq-width" => ax_freq_width,
        "freq-wavenum" => ax_freq_wavenum,
        "freq-sizes" => ax_freq_sizes,
    )
    for ax in values(dict_axs)
        hideydecorations!(ax; label=true, ticklabels=false, ticks=false, grid=false, minorticks=false, minorgrid=false)
        hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
        if col == 1
            ax.ylabelvisible = true
        end
    end
    ax_evol_sizes.xticklabelsvisible = true
    ax_freq_sizes.xticklabelsvisible = true
    ax_evol_sizes.xlabelvisible = true
    ax_freq_sizes.xlabelvisible = true
    ax_evol_sizes.xlabel = "t hold (ms)"
    ax_freq_sizes.xlabel = "freq (Hz)"
    return dict_axs
end

function draw_rotated_ellipse!(
    ax,
    center::Tuple{<:Real,<:Real},
    radii::Tuple{<:Real,<:Real},
    angle::Real;
    n=200,
    kwargs...
)
    x0, y0 = center
    rx, ry = radii
    φ = angle

    θ = range(0, 2π, length=n)

    x = x0 .+ rx .* cos.(θ) .* cos(φ) .- ry .* sin.(θ) .* sin(φ)
    y = y0 .+ rx .* cos.(θ) .* sin(φ) .+ ry .* sin.(θ) .* cos(φ)

    obj = lines!(ax, x, y; kwargs...)
    return obj
end

function draw_rotated_ellipse_corners!(
    ax,
    center::Tuple{<:Real,<:Real},
    radii::Tuple{<:Real,<:Real},
    angle::Real;
    kwargs...
)
    x0, y0 = center
    rx, ry = radii
    c = cos(angle)
    s = sin(angle)
    rot = (x, y) -> (x * c - y * s, x * s + y * c)

    l = minimum([rx, ry]) / 2
    vtx_x = []
    vtx_y = []
    for u in [+1, -1], v in [+1 -1]
        rxuv, ryuv = rot(u * rx, v * ry)
        rxuvly, ryuvly = rot(u * rx, v * (ry - l))
        rxuvlx, ryuvlx = rot(u * (rx - l), v * ry)
        push!(vtx_x, x0 + rxuvlx)
        push!(vtx_x, x0 + rxuv)
        push!(vtx_x, x0 + rxuvly)
        push!(vtx_x, NaN)
        push!(vtx_y, y0 + ryuvlx)
        push!(vtx_y, y0 + ryuv)
        push!(vtx_y, y0 + ryuvly)
        push!(vtx_y, NaN)
    end
    obj = lines!(ax, vtx_x, vtx_y; kwargs...)
    return obj
end
