using CairoMakie: extract_attributes!
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
        print("\r\033[2Kbuilding solo axis for rep $r, $t")
        gl = GridLayout()
        # fig[1, 1][t, (r-1)*n_dim_vars[3]+i] = gl
        fig[1, 3*(i-1)+1][t, r] = gl
        axs_solo[r, t, i] = panel_setter(gl)
    end
    for t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
        print("\r\033[2Kbuilding stack axis for $t")
        gl = GridLayout()
        fig[1, 3*(i-1)+2][t, 1] = gl
        axs_stacked[t, i] = panel_setter(gl)
    end
    colsize!(fig.layout, 3, Fixed(2))
    return fig, axs_solo, axs_stacked
end

function set_axis_sidepeak_nvlp!(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer}, panel_setter::Function, runinfo)
    length(n_dim_vars) == 3 || throw(ArgumentError("n_dim_vars must be a 3-tuple"))
    fig = Figure()
    axs_repeats = Array{Dict}(undef, n_dim_vars[1])
    fig[0, 1] = Label(fig, text="$(runinfo.date) $(@sprintf("run%02d", runinfo.runid)) IB=$(@sprintf("%.3f", runinfo.IB))A $(runinfo.tag_head)"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    for r in 1:n_dim_vars[1]
        fig[1, r] = Label(fig, text="repeat $r"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
        print("\r\033[2K\rbuilding axes for side peak trend for repeat $r")
        gl = GridLayout()
        fig[2, r] = gl
        axs_repeats[r] = panel_setter(gl, r)
    end
    println("\r\033[2K\raxes built for trends.")
    fig[2, n_dim_vars[1]+1] |> Box
    colsize!(fig.layout, n_dim_vars[1] + 1, Fixed(2))
    gl = GridLayout()
    fig[1, n_dim_vars[1]+2] = Label(fig, text="Processed after stacked"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+2] = gl
    axs_stacked = panel_setter(gl, 1)
    gl = GridLayout()
    fig[1, n_dim_vars[1]+3] = Label(fig, text="Reps overlayed"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+3] = gl
    axs_all = panel_setter(gl, 1; extra=true)
    return fig, Dict("repeats" => axs_repeats, "stacked" => axs_stacked, "all" => axs_all)
end

function clean_gridlayout!(gl::GridLayout)
    for obj in contents(gl)
        obj isa Axis && delete!(obj)
    end
    trim!(gl)
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
