using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using Colors: Oklch
using LaTeXStrings

hue_theme_istp = Dict("162" => 11.0, "164" => 250.0)
GLMakie.activate!()

function gen_clrmap_solo(hue; thres_alpha=0.0, alpha_base=1.0)
    clrmap = [
        begin
            alpha = thres_alpha <= 0 || abs(t) > thres_alpha ?
                    1.0 :
                    (abs(t) / thres_alpha * (1 - alpha_base) + alpha_base)
            Oklch(1 - 0.8 * t, 0.24 * t, hue) |> c -> RGBAf(c, alpha)
        end
        for t in range(0, 1; length=256)]
    return clrmap
end

function gen_clrmap_posneg_nonlin(hue_pos, hue_neg; thres_alpha=0.6, alpha_base=0.2)
    return [
        begin
            alpha = abs(t) > thres_alpha ? 1.0 : (abs(t) / thres_alpha * (1 - alpha_base) + alpha_base)
            Oklch(1 - 0.6 * abs(t), 0.4 * abs2(t), t > 0 ? hue_pos : hue_neg) |> c -> RGBAf(c, alpha)
        end
        for t in range(-1, 1; length=256)
    ]
end

function clear_axes!(axs)
    for obj in axs
        obj isa Axis && empty!(obj)
    end
end

function set_axis!(title::String; kwargs...)
    fig = Figure(
        size=(920, 620),
    )
    ax = Axis(
        fig[1, 1];
        title=title,
        xgridvisible=true,
        ygridvisible=true,
        kwargs...
    )
    return fig, ax
end

function clean_gridlayout!(gl::GridLayout)
    for obj in reverse(contents(gl))
        obj isa GridLayout && clean_gridlayout!(obj)
        obj isa Axis && delete!(obj)
        obj isa Label && delete!(obj)
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

function draw_bound!(
    ax,
    center::Tuple{<:Real,<:Real},
    smwh::Tuple{<:Integer, <:Integer},
    step;
    kwargs...
)
    x0, y0 = center
    rx, ry = (smwh .+ 0.5) .* step
    vtx_x = []
    vtx_y = []
    for (u, v) in [(+1, +1), (+1, -1), (-1, -1), (-1, +1), (+1, +1)]
        rxuv, ryuv = (u * rx, v * ry)
        push!(vtx_x, x0 + rxuv)
        push!(vtx_y, y0 + ryuv)
    end
    obj = lines!(ax, vtx_x, vtx_y; kwargs...)
    return obj
end

function matching_axes(dict_axes::Dict{String,Axis}, pattern)
    return [dict_axes[k] for k in keys(dict_axes) if occursin(pattern, k)]
end
