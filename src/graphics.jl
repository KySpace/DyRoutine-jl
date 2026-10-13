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
                    clamp(abs(t) / thres_alpha * (1 - alpha_base) + alpha_base, 0, 1)
            Oklch(1 - 0.8 * t, 0.24 * t, hue) |> c -> RGBAf(c, alpha)
        end
        for t in range(0, 1; length=256)]
    return clrmap
end

function clr_parabola(t, hue, light_maxchroma, chroma_max, light_min;
    thres_alpha=0.0, alpha_base=1.0, light_max=1.0, chroma_lightmax=0,
    hue_range=(0, 0), prescale=(t -> t))
    t = prescale(t)
    l = (1 - t) * (light_max - light_min) + light_min
    c = (chroma_lightmax - chroma_max) / (light_max - light_maxchroma)^2 * (l - light_maxchroma)^2 + chroma_max
    h = hue + (t - 0.5) * (hue_range[2] - hue_range[1])
    alpha = thres_alpha <= 0 || abs(t) > thres_alpha ?
        1.0 : clamp(abs(t) / thres_alpha * (1 - alpha_base) + alpha_base, 0, 1)
    return RGBAf(Oklch(clamp(l, 0, 1), clamp(c, 0, 1), clamp(h, 0, 360)), alpha)
end

function gen_clrmap_parabola(hue, light_maxchroma, chroma_max, light_min; kwargs...)
    return [clr_parabola(t, hue, light_maxchroma, chroma_max, light_min; kwargs...)
        for t in range(0, 1; length=256)]
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

function to_miscibility_clr(dens1, dens2, hue1, hue2;
    max=16,
    to_norm_each=false,
    thres_alpha=0.1,
    alpha_base=0.1,
)
    size(dens1) == size(dens2) ||
        throw(ArgumentError("dens1 and dens2 must have the same size"))
    norm = to_norm_each ? (d -> d ./ maximum(d)) : (d -> clamp.(d, 0, max) / max)
    alpha = t -> abs(t) > thres_alpha ? 1.0 :
        (abs(t) / thres_alpha * (1 - alpha_base) + alpha_base)
    dens_norm_1, dens_norm_2 = (dens1, dens2) |> ds -> map(norm, ds)
    shader = (a, b) -> Oklch(
        1 - (a + b) / 2,
        abs(a - b) * 0.24,
        a > b ? hue1 : hue2,
    ) |> c -> RGBAf(c, alpha(a + b))
    [shader(dens_norm_1[x, y], dens_norm_2[x, y])
        for x in axes(dens1, 1), y in axes(dens1, 2)]
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
