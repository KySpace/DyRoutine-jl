import CairoMakie
using GLMakie
using Colors: Oklch, RGB
using LaTeXStrings

Γ = 0.04
t_lim = (-2, 12)
f_lim = (-2, 3.5)
f_fade_start = 2
f_fade_stop = 3.5
period_chop = 3
δt = 0.08
alpha_min = 0.65
alpha_max = 0.96
alpha_gamma = 0.1

height_l_base = 0.98
height_l_peak = 0.82
height_chroma_base = 0.04
height_chroma_peak = 0.16
height_hue = 277

orientation_l = 0.90
orientation_chroma = 0.0
orientation_hue = 190
orientation_strength = 0 # -0.1
orientation_gamma = 0.85

profile_plane_t = 14
profile_plane_f_lim = (-2, 2)
profile_plane_z_lim = (-0.3, 1.2)
profile_line_l = 0.42
profile_line_chroma = 0.10
profile_line_hue = 277

cross_plane_f = 4.5
crosssection_f_pos = 1
crosssection_f_neg = -1
cross_line_l = 0.48
cross_line_chroma = 0.12
cross_line_hue_pos = 30
cross_line_hue_neg = 215
cross_projection_repetitions = 3.5
cross_hold_pos_pulse_idx = 4
cross_ellipsis_t_offset = 0.35
cross_ellipsis_t_gap = 0.35
cross_ellipsis_z = 0.6
cross_ellipsis_markersize = 7
cross_label_rep_idx = 2
cross_label_z_offset = 0.2
cross_label_fontsize = 18
cross_marker_z_offset = 0.1
cross_marker_tick_length = 0.08
cross_marker_arrowhead_length = 0.18
cross_marker_arrowhead_height = 0.08
track_t_lim = (0, 12)
track_z = 0.01
track_linewidth = 1.0
track_line_l = 0.48
track_line_chroma = 0.0
track_line_hue_pos = 30
track_line_hue_neg = 215
track_line_l_zero = 0.48
track_line_alpha = 0.95
track_line_alpha_zero = 0.50

plane_l = 0.94
plane_chroma = 0.015
plane_hue = 245
box_l = 0.82
box_chroma = 0.03
box_hue = 245
box_alpha = 0.40

azimuth_view = -2.3
elevation_view = 0.41
aspect_xy = (1.0, 0.45)
aspect_z_over_xy = 0.18
path_plot = joinpath(@__DIR__, "freq_chop_surface.pdf")
to_save_pdf = true
surface_rasterize_scale = 3
to_show_gl = true
to_wait_for_window = true

GLMakie.activate!()

lorentz(f::Real, f0::Real) = 1 / (1 + ((f - f0) / Γ)^2)

function ramp_gate(phase::Real, start::Real, stop::Real, δt::Real)
    start < phase < stop || return 0.0
    δt <= 0 && return 1.0

    width_rise = min(δt, stop - start)
    return clamp((phase - start) / width_rise, 0, 1)
end

function freq_chop(f::Real, t::Real)
    t < 0 && return 0.0

    phase = mod(t / period_chop, 1)
    pos = ramp_gate(phase, 0.02, 0.48, δt)
    neg = ramp_gate(phase, 0.52, 0.98, δt)

    return pos * lorentz(f, 1) + neg * lorentz(f, -1)
end

function freq_chop_cross_pos(f::Real, t::Real)
    t < 0 && return 0.0

    idx_pulse = floor(Int, t / period_chop) + 1
    phase = mod(t / period_chop, 1)
    pos = ramp_gate(phase, 0.02, 0.48, δt)
    if idx_pulse == cross_hold_pos_pulse_idx && phase >= 0.02 + δt
        pos = 1.0
    end

    return pos * lorentz(f, 1)
end

function calc_cross_projection(t::Real)
    z_pos = freq_chop_cross_pos(crosssection_f_pos, t)
    z_neg = freq_chop(crosssection_f_neg, t)
    z = max(z_pos, z_neg)
    color_id = z <= eps(Float64) ? :zero : (z_pos >= z_neg ? :pos : :neg)

    return z, color_id
end

function calc_peak_track(t::Real)
    phase = mod(t / period_chop, 1)

    if 0.02 < phase < 0.48
        return 1.0, :pos
    elseif 0.52 < phase < 0.98
        return -1.0, :neg
    elseif 0.48 <= phase <= 0.52
        return 1 - 2 * (phase - 0.48) / (0.52 - 0.48), :zero
    else
        phase_gap = phase > 0.98 ? phase - 0.98 : phase + 0.02
        return -1 + 2 * phase_gap / (0.04), :zero
    end
end

function calc_plateau_t_span(idx_rep::Integer, phase_start::Real, phase_stop::Real)
    t_start = period_chop * (idx_rep - 1 + phase_start + δt)
    t_stop = period_chop * (idx_rep - 1 + phase_stop)

    return t_start, t_stop
end

function draw_width_marker!(
    ax,
    t_span::Tuple{<:Real, <:Real},
    f::Real,
    z_plateau::Real;
    color,
)
    t_start, t_stop = t_span
    z_arrow = z_plateau + cross_marker_z_offset
    z_tick1 = z_arrow - cross_marker_tick_length / 2
    z_tick2 = z_arrow + cross_marker_tick_length / 2
    t_head = min(cross_marker_arrowhead_length, (t_stop - t_start) / 4)
    z_head = cross_marker_arrowhead_height / 2

    lines!(ax, [t_start, t_start], [f, f], [z_tick1, z_tick2]; color, linewidth=1)
    lines!(ax, [t_stop, t_stop], [f, f], [z_tick1, z_tick2]; color, linewidth=1)
    lines!(ax, [t_start, t_stop], [f, f], [z_arrow, z_arrow]; color, linewidth=1)
    mesh!(
        ax,
        [
            Point3f(t_start, f, z_arrow),
            Point3f(t_start + t_head, f, z_arrow - z_head),
            Point3f(t_start + t_head, f, z_arrow + z_head),
        ],
        [1 2 3];
        color,
    )
    mesh!(
        ax,
        [
            Point3f(t_stop, f, z_arrow),
            Point3f(t_stop - t_head, f, z_arrow - z_head),
            Point3f(t_stop - t_head, f, z_arrow + z_head),
        ],
        [1 2 3];
        color,
    )
end

freq_chop(f::AbstractArray, t::AbstractArray) = freq_chop.(f, t)
fade_alpha(t::Real) = t <= 8 ? 1.0 : t < 11 ? (11 - t) / 3 : 0.0
fade_f_alpha(f::Real) = f <= f_fade_start ? 1.0 : f < f_fade_stop ? (f_fade_stop - f) / (f_fade_stop - f_fade_start) : 0.0

function calc_gradient(arr::AbstractMatrix{<:Real}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    grad_x = similar(arr, Float64)
    grad_y = similar(arr, Float64)

    for j in axes(arr, 2)
        grad_x[1, j] = (arr[2, j] - arr[1, j]) / (xs[2] - xs[1])
        grad_x[end, j] = (arr[end, j] - arr[end - 1, j]) / (xs[end] - xs[end - 1])
        for i in 2:size(arr, 1)-1
            grad_x[i, j] = (arr[i + 1, j] - arr[i - 1, j]) / (xs[i + 1] - xs[i - 1])
        end
    end

    for i in axes(arr, 1)
        grad_y[i, 1] = (arr[i, 2] - arr[i, 1]) / (ys[2] - ys[1])
        grad_y[i, end] = (arr[i, end] - arr[i, end - 1]) / (ys[end] - ys[end - 1])
        for j in 2:size(arr, 2)-1
            grad_y[i, j] = (arr[i, j + 1] - arr[i, j - 1]) / (ys[j + 1] - ys[j - 1])
        end
    end

    return grad_x, grad_y
end

function add_rgb(c1, c2, weight::Real)
    rgb1 = RGB(c1)
    rgb2 = RGB(c2)

    return RGB(
        clamp(rgb1.r + weight * rgb2.r, 0, 1),
        clamp(rgb1.g + weight * rgb2.g, 0, 1),
        clamp(rgb1.b + weight * rgb2.b, 0, 1),
    )
end

function draw_t_plane!(
    ax,
    t::Real,
    f_span::Tuple{<:Real, <:Real},
    z_span::Tuple{<:Real, <:Real};
    color_box,
)
    f1, f2 = f_span
    z1, z2 = z_span

    lines!(ax, fill(t, 5), [f1, f2, f2, f1, f1], [z1, z1, z2, z2, z1]; color=color_box, linewidth=1.5)
end

function draw_f_plane!(
    ax,
    f::Real,
    t_span::Tuple{<:Real, <:Real},
    z_span::Tuple{<:Real, <:Real};
    color_box,
)
    return nothing
end

t_vec = range(t_lim...; length=960)
f_vec = range(f_lim...; length=280)
amp_tf = [freq_chop(f, t) for t in t_vec, f in f_vec]
grad_t, grad_f = calc_gradient(amp_tf, t_vec, f_vec)
angle_tf = @. atan(hypot(grad_t, grad_f)) / (π / 2)
clr_tf = [
    RGBAf(
        add_rgb(
            Oklch(
                (1 - amp_tf[i, j]) * height_l_base + amp_tf[i, j] * height_l_peak,
                (1 - amp_tf[i, j]) * height_chroma_base + amp_tf[i, j] * height_chroma_peak,
                height_hue,
            ),
            Oklch(orientation_l, orientation_chroma, orientation_hue),
            orientation_strength * angle_tf[i, j]^orientation_gamma,
        ),
        fade_alpha(t_vec[i]) * fade_f_alpha(f_vec[j]) * (alpha_min + (alpha_max - alpha_min) * amp_tf[i, j]^alpha_gamma),
    )
    for i in eachindex(t_vec), j in eachindex(f_vec)
]

function make_freq_chop_figure(; backgroundcolor=:white)
    fig = Figure(; size=(980, 720), backgroundcolor)
    t_axis_lim = (t_lim[1], max(t_lim[2], profile_plane_t))
    f_axis_lim = (f_lim[1], max(f_lim[2], cross_plane_f))
    z_axis_lim = profile_plane_z_lim
    ax = Axis3(
        fig[1, 1];
        xlabel="",
        ylabel="",
        zlabel="",
        limits=(t_axis_lim, f_axis_lim, z_axis_lim),
        aspect=(aspect_xy[1], aspect_xy[2], aspect_z_over_xy),
        azimuth=azimuth_view,
        elevation=elevation_view,
        perspectiveness=0.0,
        xspinesvisible=false,
        yspinesvisible=false,
        zspinesvisible=false,
        xgridvisible=false,
        ygridvisible=false,
        zgridvisible=false,
        xticksvisible=false,
        yticksvisible=false,
        zticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        zticklabelsvisible=false,
    )

    surface!(
        ax,
        t_vec,
        f_vec,
        amp_tf;
        color=clr_tf,
        transparency=true,
        shading=NoShading,
        rasterize=surface_rasterize_scale,
    )

    clr_box = RGBAf(Oklch(box_l, box_chroma, box_hue), box_alpha)
    clr_profile = RGBAf(Oklch(profile_line_l, profile_line_chroma, profile_line_hue), 0.85)
    clr_cross_pos = RGBAf(Oklch(cross_line_l, cross_line_chroma, cross_line_hue_pos), 0.95)
    clr_cross_neg = RGBAf(Oklch(cross_line_l, cross_line_chroma, cross_line_hue_neg), 0.95)
    clr_cross_zero = RGBAf(Oklch(cross_line_l, 0, cross_line_hue_pos), 0.50)
    clr_track_pos = RGBAf(Oklch(track_line_l, track_line_chroma, track_line_hue_pos), track_line_alpha)
    clr_track_neg = RGBAf(Oklch(track_line_l, track_line_chroma, track_line_hue_neg), track_line_alpha)
    clr_track_zero = RGBAf(Oklch(track_line_l_zero, 0, track_line_hue_pos), track_line_alpha_zero)

    t_track_vec = range(track_t_lim...; length=720)
    track_vals = [calc_peak_track(t) for t in t_track_vec]
    f_track = first.(track_vals)
    color_ids_track = last.(track_vals)
    for i in 1:length(t_track_vec)-1
        clr_track = if :pos in color_ids_track[i:i+1]
            clr_track_pos
        elseif :neg in color_ids_track[i:i+1]
            clr_track_neg
        else
            clr_track_zero
        end
        lines!(
            ax,
            t_track_vec[i:i+1],
            f_track[i:i+1],
            fill(track_z, 2);
            color=clr_track,
            linewidth=track_linewidth,
        )
    end

    draw_t_plane!(
        ax,
        profile_plane_t,
        profile_plane_f_lim,
        profile_plane_z_lim;
        color_box=clr_box,
    )
    f_profile = range(profile_plane_f_lim...; length=360)
    z_profile = [lorentz(f, crosssection_f_pos) + lorentz(f, crosssection_f_neg) for f in f_profile]
    lines!(ax, fill(profile_plane_t, length(f_profile)), f_profile, z_profile; color=clr_profile, linewidth=1)

    draw_f_plane!(
        ax,
        cross_plane_f,
        t_lim,
        profile_plane_z_lim;
        color_box=clr_box,
    )
    t_cross_vec = range(0, period_chop * cross_projection_repetitions; length=560)
    cross_vals = [calc_cross_projection(t) for t in t_cross_vec]
    z_cross = first.(cross_vals)
    color_ids_cross = last.(cross_vals)
    for i in 1:length(t_cross_vec)-1
        clr_cross = if :pos in color_ids_cross[i:i+1]
            clr_cross_pos
        elseif :neg in color_ids_cross[i:i+1]
            clr_cross_neg
        else
            clr_cross_zero
        end
        lines!(
            ax,
            t_cross_vec[i:i+1],
            fill(cross_plane_f, 2),
            z_cross[i:i+1];
            color=clr_cross,
            linewidth=1,
        )
    end

    t_ellipsis = period_chop * cross_projection_repetitions .+
        cross_ellipsis_t_offset .+
        cross_ellipsis_t_gap .* (0:2)
    scatter!(
        ax,
        t_ellipsis,
        fill(cross_plane_f, length(t_ellipsis)),
        fill(cross_ellipsis_z, length(t_ellipsis));
        color=clr_box,
        markersize=cross_ellipsis_markersize,
    )
    z_label = 1 + cross_label_z_offset
    t_span_pos = calc_plateau_t_span(cross_label_rep_idx, 0.02, 0.48)
    t_span_neg = calc_plateau_t_span(cross_label_rep_idx, 0.52, 0.98)
    t_label_pos = sum(t_span_pos) / 2
    t_label_neg = sum(t_span_neg) / 2
    draw_width_marker!(ax, t_span_pos, cross_plane_f, 1; color=clr_cross_pos)
    draw_width_marker!(ax, t_span_neg, cross_plane_f, 1; color=clr_cross_neg)
    text!(
        ax,
        [t_label_pos],
        [cross_plane_f],
        [z_label];
        text=[L"t_{162}"],
        color=clr_cross_pos,
        align=(:center, :bottom),
        fontsize=cross_label_fontsize,
    )
    text!(
        ax,
        [t_label_neg],
        [cross_plane_f],
        [z_label];
        text=[L"t_{164}"],
        color=clr_cross_neg,
        align=(:center, :bottom),
        fontsize=cross_label_fontsize,
    )

    return fig, ax
end

if to_save_pdf
    fig_svg, _ = make_freq_chop_figure(; backgroundcolor=:transparent)
    save(path_plot, fig_svg; backend=CairoMakie)
end

if to_show_gl
    fig_gl, ax_gl = make_freq_chop_figure()
    screen = display(fig_gl)
    to_wait_for_window && wait(screen)
    println("azimuth = ", ax_gl.azimuth[])
    println("elevation = ", ax_gl.elevation[])
end
