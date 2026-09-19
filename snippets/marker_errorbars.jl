using CairoMakie

"""
    marker_errorbars!(ax, x, y, yerror; kwargs...)
    marker_errorbars!(ax, x, y, yerror_low, yerror_high; kwargs...)

Draw scatter markers with two capless vertical error-bar segments per marker.
Each segment starts at the actual lower or upper edge of the marker and ends at
`y - yerror_low` or `y + yerror_high`. Consequently, an error shorter than the
corresponding marker half-height points inward across the marker, while a longer
error points outward.

The marker-edge positions react to changes in the axis limits, y-axis scale,
and pixel size. Linear, logarithmic, and other Makie scales with a defined
`inverse_transform` are supported. Markers must be unrotated and use pixel
space.

Keyword arguments `marker`, `markersize`, `color`, `strokecolor`, and
`strokewidth` are passed to `scatter!`. `errorlinewidth` controls the width of
the two error segments. `markeredgeoffset` adds an optional pixel adjustment
outside the visible marker stroke. Other keyword arguments are forwarded to
`scatter!`. The marker is drawn first and the segments second, so inward
segments remain visible on top of the marker.

Points follow Makie's input draw order: each later marker covers earlier points'
error segments, while its own error segments remain above it.

Returns a named tuple `(scatter, marker_overlays, errorbars, segments)`. The
first field is the combined scatter plot used for marker geometry and legends;
the two middle fields contain the plots used for ordered compositing.
"""
function marker_errorbars!(
    ax::Axis,
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    yerror::Union{Real,AbstractVector{<:Real}};
    kwargs...
)
    return marker_errorbars!(ax, x, y, yerror, yerror; kwargs...)
end

function marker_errorbars!(
    ax::Axis,
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    yerror_low::Union{Real,AbstractVector{<:Real}},
    yerror_high::Union{Real,AbstractVector{<:Real}};
    marker=:circle,
    markersize=12,
    color=:white,
    strokecolor=:black,
    strokewidth=1.5,
    errorlinewidth=strokewidth,
    markeredgeoffset=0,
    kwargs...
)
    n = length(x)
    length(y) == n || throw(DimensionMismatch(
        "x has length $n, but y has length $(length(y))"
    ))
    _check_marker_error_length(yerror_low, n, "yerror_low")
    _check_marker_error_length(yerror_high, n, "yerror_high")
    _check_marker_attribute_length(marker, n, "marker")
    _check_marker_attribute_length(markersize, n, "markersize")
    _check_marker_attribute_length(color, n, "color")
    _check_marker_attribute_length(strokecolor, n, "strokecolor")
    _check_marker_attribute_length(strokewidth, n, "strokewidth")
    _check_marker_attribute_length(errorlinewidth, n, "errorlinewidth")
    _check_marker_attribute_length(markeredgeoffset, n, "markeredgeoffset")
    _check_nonnegative_errors(yerror_low, "yerror_low")
    _check_nonnegative_errors(yerror_high, "yerror_high")

    _check_scale_domain(ax.yscale[], y, yerror_low, yerror_high)
    haskey(kwargs, :markerspace) && kwargs[:markerspace] != :pixel &&
        throw(ArgumentError("markerspace must be :pixel"))
    haskey(kwargs, :rotation) && !_is_unrotated(kwargs[:rotation]) &&
        throw(ArgumentError("rotated markers are not currently supported"))

    scatter_plot = scatter!(
        ax,
        x,
        y;
        marker,
        markersize,
        color,
        strokecolor,
        strokewidth,
        kwargs...,
    )

    segments = lift(
        ax.finallimits,
        ax.scene.viewport,
        ax.yscale,
        scatter_plot.marker,
        scatter_plot.markersize,
        scatter_plot.strokewidth,
    ) do limits, viewport, yscale, marker_value, markersize_value, strokewidth_value
        axis_height = Float64(widths(viewport)[2])
        axis_height > 0 || return Point2f[]
        y_limit_low = Float64(limits.origin[2])
        y_limit_high = y_limit_low + Float64(widths(limits)[2])
        y_scaled_low = Float64(yscale(y_limit_low))
        y_scaled_high = Float64(yscale(y_limit_high))
        scaled_per_pixel = (y_scaled_high - y_scaled_low) / axis_height
        inverse_yscale = Makie.inverse_transform(yscale)

        points = Vector{Point2f}(undef, 4n)
        for i in eachindex(x, y)
            marker_i = _marker_value(marker_value, i)
            markersize_i = _marker_size_value(markersize_value, i)
            strokewidth_i = Float64(_marker_value(strokewidth_value, i))
            edgeoffset_i = Float64(_marker_value(markeredgeoffset, i))
            y_pixel_1, y_pixel_2 = _visible_marker_y_edges(
                marker_i,
                markersize_i,
                strokewidth_i + edgeoffset_i,
            )
            y_scaled = Float64(yscale(y[i]))
            y_bottom = inverse_yscale(
                y_scaled + min(y_pixel_1, y_pixel_2) * scaled_per_pixel
            )
            y_top = inverse_yscale(
                y_scaled + max(y_pixel_1, y_pixel_2) * scaled_per_pixel
            )
            error_low_i = Float64(_marker_value(yerror_low, i))
            error_high_i = Float64(_marker_value(yerror_high, i))
            j = 4i - 3
            points[j] = Point2f(x[i], y_bottom)
            points[j + 1] = Point2f(x[i], y[i] - error_low_i)
            points[j + 2] = Point2f(x[i], y_top)
            points[j + 3] = Point2f(x[i], y[i] + error_high_i)
        end
        return points
    end

    error_plots = []
    marker_overlays = []
    overlay_kwargs = merge((; kwargs...), (; label=nothing, inspectable=false))
    segment_autolimits = (
        xautolimits=get(kwargs, :xautolimits, true),
        yautolimits=get(kwargs, :yautolimits, true),
    )
    for i in eachindex(x, y)
        segment_i = lift(points -> points[(4i - 3):(4i)], segments)
        push!(error_plots, linesegments!(
            ax,
            segment_i;
            color=_marker_value(strokecolor, i),
            linewidth=_marker_value(errorlinewidth, i),
            label=nothing,
            inspectable=false,
            segment_autolimits...,
        ))

        i == n && continue
        j = i + 1
        push!(marker_overlays, scatter!(
            ax,
            [x[j]],
            [y[j]];
            marker=_marker_value(marker, j),
            markersize=_marker_size_value(markersize, j),
            color=_marker_value(color, j),
            strokecolor=_marker_value(strokecolor, j),
            strokewidth=_marker_value(strokewidth, j),
            overlay_kwargs...,
        ))
    end
    return (
        scatter=scatter_plot,
        marker_overlays,
        errorbars=error_plots,
        segments,
    )
end

function _check_scale_domain(yscale, y, yerror_low, yerror_high)
    interval = Makie.defined_interval(yscale)
    for i in eachindex(y)
        error_low_i = _marker_value(yerror_low, i)
        error_high_i = _marker_value(yerror_high, i)
        values = (y[i], y[i] - error_low_i, y[i] + error_high_i)
        all(value -> value in interval, values) || throw(ArgumentError(
            "y and both error endpoints must lie in the y-scale domain $interval; " *
            "point $i gives $(values)"
        ))
    end
    Makie.inverse_transform(yscale)
    return nothing
end

_marker_value(value::AbstractVector, i::Integer) = value[i]
_marker_value(value, ::Integer) = value

_marker_size_value(value::Makie.VecTypes, ::Integer) = value
_marker_size_value(value::Tuple{<:Real,<:Real}, ::Integer) = value
_marker_size_value(value::AbstractVector, i::Integer) = value[i]
_marker_size_value(value, ::Integer) = value

function _visible_marker_y_edges(marker::Makie.BezierPath, markersize, padding::Real)
    bbox = Makie.bbox(marker)
    markersize_y = markersize isa Real ? markersize : markersize[2]
    y_bottom = markersize_y * bbox.origin[2] - padding
    y_top = markersize_y * (bbox.origin[2] + widths(bbox)[2]) + padding
    return Float64(y_bottom), Float64(y_top)
end

function _visible_marker_y_edges(marker, markersize, padding::Real)
    throw(ArgumentError(
        "visible marker-edge compensation requires a BezierPath marker; " *
        "got $(typeof(marker))"
    ))
end

function _check_marker_error_length(value, n::Integer, name::String)
    value isa AbstractVector && length(value) != n && throw(DimensionMismatch(
        "$name has length $(length(value)), expected $n"
    ))
    return nothing
end

function _check_marker_attribute_length(value, n::Integer, name::String)
    value isa AbstractVector &&
        !(value isa AbstractString) &&
        !(value isa Makie.VecTypes) &&
        length(value) != n &&
        throw(DimensionMismatch("$name has length $(length(value)), expected $n"))
    return nothing
end

function _check_nonnegative_errors(value, name::String)
    values = value isa AbstractVector ? value : (value,)
    for error in values
        isfinite(error) || throw(ArgumentError("$name must be finite, got $error"))
        error >= 0 || throw(ArgumentError("$name must be nonnegative, got $error"))
    end
    return nothing
end

_is_unrotated(rotation::Real) = iszero(rotation)
_is_unrotated(rotation::Makie.Billboard) = true
_is_unrotated(rotation) = false
