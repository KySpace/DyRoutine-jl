module CalibrationBenchmark

export calc_max_curve_disagreement

function find_stationary_positions(a::Real, b::Real, c::Real)::Vector{Float64}
    tolerance = 64eps(Float64) * max(abs(a), abs(b), abs(c))
    if abs(a) <= tolerance
        return abs(b) <= tolerance ? Float64[] : [-c / b]
    end
    discriminant = b^2 - 4a * c
    tolerance_discriminant = 64eps(Float64) * (b^2 + abs(4a * c))
    discriminant < -tolerance_discriminant && return Float64[]
    q = -0.5 * (b + copysign(sqrt(max(0.0, discriminant)), b))
    return iszero(q) ? [-b / (2a)] : [q / a, c / q]
end

"""
    calc_max_curve_disagreement(calc_frequency, calc_slope, ids_source, times_break)

Find the continuous maximum of `maximum(f_s(t)) - minimum(f_s(t))`.
All curves must be cubic or lower degree on each supplied time interval,
with continuous values and first derivatives. Check interval endpoints and
analytic stationary points of every pairwise cubic difference.
"""
function calc_max_curve_disagreement(
    calc_frequency::F, calc_slope::G,
    ids_source::AbstractVector{<:Integer}, times_break::AbstractVector{<:Real},
) where {F,G}
    length(ids_source) >= 2 && length(unique(ids_source)) == length(ids_source) ||
        throw(ArgumentError("Need at least two distinct source indices; got $ids_source"))
    length(times_break) >= 2 && all(isfinite, times_break) && all(>(0), diff(times_break)) ||
        throw(ArgumentError("Time breakpoints must be finite and strictly increasing"))
    result = (disagreement_mhz=-Inf, time_s=Float64(first(times_break)),
        id_upper=first(ids_source), id_lower=last(ids_source))
    for idx_left in 1:(length(times_break) - 1)
        time_left_s, time_right_s = times_break[idx_left:idx_left + 1]
        span_s = time_right_s - time_left_s
        for idx_source in 1:(length(ids_source) - 1), idx_other in (idx_source + 1):length(ids_source)
            id_source, id_other = ids_source[idx_source], ids_source[idx_other]
            difference_left = calc_frequency(id_source, time_left_s) - calc_frequency(id_other, time_left_s)
            difference_right = calc_frequency(id_source, time_right_s) - calc_frequency(id_other, time_right_s)
            slope_left = calc_slope(id_source, time_left_s) - calc_slope(id_other, time_left_s)
            slope_right = calc_slope(id_source, time_right_s) - calc_slope(id_other, time_right_s)
            all(isfinite, (difference_left, difference_right, slope_left, slope_right)) ||
                throw(ArgumentError("Nonfinite curve value or slope for sources $id_source and $id_other"))
            # Cubic Hermite coefficients on the normalized interval u ∈ [0, 1].
            c1 = span_s * slope_left
            c2 = 3(difference_right - difference_left) - span_s * (2slope_left + slope_right)
            c3 = 2(difference_left - difference_right) + span_s * (slope_left + slope_right)
            positions = [0.0, 1.0]
            append!(positions, filter(u -> 0 < u < 1, find_stationary_positions(3c3, 2c2, c1)))
            for u in positions
                time_s = time_left_s + u * span_s
                difference = calc_frequency(id_source, time_s) - calc_frequency(id_other, time_s)
                if abs(difference) > result.disagreement_mhz
                    result = (; disagreement_mhz=abs(difference), time_s=Float64(time_s),
                        id_upper=difference >= 0 ? id_source : id_other,
                        id_lower=difference >= 0 ? id_other : id_source)
                end
            end
        end
    end
    return result
end

end
