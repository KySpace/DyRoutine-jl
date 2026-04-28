using LsqFit: curve_fit
using NaNStatistics: movmean

function subtract_corner_mean(arr::AbstractMatrix, wh_corner::Tuple{<:Integer,<:Integer})
    (h_corner, w_corner) = wh_corner
    h_corner > 0 || throw(ArgumentError("corner_height must be positive."))
    w_corner > 0 || throw(ArgumentError("corner_width must be positive."))

    h_im, w_im = size(arr)
    2 * h_corner <= h_im || throw(ArgumentError("corner_height is too large for array height $h_im."))
    2 * w_corner <= w_im || throw(ArgumentError("corner_width is too large for array width $w_im."))

    tl = @view arr[1:h_corner, 1:w_corner]
    tr = @view arr[1:h_corner, w_im-w_corner+1:w_im]
    bl = @view arr[h_im-h_corner+1:h_im, 1:w_corner]
    br = @view arr[h_im-h_corner+1:h_im, w_im-w_corner+1:w_im]

    corner_mean = (sum(tl) + sum(tr) + sum(bl) + sum(br)) / (4 * h_corner * w_corner)
    return arr .- corner_mean
end

function crop_center(
    arr::AbstractMatrix,
    xy::Tuple{<:Integer,<:Integer},
    smwh::Tuple{<:Integer,<:Integer},
)::AbstractMatrix{<:Real}
    x, y = xy
    smw, smh = smwh

    smw >= 0 || throw(ArgumentError("smw must be nonnegative."))
    smh >= 0 || throw(ArgumentError("smh must be nonnegative."))

    height, width = size(arr)
    1 <= x <= width || throw(ArgumentError("x=$x is out of bounds for array width $width."))
    1 <= y <= height || throw(ArgumentError("y=$y is out of bounds for array height $height."))

    left = x - smw
    right = x + smw
    top = y - smh
    bottom = y + smh

    left >= 1 || throw(ArgumentError("Crop extends past the left edge: x=$x, smw=$smw."))
    right <= width || throw(ArgumentError("Crop extends past the right edge: x=$x, smw=$smw, width=$width."))
    top >= 1 || throw(ArgumentError("Crop extends past the top edge: y=$y, smh=$smh."))
    bottom <= height || throw(ArgumentError("Crop extends past the bottom edge: y=$y, smh=$smh, height=$height."))

    return @view arr[top:bottom, left:right]
end

function calc_dens_sum(dens::AbstractMatrix{<:Real})
    return sum(dens; dims=(1, 2))
end

function moving_average_with_positions(prfl::AbstractVector, len_avg::Integer)::Tuple{AbstractVector{<:Real},AbstractVector{<:Real}}
    len_avg > 0 || throw(ArgumentError("len_avg must be positive."))
    n = length(prfl)
    len_avg <= n || throw(ArgumentError("len_avg=$len_avg exceeds profile length $n."))

    prfl_avg = movmean(Float64.(prfl), len_avg)
    pos_avg = collect((len_avg+1)/2:(n-(len_avg-1)/2))

    return prfl_avg, pos_avg
end

function find_peak_position_moving(prfl::AbstractVector; len_avg::Integer=10)::Integer
    prfl_avg, pos_avg = moving_average_with_positions(prfl, len_avg)
    return round(Int, pos_avg[argmax(prfl_avg)])
end

gaussian_1d(x, p) = @. p[1] * exp(-((x - p[2])^2) / (2 * p[3]^2))

function gaussian_fit_center_1d(prfl::AbstractVector)
    n = length(prfl)
    n > 0 || throw(ArgumentError("Profile must be nonempty."))

    x = collect(1.0:n)
    y = Float64.(prfl)
    amp0 = maximum(y)
    amp0 > 0 || throw(ArgumentError("Profile must contain a positive peak for Gaussian fitting."))

    center0 = Float64(argmax(y))
    sigma0 = clamp(n / 4, 2.0, float(n))
    p0 = [amp0, center0, sigma0]
    lower = [amp0 / 100, 0.0, min(2.0, float(n))]
    upper = [amp0, float(n), float(n)]

    fit = curve_fit(gaussian_1d, x, y, p0; lower=lower, upper=upper)
    return fit.param[2]
end

function find_positive_cluster_center(
    arr::AbstractMatrix;
    len_avg::Integer=10,
    smwh::Tuple{Integer,Integer}=(),
)::Tuple{<:Real,<:Real}
    smw, smh = smwh
    cx_coarse = find_peak_position_moving(vec(sum(arr; dims=1)); len_avg=len_avg)
    cy_coarse = find_peak_position_moving(vec(sum(arr; dims=2)); len_avg=len_avg)

    cropped = crop_center(arr, (cx_coarse, cy_coarse), (smw, smh))
    left = cx_coarse - smw
    top = cy_coarse - smh

    cx_local = gaussian_fit_center_1d(vec(sum(cropped; dims=1)))
    cy_local = gaussian_fit_center_1d(vec(sum(cropped; dims=2)))

    return left - 1 + cx_local, top - 1 + cy_local
end
