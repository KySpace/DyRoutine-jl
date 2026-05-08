# Analysis per-solo image
# When each shot contains one image, it is per-shot
# can also be applied to composite images from multiple solos, like stacked images
using LsqFit: curve_fit
using NaNStatistics: movmean
using FFTW
using DSP.Windows: hanning
using Printf
using Colors: Oklch

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

function fold_symmetric(arr::AbstractVector)::AbstractVector{<:Real}
    len = length(arr)
    tail_left, head_right = isodd(len) ? ((len + 1) / 2, (len + 1) / 2) : (len / 2, len / 2 + 1)
    return (arr[1:Int(tail_left)] .+ arr[end:-1:Int(head_right)]) |> a -> a ./ 2 |> reverse
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

function gen_win_hann_2d(smwh::Tuple{<:Real,<:Real})
    smw, smh = smwh
    win_x = hanning(2 * smw + 1)
    win_y = hanning(2 * smh + 1)
    win = win_y * win_x'
    size(win) == (2 * smh + 1, 2 * smw + 1) || throw(DimensionMismatch("Window size $(size(win)) does not match smwh $smwh"))
    return win
end

struct SoloEssentials
    dens2d::AbstractMatrix
    modl2d::AbstractMatrix
    prfl_modl::AbstractVector
    prfl_modl_norm_px::AbstractVector
    smwh::Tuple{<:Real,<:Real}
    smw_modl::Integer
    step_posi::Real
    step_modl::Real
end

struct SoloExtract
    essentials::SoloEssentials
    info_solo::Dict{String,Any}
    prfl_modl_norm_net_px::AbstractVector
end

function calc_solo_essn_2d(dens::AbstractMatrix, cent::Tuple{<:Real,<:Real}, smwh::Tuple{<:Real,<:Real}, smw_modl::Integer, px_in_um::Real)
    dens_roi = crop_center(dens, cent, smwh)
    x_cent = smwh[1] + 1
    step_posi = px_in_um
    step_modl = 1 / (2 * smwh[2] * px_in_um)
    modl_roi = dens_roi .* gen_win_hann_2d(smwh) |> fft |> fftshift |> c -> abs.(c)
    prfl_modl = modl_roi[:, x_cent-smw_modl:x_cent+smw_modl] |> m -> sum(m, dims=2) ./ (smw_modl * 2 + 1) |> vec
    prfl_modl_norm_px = prfl_modl ./ (sum(prfl_modl) * step_modl / 2)
    return SoloEssentials(dens_roi, modl_roi, prfl_modl, prfl_modl_norm_px, smwh, smw_modl, step_posi, step_modl)
end

function draw_solo_modl!(axs::Dict{String,Axis}, essn::SoloEssentials, info_solo)
    foreach(empty!, values(axs))
    modl2d_norm = essn.modl2d |> m -> m ./ (sum(m) * (essn.step_modl / 2)^2)
    x, y = essn.smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x, y) .* essn.step_posi
    x_modl, y_modl = (x, y) .* essn.step_modl
    y_modl_sm = (0:1:essn.smwh[2]) * essn.step_modl
    clrmap = gen_clrmap_solo(hue_theme_istp[info_solo["istp"]])
    heatmap!(axs["modl"], y_modl_sm, x_modl, modl2d_norm[essn.smwh[2]+1:end, :]; colorrange=(0, 10.0), colormap=clrmap, rasterize=true)
    lines!(axs["upright"], y_modl_sm, essn.prfl_modl_norm_px[essn.smwh[2]+1:end], color=:black, linewidth=1)
    lines!(axs["sideway"], essn.prfl_modl_norm_px[essn.smwh[2]+1:end], y_modl_sm, color=:black, linewidth=1)
    axs["sideway"] |> hidedecorations!
    axs["modl"] |> hidedecorations!
    axs["upright"].yticklabelsvisible = false
    axs["upright"].xticklabelsvisible = false
    xlims!(axs["upright"], 0, 0.8)
    xlims!(axs["modl"], 0, 0.8)
    ylims!(axs["upright"], 0, 2.0)
    ylims!(axs["modl"], (-10.5, 10.5) .* essn.step_modl)
    ylims!(axs["sideway"], 0.2, 0.5)
    xlims!(axs["sideway"], 0.0, 2.0)
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    vlines!(axs["upright"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    hlines!(axs["sideway"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    text!(axs["modl"], 0.55, 0.16; text="$(info_solo["t_hold"]) ms | rep $(info_solo["repeat"])", color=:black, fontsize=16, align=(:center, :bottom))
end

function draw_solo_essn_2d!(axs::Dict{String,Axis}, essn::SoloEssentials, info_solo)
    foreach(empty!, values(axs))
    modl2d_norm = essn.modl2d |> m -> m ./ (sum(m) * (essn.step_modl / 2)^2)
    x, y = essn.smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x, y) .* essn.step_posi
    x_modl, y_modl = (x, y) .* essn.step_modl
    y_modl_sm = (0:1:essn.smwh[2]) * essn.step_modl
    clrmap = gen_clrmap_solo(hue_theme_istp[info_solo["istp"]])
    heatmap!(axs["dens"], x_posi, y_posi, essn.dens2d'; colorrange=(0, 16.0), colormap=clrmap)
    heatmap!(axs["modl"], y_modl_sm, x_modl, modl2d_norm[essn.smwh[2]+1:end, :]; colorrange=(0, 10.0), colormap=:binary)
    axs["dens"].aspect = DataAspect()
    # axs["modl"].aspect = DataAspect()
    ylims!(axs["prfl_ft"], 0, 2.5)
    xlims!(axs["prfl_ft"], 0, 0.8)
    xlims!(axs["modl"], 0, 0.8)
    ylims!(axs["modl"], -0.5, 0.5)
    axs["prfl_ft"].yticksvisible = false
    axs["prfl_ft"].yticklabelsvisible = false
    axs["modl"] |> hidedecorations!
    # axs["dens"] |> hidedecorations!
    lines!(axs["prfl_ft"], y_modl_sm, essn.prfl_modl_norm_px |> fold_symmetric; color=:black)
    vlines!(axs["prfl_ft"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    hlines!(axs["modl"], [-10.5, 10.5] .* essn.step_modl; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    text!(axs["dens"], 0, 14; text=@sprintf("%i ms | rep %i", info_solo["t_hold"], info_solo["repeat"]), color=:black, fontsize=24, align=(:center, :bottom))
end

function fit_prfl_modl_twinpeak_decay_1d(coor, prfl)
    # parameters: [MainPeak.Height MainPeak.Width SidePeak.Height SidePeak.Width SidePeak.Pos Decay.Height Decay.Length]
    model(k, p) = p[1] .* exp.(-k .^ 2 ./ (2 .* p[2] .^ 2)) .+ p[3] .* exp.(-(k .- p[5]) .^ 2 ./ (2 .* p[4] .^ 2)) .+ p[6] .* exp.(-abs.(k) ./ p[7])
    p_init = [3.0, 0.1, 0.5, 0.05, 0.3, 0.5, 0.8]
    p_upper = [Inf, 0.3, 2.0, 0.100, 0.37, Inf, 5.0]
    p_lower = [2.0, 0.0, 0.0, 0.018, 0.23, 0.0, 0.5]
    fit = curve_fit(model, coor, prfl, p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    return Dict(
        "fit" => fit,
        "model" => model,
        "params" => params_fit,
        "tail" => k -> params_fit[6] .* exp.(-abs.(k) ./ params_fit[7])
    )
end

function fit_prfl_modl_twinpeak_1d(coor, prfl)
    # parameters: [MainPeak.Height MainPeak.Width SidePeak.Height SidePeak.Width SidePeak.Pos Decay.Height Decay.Length]
    model(k, p) = p[1] .* exp.(-k .^ 2 ./ (2 .* p[2] .^ 2)) .+ p[3] .* exp.(-(k .- p[5]) .^ 2 ./ (2 .* p[4] .^ 2))
    p_init = [3.0, 0.1, 0.5, 0.05, 0.3]
    p_upper = [Inf, 0.3, 2.0, 0.100, 0.37]
    p_lower = [2.0, 0.0, 0.0, 0.018, 0.23]
    fit = curve_fit(model, coor, prfl, p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    return Dict(
        "fit" => fit,
        "model" => model,
        "params" => params_fit,
    )
end
