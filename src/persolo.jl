# Analysis per-solo image
# When each shot contains one image, it is per-shot
# can also be applied to composite images from multiple solos, like stacked images
using LsqFit
using NaNStatistics: movmean
using FFTW
using DSP.Windows: hanning
using Printf
using Colors: Oklch
using NumericalIntegration: integrate

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
    prfl_modl_norm_tailess_px::AbstractVector
    sidepeak::Dict{String,Real}
    fit_tailess::Dict
    moments_modl::Dict{String,Real}
    fit_dens_2d::Dict
    envelope::Dict{String}
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

function calc_solo_extr(essn::SoloEssentials, fit_stack::Dict)
    x, y = essn.smwh |> s -> map(u -> (-u:1:u), s)
    x_modl, y_modl = (x, y) .* essn.step_modl
    x_posi, y_posi = (x, y) .* essn.step_posi
    sel_moment = y -> (y .> 0.2) .& (y .< 0.4)
    mask_mmt = sel_moment(y_modl)
    prfl_tailess = essn.prfl_modl_norm_px - fit_stack["tail"](y_modl)
    fit_tailess = fit_prfl_modl_twinpeak_1d(y_modl, prfl_tailess, (y_modl .> 0.1) .& (y_modl .< 0.5))
    sidepeak = Dict(
        "height" => fit_tailess["params"][3],
        "width" => fit_tailess["params"][4],
        "wavenum" => fit_tailess["params"][5],
        "weight" => sqrt(2 * pi) * fit_tailess["params"][3] * fit_tailess["params"][4]
    )
    moments = calc_prfl_moment(y_modl[mask_mmt], prfl_tailess[mask_mmt])
    fit_dens = fit_dens2d_gaussian_elliptic_disk(x_posi, y_posi, essn.dens2d, :)
    envelope = Dict(
        "max" => fit_dens["params"][1],
        "cent" => (fit_dens["params"][2], fit_dens["params"][3]),
        "size" => (fit_dens["params"][4], fit_dens["params"][5]),
        "rotation" => fit_dens["params"][6]
    )
    return SoloExtract(essn, prfl_tailess, sidepeak, fit_tailess, moments, fit_dens, envelope)
end

function calc_prfl_moment(coor, prfl)
    @assert length(coor) == length(prfl) "coordinate and profile length mismatch"
    prfl = prfl .- minimum(prfl)
    ntgr_over_coor = y -> integrate(coor, y)
    weight = prfl |> ntgr_over_coor
    height = weight / (coor[end] - coor[1])
    expval = coor .* prfl |> ntgr_over_coor |> u -> u ./ weight
    var = (coor .- expval) .^ 2 .* prfl |> ntgr_over_coor |> sqrt |> u -> u ./ weight
    return Dict(
        "weight" => weight,
        "wavenum" => expval,
        "width" => var,
        "height" => height,
    )
end

function draw_solo_modl!(axs::Dict{String,Axis}, extr::SoloExtract, info_solo)
    foreach(empty!, values(axs))
    essn = extr.essentials
    modl2d_norm = essn.modl2d |> m -> m ./ (sum(m) * (essn.step_modl / 2)^2)
    x, y = essn.smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x, y) .* essn.step_posi
    x_modl, y_modl = (x, y) .* essn.step_modl
    y_modl_sm = (0:1:essn.smwh[2]) * essn.step_modl
    clrmap = gen_clrmap_solo(hue_theme_istp[info_solo["istp"]])

    nvlp = extr.envelope
    shade_mainpeak = extr.fit_tailess["fitfn_main"](y_modl_sm)
    shade_peaks = extr.fit_tailess["fitfn"](y_modl_sm)
    band!(axs["upright"], y_modl_sm, 0, shade_mainpeak, color=(:gray, 0.1))
    band!(axs["upright"], y_modl_sm, shade_mainpeak, shade_peaks, color=(:darkseagreen1, 0.5))

    heatmap!(axs["dens"], x_posi, y_posi, essn.dens2d'; colorrange=(0, 16.0), colormap=clrmap, rasterize=true)
    draw_rotated_ellipse!(axs["dens"], nvlp["cent"], nvlp["size"], nvlp["rotation"]; color=(:darkseagreen1, 0.5))

    heatmap!(axs["modl"], y_modl_sm, x_modl, modl2d_norm[essn.smwh[2]+1:end, :]; colorrange=(0, 10.0), colormap=clrmap, rasterize=true)
    lines!(axs["upright"], y_modl_sm, essn.prfl_modl_norm_px[essn.smwh[2]+1:end], color=(:black, 0.4), linewidth=1)
    lines!(axs["sideway"], essn.prfl_modl_norm_px[essn.smwh[2]+1:end], y_modl_sm, color=(:black, 0.4), linewidth=1)
    lines!(axs["upright"], y_modl_sm, extr.prfl_modl_norm_tailess_px[essn.smwh[2]+1:end], color=:black, linewidth=1)
    lines!(axs["sideway"], extr.prfl_modl_norm_tailess_px[essn.smwh[2]+1:end], y_modl_sm, color=:black, linewidth=1)
    axs["sideway"].yreversed = true
    axs["sideway"] |> hidedecorations!
    axs["modl"] |> hidedecorations!
    axs["dens"] |> hidedecorations!
    axs["upright"].yticklabelsvisible = false
    axs["upright"].xticklabelsvisible = false
    axs["dens"].aspect = DataAspect()
    xlims!(axs["upright"], 0, 0.6)
    xlims!(axs["modl"], 0, 0.6)
    xlims!(axs["dens"], -5, 5)
    ylims!(axs["dens"], -10, 10)
    ylims!(axs["upright"], -0.2, 1.8)
    ylims!(axs["modl"], (-10.5, 10.5) .* essn.step_modl)
    ylims!(axs["sideway"], 0.15, 0.45)
    xlims!(axs["sideway"], 0.0, 1.5)
    vlines!(axs["modl"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.2))
    vlines!(axs["upright"], 0.3; color=RGBAf(Oklch(0.3, 0, 0), 0.4))
    hlines!(axs["upright"], 0.0; color=(:darkseagreen1, 0.5))
    vlines!(axs["upright"], extr.sidepeak["wavenum"]; color=(:mediumspringgreen, 1.0))
    vlines!(axs["sideway"], extr.sidepeak["height"]; color=(:mediumspringgreen, 1.0))
    mmt = extr.moments_modl
    errorbars!(axs["upright"], [mmt["wavenum"]], [1.5], [mmt["width"] / 2], [mmt["width"] / 2]; direction=:x, color=:sienna2, whiskerwidth=8)
    lines!(axs["sideway"], [mmt["height"], mmt["height"]], [0.2, 0.4]; color=(:sienna2, 1.0))
    band!(axs["sideway"], [0, mmt["height"]], [0.2, 0.2], [0.4, 0.4]; color=(:sienna2, 0.2))

    text!(axs["modl"], 0.35, -0.16; text="$(info_solo["t_hold"]) ms | rep $(info_solo["repeat"])", color=:black, strokewidth=0.5, strokecolor=:white, fontsize=16, align=(:center, :top))
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

function fit_prfl_modl_twinpeak_decay_1d(coor, prfl, mask)
    # parameters: [MainPeak.Height MainPeak.Width SidePeak.Height SidePeak.Width SidePeak.Pos Decay.Height Decay.Length]
    model(k, p) = p[1] .* exp.(-k .^ 2 ./ (2 .* p[2] .^ 2)) .+ p[3] .* exp.(-(k .- p[5]) .^ 2 ./ (2 .* p[4] .^ 2)) .+ p[6] .* exp.(-abs.(k) ./ p[7])
    p_init = [3.0, 0.1, 0.5, 0.05, 0.3, 0.5, 0.8]
    p_upper = [Inf, 0.30, 2.0, 0.100, 0.37, Inf, 5.0]
    p_lower = [2.0, 0.02, 0.0, 0.018, 0.23, 0.0, 0.5]
    fit = curve_fit(model, coor[mask], prfl[mask], p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    return Dict(
        "fit" => fit,
        "model" => model,
        "params" => params_fit,
        "tail" => k -> params_fit[6] .* exp.(-abs.(k) ./ params_fit[7])
    )
end

function fit_dens2d_gaussian_elliptic_disk(xs, ys, dens, mask)
    X = [x for y in ys, x in xs]
    Y = [y for y in ys, x in xs]
    xydata = hcat(vec(X[mask]), vec(Y[mask]))
    zdata = vec(dens[mask])
    model(coords, p) = begin
        x = coords[:, 1]
        y = coords[:, 2]
        A, x0, y0, σx, σy, θ = p
        c = cos(θ)
        s = sin(θ)
        dx = x .- x0
        dy = y .- y0
        xp = c .* dx .+ s .* dy
        yp = (-s) .* dx .+ c .* dy
        return A .* exp.(-(xp .^ 2 ./ (2σx^2) .+ yp .^ 2 ./ (2σy^2)))
    end
    params_init = Float64[10, 0, 0, 2, 5, -15/180*π]
    params_upper = Float64[25, 5, 10, 10, 20, 45/180*π]
    params_lower = Float64[0, -5, -10, 1, 2, -45/180*π]
    fit = curve_fit(
        model, xydata, zdata,
        params_init;
        lower=params_lower,
        upper=params_upper,
    )
    params_fit = coef(fit)
    fitfn(coords) = model(coords, params_fit)
    return Dict(
        "fit" => fit,
        "model" => model,
        "params" => params_fit,
        "fitfn" => fitfn
    )
end

function fit_prfl_modl_twinpeak_1d(coor, prfl, mask)
    # parameters: [MainPeak.Height MainPeak.Width SidePeak.Height SidePeak.Width SidePeak.Pos]
    model(k, p) = @. p[1] * exp(-k^2 / (2 * p[2]^2)) + p[3] * exp(-(k - p[5])^2 / (2 * p[4]^2))
    p_init = [3.0, 0.1, 0.5, 0.05, 0.3]
    p_upper = [Inf, 0.30, 2.0, 0.100, 0.37]
    p_lower = [2.0, 0.02, 0.0, 0.018, 0.23]
    fit = curve_fit(model, coor[mask], prfl[mask], p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    fitfn_main(k) = params_fit |> p -> (@. p[1] * exp(-k^2 / (2 * p[2]^2)))
    fitfn(k) = model(k, params_fit)
    return Dict(
        "fit" => fit,
        "model" => model,
        "params" => params_fit,
        "fitfn_main" => fitfn_main,
        "fitfn" => fitfn,
    )
end
