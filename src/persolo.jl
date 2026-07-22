# Analysis per-solo image
# When each shot contains one image, it is per-shot
# can also be applied to composite images from multiple solos, like stacked images
using LsqFit
using NaNStatistics: movmean
using FFTW
using DSP.Windows: hanning
using NumericalIntegration: integrate
using StatsBase: geomean

isdefined(Main, :fit_prfl_modl_sidepeak_1d_model) || include(joinpath(@__DIR__, "fitmodels.jl"))

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
    arr::AbstractMatrix,
    smwh::Tuple{Integer,Integer};
    len_avg::Integer=10,
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

function calc_prfl_moment(coor, prfl)
    @assert length(coor) == length(prfl) "coordinate and profile length mismatch"
    prfl = clamp.(prfl, 0, Inf)
    ntgr_over_coor = y -> integrate(coor, y)
    weight = prfl |> ntgr_over_coor
    height = weight / (coor[end] - coor[1])
    expval = coor .* prfl |> ntgr_over_coor |> u -> u ./ weight
    var = (coor .- expval) .^ 2 .* prfl |> ntgr_over_coor |> sqrt |> u -> u / sqrt(weight)
    return (; weight, wavenum=expval, width=var, height, coor)
end

function fit_prfl_modl_twinpeak_decay_1d(
    coor, prfl, mask;
    M_hint=(max=Inf, min=2.0, init=3.0),
    σ0_hint=(max=0.30, min=0.02, init=0.1),
    P_hint=(max=2.0, min=0.0, init=0.5),
    σ_hint=(max=0.100, min=0.018, init=0.05),
    p_hint=(max=0.37, min=0.23, init=0.3),
    D_hint=(max=Inf, min=0.0, init=0.5),
    λ_hint=(max=5.0, min=0.5, init=0.8),
)
    # parameters: [MainPeak.Height MainPeak.Width SidePeak.Height SidePeak.Width SidePeak.Pos Decay.Height Decay.Length]
    p_init = Float64[M_hint.init, σ0_hint.init, P_hint.init, σ_hint.init, p_hint.init, D_hint.init, λ_hint.init]
    p_upper = Float64[M_hint.max, σ0_hint.max, P_hint.max, σ_hint.max, p_hint.max, D_hint.max, λ_hint.max]
    p_lower = Float64[M_hint.min, σ0_hint.min, P_hint.min, σ_hint.min, p_hint.min, D_hint.min, λ_hint.min]
    fit = curve_fit(fit_prfl_modl_twinpeak_decay_1d_model, coor[mask], prfl[mask], p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    rss_rel = (fit |> residuals |> r -> sqrt(sum(abs2, r))) / (prfl[mask] |> d -> sqrt(sum(abs2, d)))
    return (; params=params_fit, rss_rel)
end

function fit_prfl_modl_sidepeak_decay_1d(
    coor, prfl, mask;
    P_hint=(max=2.0, min=0.0, init=0.5),
    σ_hint=(max=0.100, min=0.018, init=0.05),
    p_hint=(max=0.37, min=0.23, init=0.3),
    D_hint=(max=Inf, min=0.0, init=0.5),
    λ_hint=(max=5.0, min=0.5, init=0.8),
)
    # parameters: [SidePeak.Height SidePeak.Width SidePeak.Pos Decay.Height Decay.Length]
    p_init = Float64[P_hint.init, σ_hint.init, p_hint.init, D_hint.init, λ_hint.init]
    p_upper = Float64[P_hint.max, σ_hint.max, p_hint.max, D_hint.max, λ_hint.max]
    p_lower = Float64[P_hint.min, σ_hint.min, p_hint.min, D_hint.min, λ_hint.min]
    fit = curve_fit(fit_prfl_modl_sidepeak_decay_1d_model, coor[mask], prfl[mask], p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    rss_rel = (fit |> residuals |> r -> sqrt(sum(abs2, r))) / (prfl[mask] |> d -> sqrt(sum(abs2, d)))
    return (; params=params_fit, rss_rel)
end

function copy_symmetric_2d(mask::AbstractMatrix{Bool})
    size(mask, 1) |> isodd || throw(ArgumentError("mask height must be odd, got $(size(mask, 1))."))
    return mask .| reverse(mask; dims=(1, 2))
end

function copy_symmetric_2d(modl::AbstractMatrix{<:Real})
    size(modl, 1) |> isodd || throw(ArgumentError("modl height must be odd, got $(size(modl, 1))."))
    pos_cent = (size(modl, 1) + 1) ÷ 2
    modl_symm = modl .+ reverse(modl; dims=(1, 2))
    modl_symm[pos_cent, :] ./= 2
    return modl_symm
end

function calc_masked_prfl_modl(modl2d::AbstractMatrix{<:Real}, mask::AbstractMatrix{Bool})
    size(modl2d) == size(mask) || throw(DimensionMismatch("modl2d size $(size(modl2d)) does not match mask size $(size(mask))."))
    sum_mask = sum(Int.(mask); dims=2)
    prfl_sum = sum(modl2d .* mask; dims=2)
    return map((val, n) -> n == 0 ? 0.0 : val / n, vec(prfl_sum), vec(sum_mask))
end

function build_modl_mask(mask_spec, x_modl::AbstractVector, y_modl::AbstractVector)
    mask =
        if mask_spec isa AbstractMatrix{Bool}
            mask_spec
        elseif mask_spec isa Function
            [mask_spec(x, y) for y in y_modl, x in x_modl]
        else
            throw(ArgumentError("modulation mask entries must be Bool matrices or functions of (x, y); got $(typeof(mask_spec))."))
        end
    size(mask) == (length(y_modl), length(x_modl)) ||
        throw(DimensionMismatch("mask size $(size(mask)) must match modulation grid $((length(y_modl), length(x_modl)))."))
    return copy_symmetric_2d(mask)
end

function build_modl_masks(mask_modl::NamedTuple, x_modl::AbstractVector, y_modl::AbstractVector)
    required = (:main, :sidepeak)
    all(hasproperty(mask_modl, key) for key in required) ||
        throw(ArgumentError("mask_modl must include at least :main and :sidepeak masks."))
    mask_base = build_modl_mask(mask_modl.sidepeak, x_modl, y_modl)
    mask_main_base = build_modl_mask(mask_modl.main, x_modl, y_modl)
    mask_fringe = hasproperty(mask_modl, :fringe) ? build_modl_mask(mask_modl.fringe, x_modl, y_modl) : falses(length(y_modl), length(x_modl))
    mask_center = hasproperty(mask_modl, :center) ? build_modl_mask(mask_modl.center, x_modl, y_modl) : falses(length(y_modl), length(x_modl))
    mask_side = @. mask_base & !(mask_fringe | mask_center)
    mask_main = @. mask_main_base & !mask_fringe
    return (; fringe=mask_fringe, center=mask_center, side=mask_side, main=mask_main)
end

function build_default_modl_masks(x_modl::AbstractVector, y_modl::AbstractVector)
    mask_all = trues(length(y_modl), length(x_modl))
    mask_none = falses(length(y_modl), length(x_modl))
    return (; fringe=mask_none, center=mask_none, side=mask_all, main=mask_all)
end

function calc_prfl_norm_px_masked(modl2d::AbstractMatrix{<:Real}, masks::NamedTuple, step_modl::Tuple{<:Real,<:Real})
    main_raw = calc_masked_prfl_modl(modl2d, masks.main)
    side_raw = calc_masked_prfl_modl(modl2d, masks.side)
    norm_main = sum(main_raw) * step_modl[2] / 2
    norm_main > 0 || throw(ArgumentError("masked main modulation profile has nonpositive normalization $norm_main."))
    return (;
        main=(; norm=norm_main, raw=main_raw, normed_px=main_raw ./ norm_main),
        side=(; raw=side_raw, normed_px=side_raw ./ norm_main),
        mask=masks,
    )
end

function fit_dens2d_gaussian_elliptic_disk(
    xs, ys, dens, mask;
    θ_hint=(max=20.0 / 180 * π, min=-10.0 / 180 * π, init=10.0 / 180 * π),
    A_hint=(max=25.0, min=0, init=10.0),
    preprocess = (x -> x),
)
    X = [x for y in ys, x in xs]
    Y = [y for y in ys, x in xs]
    xydata = hcat(vec(X[mask]), vec(Y[mask]))
    zdata = dens |> preprocess |> ds -> ds[mask] |> vec
    x_min, x_max = [minimum(xs), maximum(xs)]
    y_min, y_max = [minimum(ys), maximum(ys)]
    x_mid, y_mid = [(x_min + x_max) / 2, (y_min + y_max) / 2]
    x_scale, y_scale = [(x_max - x_min) / 2, (y_max - y_min) / 2] ./ 3
    params_init = Float64[A_hint.init, x_mid, y_mid, x_scale, y_scale, θ_hint.init]
    params_upper = Float64[A_hint.max, x_max, y_max, x_scale*10, y_scale*10, θ_hint.max]
    params_lower = Float64[A_hint.min, x_min, y_min, x_scale/10, y_scale/10, θ_hint.min]
    fit = curve_fit(
        fit_dens2d_gaussian_elliptic_disk_model, xydata, zdata,
        params_init;
        lower=params_lower,
        upper=params_upper,
    )
    params_fit = coef(fit)
    rss_rel = (fit |> residuals |> r -> sqrt(sum(abs2, r))) / (dens[mask] |> d -> sqrt(sum(abs2, d)))
    return (; params=params_fit, rss_rel)
end

function fit_dens2d_gaussian_round_disk(
    xs, ys, dens, mask;
    A_hint=(max=25.0, min=0, init=10.0),
)
    X = [x for y in ys, x in xs]
    Y = [y for y in ys, x in xs]
    xydata = hcat(vec(X[mask]), vec(Y[mask]))
    zdata = vec(dens[mask])
    x_min, x_max = [minimum(xs), maximum(xs)]
    y_min, y_max = [minimum(ys), maximum(ys)]
    x_mid, y_mid = [(x_min + x_max) / 2, (y_min + y_max) / 2]
    scale = [(x_max - x_min) / 2, (y_max - y_min) / 2] ./ 3 |> geomean
    params_init = Float64[A_hint.init, x_mid, y_mid, scale]
    params_upper = Float64[A_hint.max, x_max, y_max, scale*10]
    params_lower = Float64[A_hint.min, x_min, y_min, scale/10]
    fit = curve_fit(
        fit_dens2d_gaussian_round_disk_model, xydata, zdata,
        params_init;
        lower=params_lower,
        upper=params_upper,
    )
    params_fit = coef(fit)
    rss_rel = (fit |> residuals |> r -> sqrt(sum(abs2, r))) / (dens[mask] |> d -> sqrt(sum(abs2, d)))
    return (; params=params_fit, rss_rel)
end

function fit_prfl_modl_twinpeak_1d(
    coor, prfl, mask;
    M_hint=(max=Inf, min=2.0, init=3.0),
    σ0_hint=(max=0.30, min=0.05, init=0.1),
    P_hint=(max=10.0, min=0.0, init=0.5),
    σ_hint=(max=0.200, min=0.018, init=0.05),
    p_hint=(max=0.37, min=0.23, init=0.3),
)
    # parameters: [MainPeak.Height MainPeak.Width SidePeak.Height SidePeak.Width SidePeak.Pos]
    p_init = Float64[M_hint.init, σ0_hint.init, P_hint.init, σ_hint.init, p_hint.init]
    p_upper = Float64[M_hint.max, σ0_hint.max, P_hint.max, σ_hint.max, p_hint.max]
    p_lower = Float64[M_hint.min, σ0_hint.min, P_hint.min, σ_hint.min, p_hint.min]
    fit = curve_fit(fit_prfl_modl_twinpeak_1d_model, coor[mask], prfl[mask], p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    rss_rel = (fit |> residuals |> r -> sqrt(sum(abs2, r))) / (prfl[mask] |> d -> sqrt(sum(abs2, d)))
    return (; params=params_fit, rss_rel)
end

function fit_prfl_modl_sidepeak_1d(
    coor, prfl, mask;
    P_hint=(max=10.0, min=0.0, init=0.5),
    σ_hint=(max=0.200, min=0.018, init=0.05),
    p_hint=(max=0.37, min=0.23, init=0.3),
)
    # parameters: [SidePeak.Height SidePeak.Width SidePeak.Pos]
    p_init = Float64[P_hint.init, σ_hint.init, p_hint.init]
    p_upper = Float64[P_hint.max, σ_hint.max, p_hint.max]
    p_lower = Float64[P_hint.min, σ_hint.min, p_hint.min]
    fit = curve_fit(fit_prfl_modl_sidepeak_1d_model, coor[mask], prfl[mask], p_init; lower=p_lower, upper=p_upper)
    params_fit = coef(fit)
    rss_rel = (fit |> residuals |> r -> sqrt(sum(abs2, r))) / (prfl[mask] |> d -> sqrt(sum(abs2, d)))
    return (; params=params_fit, rss_rel)
end

struct SoloEssentials
    dens2d::AbstractMatrix
    modl2d::AbstractMatrix
    dens2d_core::AbstractMatrix
    offset_cent_core::Tuple{<:Real,<:Real}
    smwh_core::Tuple{<:Real,<:Real}
    prfl_strip::AbstractVector
    prfl_modl::NamedTuple
    prfls_core::NamedTuple
    smwh::Tuple{<:Real,<:Real}
    smwh_strip::Tuple{<:Real,<:Real}
    step_posi::Tuple{Real,Real}
    step_modl::Tuple{Real,Real}
    sum_dens_full::Real
end

struct SoloSidepeak
    prfl_norm_tailess_px::AbstractVector
    params_tailess::NamedTuple
    fit_tailess::NamedTuple
    moments::NamedTuple
end

struct SoloEnvelope
    fit_asymm_2d::NamedTuple
    params_asymm::NamedTuple
    fit_round_2d::NamedTuple
    params_round::NamedTuple
end

struct SoloExtract
    essentials::SoloEssentials
    sidepeak::Union{SoloSidepeak,Nothing}
    envelope::Union{SoloEnvelope,Nothing}
end

function calc_solo_essn_2d(
    dens::AbstractMatrix,
    cent::Tuple{<:Real,<:Real},
    smwh::Tuple{<:Real,<:Real},
    px_in_um::Union{Real,Tuple{<:Real,<:Real}},
    cent_core::Tuple{<:Real,<:Real},
    smwh_core::Tuple{<:Real,<:Real};
    smwh_strip::Tuple{<:Real,<:Real}=smwh,
    mask_modl::NamedTuple=NamedTuple(),
)
    px_in_um = length(px_in_um) == 1 ? (px_in_um, px_in_um) : px_in_um
    dens_roi = crop_center(dens, cent, smwh) |> copy
    sum_dens = sum(dens)
    step_posi = px_in_um
    step_modl = 1 ./ (2 .* smwh_core .* px_in_um)
    x_posi, y_posi = map(u -> (-u:1:u), smwh) .* step_posi
    dens2d_core = crop_center(dens, cent_core, smwh_core) |> copy
    prfl_strip = crop_center(dens, cent_core, smwh_strip) |> m -> mean(m, dims=2) |> vec
    modl_roi = dens2d_core .* gen_win_hann_2d(smwh_core) |> fft |> fftshift |> c -> abs.(c)
    x_modl, y_modl = smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_modl
    masks = isempty(mask_modl) ? build_default_modl_masks(x_modl, y_modl) : build_modl_masks(mask_modl, x_modl, y_modl)
    prfls_modl = calc_prfl_norm_px_masked(modl_roi, masks, step_modl)
    prfls_core = dens2d_core |> ds -> (; axial = vec(mean(ds; dims=2)), radial = vec(mean(ds; dims=1)))
    return SoloEssentials(
        dens_roi,
        modl_roi,
        dens2d_core,
        (x_posi[cent_core[1]], y_posi[cent_core[2]]),
        smwh_core,
        prfl_strip,
        prfls_modl,
        prfls_core,
        smwh,
        smwh_strip,
        step_posi,
        step_modl,
        sum_dens
    )
end

function calc_solo_extr(
    essn::SoloEssentials,
    fit_stack::Union{NamedTuple,Nothing};
    proc_sidepeak::Bool=false,
    proc_envelope::Bool=false,
    selector_moment::Function=y -> (y .> 0.10) .& (y .< 0.50),
    selector_sidepeak::Function=y -> (y .> 0.1) .& (y .< 0.5),
    fit_tailess_kwargs::NamedTuple=NamedTuple(),
    fit_asymm_kwargs::NamedTuple=NamedTuple(),
    fit_round_kwargs::NamedTuple=NamedTuple(),
    fit_stack_tail_model::Function=fit_prfl_modl_sidepeak_decay_1d_tail,
)
    x_modl, y_modl = essn.smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_modl
    x_posi, y_posi = essn.smwh |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* essn.step_posi
    sel_sidepeak = selector_sidepeak(y_modl)
    sidepeak = proc_sidepeak ?
    begin
        mask_mmt = selector_moment(y_modl)
        prfl_tailess = essn.prfl_modl.side.normed_px - fit_stack_tail_model(y_modl, fit_stack.params)
        fit_tailess = fit_prfl_modl_sidepeak_1d(y_modl, prfl_tailess, sel_sidepeak; fit_tailess_kwargs...)
        params_tailess = (;
            height=fit_tailess.params[1],
            width=fit_tailess.params[2],
            wavenum=fit_tailess.params[3],
            weight=sqrt(2 * pi) * fit_tailess.params[1] * fit_tailess.params[2],
            rel_residue=fit_tailess.rss_rel,
        )
        moments = calc_prfl_moment(y_modl[mask_mmt], prfl_tailess[mask_mmt])
        SoloSidepeak(prfl_tailess, params_tailess, fit_tailess, moments)
    end : nothing
    envelope = proc_envelope ?
    begin
        fit_asymm = fit_dens2d_gaussian_elliptic_disk(x_posi, y_posi, essn.dens2d, :; fit_asymm_kwargs...)
        params_asymm = (;
            max=fit_asymm.params[1],
            cent=(fit_asymm.params[2], fit_asymm.params[3]),
            size=(fit_asymm.params[4], fit_asymm.params[5]),
            rotation=fit_asymm.params[6],
            rel_residue=fit_asymm.rss_rel,
        )
        fit_round = fit_dens2d_gaussian_round_disk(x_posi, y_posi, essn.dens2d, :; fit_round_kwargs...)
        params_round = (;
            max=fit_round.params[1],
            cent=(fit_round.params[2], fit_round.params[3]),
            size=fit_round.params[4],
            rel_residue=fit_round.rss_rel,
        )
        SoloEnvelope(fit_asymm, params_asymm, fit_round, params_round)
    end : nothing
    return SoloExtract(essn, sidepeak, envelope)
end
