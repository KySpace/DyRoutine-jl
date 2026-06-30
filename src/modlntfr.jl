using FFTW
using ImageFiltering
using LinearAlgebra
using LsqFit
using Statistics

function orient_ntfr2d_axes(
    ntfr2d::AbstractArray{<:Real,4},
    x_dens::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
)
    n_x = length(x_dens)
    n_IB = length(val_IB)
    n_istp = length(val_istp)

    ntfr2d_fmt = Array{Matrix{Float64}}(undef, n_IB, n_istp)
    if size(ntfr2d) == (n_x, n_x, n_istp, n_IB)
        for idx_IB in 1:n_IB, idx_istp in 1:n_istp
            ntfr2d_fmt[idx_IB, idx_istp] = Float64.(@view ntfr2d[:, :, idx_istp, idx_IB])
        end
        return ntfr2d_fmt
    end
    if size(ntfr2d) == (n_IB, n_istp, n_x, n_x)
        for idx_IB in 1:n_IB, idx_istp in 1:n_istp
            ntfr2d_fmt[idx_IB, idx_istp] = Float64.(@view ntfr2d[idx_IB, idx_istp, :, :])
        end
        return ntfr2d_fmt
    end

    throw(DimensionMismatch(
        "ntfr2d_mean size $(size(ntfr2d)) must be either (x, y, istp, IB) " *
        "$((n_x, n_x, n_istp, n_IB)) or (IB, istp, x, y) $((n_IB, n_istp, n_x, n_x)).",
    ))
end

function orient_prfl_axes(
    prfl::AbstractArray{<:Real,3},
    x_modl::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
)
    n_x = length(x_modl)
    n_IB = length(val_IB)
    n_istp = length(val_istp)

    size(prfl) == (n_x, n_istp, n_IB) && return prfl
    size(prfl) == (n_IB, n_istp, n_x) && return permutedims(prfl, (3, 2, 1))

    throw(DimensionMismatch(
        "profile size $(size(prfl)) must be either (x_modl, istp, IB) " *
        "$((n_x, n_istp, n_IB)) or (IB, istp, x_modl) $((n_IB, n_istp, n_x)).",
    ))
end

function calc_grid_center_profile(dens2d::AbstractMatrix{<:Real})
    idx_center = cld(size(dens2d, 1), 2)
    return vec(@view dens2d[idx_center, :])
end

function calc_symmetric_grid_column_profile(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
)
    profile = calc_grid_center_profile(dens2d)
    idx_center = cld(length(x_dens), 2)
    x_half = abs.(x_dens[idx_center:end])
    profile_half = (profile[idx_center:end] .+ reverse(profile[1:idx_center])) ./ 2
    return x_half, profile_half
end

function gaussian_tail_1d_model(r, params)
    A_wide, sigma_wide = params[1:2]
    return @. A_wide * exp(-r^2 / (2 * sigma_wide^2))
end

function gaussian_narrow_1d_model(r, params)
    A_narrow, sigma_narrow = params
    return @. A_narrow * exp(-r^2 / (2 * sigma_narrow^2))
end

function lorentzian_narrow_1d_model(r, params)
    A_narrow, gamma_narrow = params
    return @. A_narrow / (1 + (r / gamma_narrow)^2)
end

function get_narrow_1d_model(model_center::Symbol)
    model_center == :gaussian && return gaussian_narrow_1d_model
    model_center == :lorentzian && return lorentzian_narrow_1d_model
    throw(ArgumentError("model_center must be :gaussian or :lorentzian, got $model_center."))
end

function fit_two_gaussian_1d_guess(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    r_tail_min::Real,
    model_center::Symbol=:gaussian,
)
    r_profile, profile = calc_symmetric_grid_column_profile(x_dens, dens2d)
    mask_tail = r_profile .> r_tail_min
    any(mask_tail) || throw(ArgumentError("No profile coordinates found above r_tail_min=$r_tail_min."))
    max_r = maximum(r_profile)
    step_r = minimum(diff(r_profile))

    y_tail = Float64.(profile[mask_tail])
    p_init_tail = [max(maximum(y_tail), eps(Float64)), max_r / 2]
    p_lower_tail = [0.0, step_r]
    p_upper_tail = [Inf, max_r * 3]
    fit_tail = curve_fit(
        gaussian_tail_1d_model,
        r_profile[mask_tail],
        y_tail,
        p_init_tail;
        lower=p_lower_tail,
        upper=p_upper_tail,
        maxIter=20_000,
    )
    params_tail = coef(fit_tail)
    wide = gaussian_tail_1d_model(r_profile, params_tail)

    y_narrow = max.(Float64.(profile) .- wide, 0.0)
    mask_narrow = r_profile .<= r_tail_min
    p_init_narrow = [max(maximum(y_narrow[mask_narrow]), eps(Float64)), max(r_tail_min / 4, step_r)]
    p_lower_narrow = [0.0, step_r]
    p_upper_narrow = [Inf, r_tail_min]
    fit_narrow = curve_fit(
        get_narrow_1d_model(model_center),
        r_profile[mask_narrow],
        y_narrow[mask_narrow],
        p_init_narrow;
        lower=p_lower_narrow,
        upper=p_upper_narrow,
        maxIter=20_000,
    )
    params_narrow = coef(fit_narrow)

    return (;
        A_narrow=params_narrow[1],
        sigma_narrow=params_narrow[2],
        A_wide=params_tail[1],
        sigma_wide=params_tail[2],
        B=0.0,
        r_profile,
        profile,
        wide,
        model_center,
    )
end

function two_gaussian_2d_model(coords, params)
    x0, y0, A_narrow, sigma_narrow, A_wide, sigma_wide = params[1:6]
    x = @view coords[1, :]
    y = @view coords[2, :]
    r2 = @. (x - x0)^2 + (y - y0)^2
    return @. A_narrow * exp(-r2 / (2 * sigma_narrow^2)) +
              A_wide * exp(-r2 / (2 * sigma_wide^2))
end

function lorentzian_gaussian_2d_model(coords, params)
    x0, y0, A_narrow, gamma_narrow, A_wide, sigma_wide = params[1:6]
    x = @view coords[1, :]
    y = @view coords[2, :]
    r2 = @. (x - x0)^2 + (y - y0)^2
    return @. A_narrow / (1 + r2 / gamma_narrow^2) +
              A_wide * exp(-r2 / (2 * sigma_wide^2))
end

function get_density_2d_model(model_center::Symbol)
    model_center == :gaussian && return two_gaussian_2d_model
    model_center == :lorentzian && return lorentzian_gaussian_2d_model
    throw(ArgumentError("model_center must be :gaussian or :lorentzian, got $model_center."))
end

function double_gaussian_disk_2d_model(coords, params)
    x0, y0, A_narrow, sigma_x_narrow, sigma_y_narrow, A_wide, sigma_wide = params
    x = @view coords[1, :]
    y = @view coords[2, :]
    dx2 = @. (x - x0)^2
    dy2 = @. (y - y0)^2
    r2 = @. dx2 + dy2
    return @. A_narrow * exp(-dx2 / (2 * sigma_x_narrow^2) - dy2 / (2 * sigma_y_narrow^2)) +
              A_wide * exp(-r2 / (2 * sigma_wide^2))
end

function double_gaussian_disk_2d_model_abrr(coords, params)
    length(params) == 8 || throw(ArgumentError(
        "double_gaussian_disk_2d_model_abrr expects 8 params " *
        "(x0, y0, A_narrow, sigma_x_narrow, sigma_y_narrow, A_wide, sigma_wide, beta), got $(length(params)).",
    ))
    x0, y0, A_narrow, sigma_x_narrow, sigma_y_narrow, A_wide, sigma_wide, beta = params
    x = @view coords[1, :]
    y = @view coords[2, :]
    dx2 = @. (x - x0)^2
    dy2 = @. (y - y0)^2
    r2 = @. dx2 + dy2
    narrow = @. A_narrow * exp(-dx2 / (2 * sigma_x_narrow^2) - dy2 / (2 * sigma_y_narrow^2))
    tail = @. A_wide * exp(-r2 / (2 * sigma_wide^2))
    return @. tail + narrow + beta * narrow^2
end

function log_double_gaussian_disk_2d_model(coords, params)
    return log.(max.(double_gaussian_disk_2d_model(coords, params), eps(Float64)))
end

function log_double_gaussian_disk_2d_model_abrr(coords, params)
    return log.(max.(double_gaussian_disk_2d_model_abrr(coords, params), eps(Float64)))
end

function calc_fit_coords_2d(x_dens::AbstractVector{<:Real}, stride::Integer)
    idx = 1:stride:length(x_dens)
    x_fit = Float64.(x_dens[idx])
    coords = Matrix{Float64}(undef, 2, length(x_fit)^2)
    coords[1, :] .= repeat(x_fit; outer=length(x_fit))
    coords[2, :] .= repeat(x_fit; inner=length(x_fit))
    return idx, coords
end

function fit_two_gaussian_2d(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    guess_1d;
    center_bound::Real,
    stride::Integer,
    threshold::Real,
    sigma_wide_min::Real,
    r_narrow_max::Real,
    maxiter::Integer=30_000,
    model_center::Symbol=:gaussian,
)
    idx_fit, coords = calc_fit_coords_2d(x_dens, stride)
    z_fit = vec(Float64.(@view dens2d[idx_fit, idx_fit]))
    mask_fit = @. isfinite(z_fit) & (z_fit > threshold)
    any(mask_fit) || throw(ArgumentError("No 2D density points above threshold=$threshold."))
    coords_fit = @view coords[:, mask_fit]
    z_fit_sel = z_fit[mask_fit]
    step_x = minimum(diff(x_dens))
    max_x = maximum(abs, x_dens)
    p_init = [
        0.0,
        0.0,
        max(guess_1d.A_narrow, eps(Float64)),
        max(guess_1d.sigma_narrow, step_x),
        max(guess_1d.sigma_narrow, step_x),
        max(guess_1d.A_wide, eps(Float64)),
        max(guess_1d.sigma_wide, sigma_wide_min),
        0.01,
    ]
    p_lower = [-center_bound, -center_bound, eps(Float64), step_x, step_x, eps(Float64), sigma_wide_min, 0.0]
    p_upper = [center_bound, center_bound, Inf, r_narrow_max, r_narrow_max, Inf, max_x * 3, 2.0]
    fit = curve_fit(
        double_gaussian_disk_2d_model_abrr,
        coords_fit,
        z_fit_sel,
        p_init;
        lower=p_lower,
        upper=p_upper,
        maxIter=maxiter,
    )
    params = coef(fit)
    rss_rel = norm(residuals(fit)) / max(norm(z_fit_sel), eps(Float64))
    maxiter_reached = !fit.converged
    return (; params, rss_rel, maxiter_reached, guess_1d, model_center, threshold)
end

function interp2_bilinear(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    xq::Real,
    yq::Real,
)
    (xq < first(x_dens) || xq > last(x_dens) || yq < first(x_dens) || yq > last(x_dens)) && return NaN
    idx_x_hi = searchsortedfirst(x_dens, xq)
    idx_y_hi = searchsortedfirst(x_dens, yq)
    idx_x_hi == 1 && (idx_x_hi = 2)
    idx_y_hi == 1 && (idx_y_hi = 2)
    idx_x_hi > length(x_dens) && (idx_x_hi = length(x_dens))
    idx_y_hi > length(x_dens) && (idx_y_hi = length(x_dens))
    idx_x_lo = idx_x_hi - 1
    idx_y_lo = idx_y_hi - 1
    tx = (xq - x_dens[idx_x_lo]) / (x_dens[idx_x_hi] - x_dens[idx_x_lo])
    ty = (yq - x_dens[idx_y_lo]) / (x_dens[idx_y_hi] - x_dens[idx_y_lo])
    z00 = dens2d[idx_x_lo, idx_y_lo]
    z10 = dens2d[idx_x_hi, idx_y_lo]
    z01 = dens2d[idx_x_lo, idx_y_hi]
    z11 = dens2d[idx_x_hi, idx_y_hi]
    return (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 + (1 - tx) * ty * z01 + tx * ty * z11
end

function calc_centered_cross_profile(
    x_dens::AbstractVector{<:Real},
    dens2d::AbstractMatrix{<:Real},
    fit_density;
    axis::Symbol,
)
    x0, y0, A_narrow, sigma_x_narrow, sigma_y_narrow, A_wide, sigma_wide, beta = fit_density.params
    s_profile = Float64.(x_dens)
    profile =
        axis == :column ? [interp2_bilinear(x_dens, dens2d, x0, y0 + s) for s in s_profile] :
        axis == :row ? [interp2_bilinear(x_dens, dens2d, x0 + s, y0) for s in s_profile] :
        throw(ArgumentError("axis must be :column or :row, got $axis."))
    if axis == :column
        dx2 = zero.(s_profile)
        dy2 = @. s_profile^2
    elseif axis == :row
        dx2 = @. s_profile^2
        dy2 = zero.(s_profile)
    else
        throw(ArgumentError("axis must be :column or :row, got $axis."))
    end
    r2 = @. dx2 + dy2
    tail_raw = @. A_wide * exp(-r2 / (2 * sigma_wide^2))
    narrow_raw = @. A_narrow * exp(-dx2 / (2 * sigma_x_narrow^2) - dy2 / (2 * sigma_y_narrow^2))
    tail = tail_raw
    narrow = @. narrow_raw + beta * narrow_raw^2
    narrow_abrr = narrow
    tailess = profile .- tail
    return (; axis, s_profile, profile, tail, narrow, narrow_abrr, tailess, tail_raw, narrow_raw, beta, fit_density)
end

function fit_centered_density_profiles(
    x_dens::AbstractVector{<:Real},
    ntfr2d_fmt::AbstractMatrix{<:AbstractMatrix},
    r_tail_min::Real;
    center_bound::Real,
    stride::Integer,
    threshold::Real,
    sigma_wide_min::Real,
    maxiter::Integer=30_000,
    model_center::Symbol=:gaussian,
    log_tag=nothing,
    val_IB=nothing,
    val_istp=nothing,
)
    get_narrow_1d_model(model_center)
    get_density_2d_model(model_center)
    fit_density = Array{NamedTuple}(undef, size(ntfr2d_fmt))
    profile_column = Array{NamedTuple}(undef, size(ntfr2d_fmt))
    profile_row = similar(profile_column)

    for idx_IB in axes(ntfr2d_fmt, 1), idx_istp in axes(ntfr2d_fmt, 2)
        if !isnothing(log_tag)
            label_IB = isnothing(val_IB) ? idx_IB : val_IB[idx_IB]
            label_istp = isnothing(val_istp) ? idx_istp : val_istp[idx_istp]
            println("  [$log_tag] fitting 2D density IB_idx=$idx_IB IB=$label_IB istp=$label_istp")
            flush(stdout)
        end
        dens2d = ntfr2d_fmt[idx_IB, idx_istp]
        guess_1d = fit_two_gaussian_1d_guess(x_dens, dens2d, r_tail_min, model_center)
        fit_2d = fit_two_gaussian_2d(
            x_dens,
            dens2d,
            guess_1d;
            center_bound,
            stride,
            threshold,
            sigma_wide_min,
            r_narrow_max=r_tail_min,
            maxiter,
            model_center,
        )
        if !isnothing(log_tag)
            params = fit_2d.params
            println(
                "  [$log_tag] fit done IB_idx=$idx_IB istp_idx=$idx_istp " *
                "rss=$(round(fit_2d.rss_rel; digits=4)) " *
                (fit_2d.maxiter_reached ? "maxiter=true " : "") *
                "β=$(round(params[8]; digits=4))",
            )
            flush(stdout)
        end
        fit_density[idx_IB, idx_istp] = fit_2d
        profile_column[idx_IB, idx_istp] = calc_centered_cross_profile(x_dens, dens2d, fit_2d; axis=:column)
        profile_row[idx_IB, idx_istp] = calc_centered_cross_profile(x_dens, dens2d, fit_2d; axis=:row)
    end

    return (; fit_density, profile_fits=(profile_column, profile_row))
end

function tucky1d(sml; alpha=0.2)
    [abs(x) - 1 + alpha |> x -> x < 0 ? 1.0 : (1 + cos(pi * x / alpha)) / 2
     for x in (-sml:sml) ./ sml]
end

function calc_prfl_modl_1d(dens::AbstractMatrix{<:Real}, smwh::Tuple{<:Integer,<:Integer}; step_modl::Real=1)
    smwh .* 2 .+ 1 == size(dens) || throw(DimensionMismatch("smwh $smwh expects density size $(smwh .* 2 .+ 1), got $(size(dens))."))
    smw, smh = smwh
    smw_dens_strip = 20
    smh_modl = 120
    tucky = tucky1d(smh; alpha=0.2)
    idx_strip = (smw_dens_strip |> s -> (-s:1:s) .+ smw .+ 1)
    idx_modl = (smh_modl |> s -> (-s:1:s) .+ smh .+ 1)
    dens_strip = @view dens[:, idx_strip]
    dens_mean = imfilter(dens_strip, Kernel.gaussian(2.5)) |> ds -> vec(mean(ds; dims=2))
    modl_full = abs.(fftshift(fft(dens_mean .* tucky)))
    modl = modl_full[idx_modl]
    return modl / (sum(modl) * step_modl / 2)
end

function calc_reconstructed_ntfr2d(
    x_dens::AbstractVector{<:Real},
    fit_density::AbstractMatrix,
)
    coords = Matrix{Float64}(undef, 2, length(x_dens)^2)
    x_grid = Float64.(x_dens)
    coords[1, :] .= repeat(x_grid; outer=length(x_grid))
    coords[2, :] .= repeat(x_grid; inner=length(x_grid))

    return map(fit_density) do fit
        reshape(
            double_gaussian_disk_2d_model_abrr(coords, fit.params),
            length(x_dens),
            length(x_dens),
        )
    end
end

function calc_reconstructed_prfl_modl(
    ntfr2d_fit::AbstractMatrix{<:AbstractMatrix},
    smwh::Tuple{<:Integer,<:Integer};
    step_modl::Real,
)
    return map(dens2d -> calc_prfl_modl_1d(dens2d, smwh; step_modl), ntfr2d_fit)
end

function fit_reconstructed_prfl_modl(
    x_dens::AbstractVector{<:Real},
    ntfr2d_fmt::AbstractMatrix{<:AbstractMatrix},
    smwh::Tuple{<:Integer,<:Integer};
    step_modl::Real,
    r_tail_min::Real,
    center_bound::Real,
    stride::Integer,
    threshold::Real,
    sigma_wide_min::Real,
    maxiter::Integer=30_000,
    model_center::Symbol=:gaussian,
)
    fit_centered = fit_centered_density_profiles(
        x_dens,
        ntfr2d_fmt,
        r_tail_min;
        center_bound,
        stride,
        threshold,
        sigma_wide_min,
        maxiter,
        model_center,
    )
    ntfr2d_fit = calc_reconstructed_ntfr2d(x_dens, fit_centered.fit_density)
    prfl_modl_fit = calc_reconstructed_prfl_modl(ntfr2d_fit, smwh; step_modl)
    return (; fit_centered, ntfr2d_fit, prfl_modl_fit)
end

function pack_prfl_modl_fit(
    prfl_modl_fit::AbstractMatrix{<:AbstractVector},
    x_modl::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
)
    size(prfl_modl_fit) == (length(val_IB), length(val_istp)) || throw(DimensionMismatch(
        "prfl_modl_fit size $(size(prfl_modl_fit)) must match (IB, istp) $((length(val_IB), length(val_istp))).",
    ))
    n_modl = length(prfl_modl_fit[1])
    length(x_modl) == n_modl || throw(DimensionMismatch(
        "x_modl length $(length(x_modl)) must match reconstructed profile length $n_modl.",
    ))
    prfl_modl = Array{Float64}(undef, n_modl, length(val_istp), length(val_IB))
    for idx_IB in axes(prfl_modl_fit, 1), idx_istp in axes(prfl_modl_fit, 2)
        length(prfl_modl_fit[idx_IB, idx_istp]) == n_modl || throw(DimensionMismatch(
            "prfl_modl_fit[$idx_IB, $idx_istp] length $(length(prfl_modl_fit[idx_IB, idx_istp])) must match $n_modl.",
        ))
        prfl_modl[:, idx_istp, idx_IB] .= prfl_modl_fit[idx_IB, idx_istp]
    end
    return prfl_modl
end
