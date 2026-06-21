# Pure fit model evaluators. Fitting helpers store only numeric params and
# callers reconstruct curves through these named functions.

function fit_prfl_modl_twinpeak_decay_1d_model(k, params)
    M, σ0, P, σ, p, D, λ = params
    return @. M * exp(-k^2 / (2 * σ0^2)) + P * exp(-(k - p)^2 / (2 * σ^2)) + D * exp(-abs(k) / λ)
end

function fit_prfl_modl_twinpeak_decay_1d_tail(k, params)
    _, _, _, _, _, D, λ = params
    return @. D * exp(-abs(k) / λ)
end

function fit_prfl_modl_sidepeak_decay_1d_model(k, params)
    P, σ, p, D, λ = params
    return @. P * exp(-(k - p)^2 / (2 * σ^2)) + D * exp(-abs(k) / λ)
end

function fit_prfl_modl_sidepeak_decay_1d_tail(k, params)
    _, _, _, D, λ = params
    return @. D * exp(-abs(k) / λ)
end

function fit_prfl_modl_twinpeak_1d_model(k, params)
    M, σ0, P, σ, p = params
    return @. M * exp(-k^2 / (2 * σ0^2)) + P * exp(-(k - p)^2 / (2 * σ^2))
end

function fit_prfl_modl_twinpeak_1d_main(k, params)
    M, σ0, _, _, _ = params
    return @. M * exp(-k^2 / (2 * σ0^2))
end

function fit_prfl_modl_sidepeak_1d_model(k, params)
    P, σ, p = params
    return @. P * exp(-(k - p)^2 / (2 * σ^2))
end

function fit_evol_oscillation_decay_model(t, params)
    A, C, λ, ν, φ = params
    scaling = 1000
    return @. A * exp(- t / λ) * cos(2π * ν * t / scaling + φ) + C
end

function fit_dens2d_gaussian_elliptic_disk_model(coords, params)
    x = coords[:, 1]
    y = coords[:, 2]
    A, x0, y0, σx, σy, θ = params
    c = cos(θ)
    s = sin(θ)
    dx = x .- x0
    dy = y .- y0
    xp = c .* dx .+ s .* dy
    yp = (-s) .* dx .+ c .* dy
    return @. A * exp(-(xp^2 / (2σx^2) + yp^2 / (2σy^2)))
end

function fit_dens2d_gaussian_round_disk_model(coords, params)
    x = coords[:, 1]
    y = coords[:, 2]
    A, x0, y0, σ = params
    dx = x .- x0
    dy = y .- y0
    return @. A * exp(-(dx^2 / (2σ^2) + dy^2 / (2σ^2)))
end
