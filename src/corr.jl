using MultivariateStats: PCA, fit, predict, projection
using Statistics: mean
using Peaks
using LsqFit: curve_fit, coef, residuals

struct ModeWeight{TProfile<:AbstractArray,TWeight<:AbstractArray}
    profile::TProfile
    weight::TWeight
end

function get_spectrum_peaks(freq, spct; min_prom=0.2)
    pks = spct |> findmaxima |> peakproms!(; min=min_prom) |> peakwidths!
    height_max = sum(pks.heights)
    pks_record = map(
        (idx, height) -> (freq=freq[idx], value=height, value_reduced=height / height_max), pks.indices, pks.heights
    ) |> p -> sort(p; by=x -> x.value_reduced, rev=true)
    return pks_record
end

function build_pca_matrix(samples::AbstractArray)
    n_sample = length(samples)
    n_sample > 0 || throw(ArgumentError("samples must contain at least one sample array."))
    sample_first = first(samples)
    sample_flat, rebuilder = flatten_with_rebuilder(sample_first)
    n_feature = length(sample_flat)
    n_feature > 0 || throw(ArgumentError("sample arrays must contain at least one value."))
    mat_sample = Matrix{Float64}(undef, n_feature, n_sample)
    for (idx_sample, sample) in enumerate(samples)
        sample_flat, _ = flatten_with_rebuilder(sample)
        mat_sample[:, idx_sample] .= sample_flat
    end
    return mat_sample, rebuilder
end

"""
`function fit_pca_modes(n_mode::Int, samples::AbstractArray)`
    samples can have an outer dimension (nd array) and an inner dimension (inner being an array or an nd array of arrays),
     - the outer dimension is where the variation is queried
     - the inner dimension is bunched together as one vector
    returns a vector of modes, where
     - the profile is reshaped into the format of the inner dimension
     - the weight is reshaped into the format of the outer dimension
"""
function fit_pca_modes(n_mode::Int, samples::AbstractArray)
    n_mode > 0 || throw(ArgumentError("n_mode=$n_mode must be positive."))

    mat_samples, rebuilder = build_pca_matrix(samples)
    n_feature, n_sample = size(mat_samples)
    n_sample > 1 || throw(ArgumentError("PCA requires at least two samples; got n_sample=$n_sample."))
    n_mode <= min(n_feature, n_sample - 1) || throw(ArgumentError("n_mode=$n_mode exceeds min(n_feature=$n_feature, n_sample - 1=$(n_sample - 1))."))

    pca_fit_lin = fit(PCA, mat_samples; maxoutdim=n_mode, pratio=1.0)
    mat_profile_lin = projection(pca_fit_lin) # (fmt of flattened sample, n_mode)
    mat_weight_lin = predict(pca_fit_lin, mat_samples) # (n_mode, n_sample)

    # returns a list of modes
    return [
        ModeWeight(
            copy(@view mat_profile_lin[:, idx_mode]) |> rebuilder,
            reshape(copy(@view mat_weight_lin[idx_mode, :]), size(samples)),
        )
        for idx_mode in 1:n_mode
    ]
end

function calc_spct_rep_evol(evols::AbstractVector{<:AbstractVector}, val_t::AbstractVector, freq_query::AbstractVector; sel_evol::Function=(_ -> true))
    n_rep = length(evols)
    mask_evol = map(sel_evol, val_t)
    val_t_sel = val_t[mask_evol]
    evol_mean = mean(evols)
    spct_mean_full = evol_mean |> evo -> query_weight(evo, :, val_t, freq_query)
    spct_mean_mask = evol_mean |> evo -> query_weight(evo, mask_evol, val_t, freq_query)
    spectra_reps_mask = [
        evols[r] |> ev -> query_weight(ev, mask_evol, val_t, freq_query)
        for r in 1:n_rep
    ]
    return (; n_rep, val_t, val_t_sel, freq_query, mask_evol, evols, evol_mean, spectra_reps_mask, spct_mean_full, spct_mean_mask)
end

function query_weight(evol, mask, t_vec, freq_query; scaling::Real=1000.0, weight=nothing)
    mask_use = mask isa Colon ? trues(length(t_vec)) : mask
    evol_sel = evol[mask_use]
    t_sel = t_vec[mask_use]
    if isnothing(weight)
        weight_sel = ones(Float64, length(evol_sel))
        evol_centered = evol_sel .- mean(evol_sel)
    else
        weight_sel = Float64.(weight[mask_use])
        sum_weight = sum(weight_sel)
        mean_weighted = sum_weight > 0 ? sum(evol_sel .* weight_sel) / sum_weight : mean(evol_sel)
        evol_centered = evol_sel .- mean_weighted
    end
    spct = [
        sum(@. weight_sel * evol_centered * exp(-2im * pi * freq_query[f] * t_sel / scaling))
        for f in freq_query] |> e -> abs2.(e)
    spct_max = maximum(spct)
    return isfinite(spct_max) && spct_max > 0 ? spct ./ spct_max : zeros(eltype(spct), size(spct))
end

function get_hint_value(hints::NamedTuple, name::Symbol, fallback)
    return haskey(hints, name) ? getproperty(hints, name) : fallback
end

function fit_evol_oscillation_decay(
    t_vec::AbstractVector{<:Real},
    evol::AbstractVector{<:Real};
    weight::Union{Nothing,AbstractVector{<:Real}}=nothing,
    A_hint=(max=maximum(evol) - minimum(evol), min=0.0, init=(maximum(evol) - minimum(evol)) / 2),
    C_hint=(max=maximum(evol), min=minimum(evol), init=mean(evol)),
    λ_hint=(max=1000, min=20, init=50),
    ν_hint=(max=85, min=10, init=50),
    φ_hint=(max=π, min=-π, init=0),
)
    length(t_vec) == length(evol) ||
        throw(DimensionMismatch("t_vec length $(length(t_vec)) must match evol length $(length(evol))"))
    if !isnothing(weight)
        length(weight) == length(evol) ||
            throw(DimensionMismatch("weight length $(length(weight)) must match evol length $(length(evol))"))
    end
    p_init = Float64[A_hint.init, C_hint.init, λ_hint.init, ν_hint.init, φ_hint.init]
    p_lower = Float64[A_hint.min, C_hint.min, λ_hint.min, ν_hint.min, φ_hint.min]
    p_upper = Float64[A_hint.max, C_hint.max, λ_hint.max, ν_hint.max, φ_hint.max]
    if isnothing(weight)
        fit = curve_fit(fit_evol_oscillation_decay_model, t_vec, evol, p_init; lower=p_lower, upper=p_upper)
        weight_fit = ones(Float64, length(evol))
    else
        weight_fit = Float64.(weight)
        any(>(0), weight_fit) || throw(ArgumentError("fit weights must contain at least one positive value."))
        fit = curve_fit(fit_evol_oscillation_decay_model, t_vec, evol, weight_fit, p_init; lower=p_lower, upper=p_upper)
    end
    params_fit = coef(fit)
    evol_fit = fit_evol_oscillation_decay_model(t_vec, params_fit)
    rss = sqrt(sum(@. weight_fit * (evol - evol_fit)^2))
    norm_evol = sqrt(sum(@. weight_fit * evol^2))
    rel_residue = isfinite(norm_evol) && norm_evol > 0 ? rss / norm_evol : Inf
    return (;
        model=:oscillation_decay,
        model_function="fit_evol_oscillation_decay_model",
        params=params_fit,
        A=params_fit[1],
        C=params_fit[2],
        λ=params_fit[3],
        ν=params_fit[4],
        φ=params_fit[5],
        rel_residue,
    )
end

function fit_evol_property_variant(
    trend_reps::AbstractVector{<:AbstractDict},
    spec,
    variant;
    fit_evol,
)
    isnothing(fit_evol) && return nothing
    fit_evol.model == :oscillation_decay ||
        throw(ArgumentError("unsupported evol fit model $(fit_evol.model) for property $(spec.name)"))
    evol_kind, _ = trend_variant_evol_spct(variant)
    key_evol = "evol-$evol_kind-$(variant.name)"
    fidl_key = trend_variant_fidl_key(spec.name, variant.name)
    hints_user = hasproperty(fit_evol, :kwargs) ? fit_evol.kwargs : NamedTuple()

    t_all = Float64[]
    evol_all = Float64[]
    weight_all = Float64[]
    for trend in trend_reps
        for key in (key_evol, spec.selection_key, "t_vec")
            haskey(trend, key) || throw(KeyError(key))
        end
        mask_sel = in.(trend["t_vec"], Ref(trend[spec.selection_key]))
        append!(t_all, Float64.(trend["t_vec"][mask_sel]))
        append!(evol_all, Float64.(trend[key_evol][mask_sel]))
        if isnothing(fidl_key)
            append!(weight_all, ones(Float64, count(mask_sel)))
        else
            haskey(trend, fidl_key) || throw(KeyError(fidl_key))
            append!(weight_all, Float64.(trend[fidl_key][mask_sel]))
        end
    end
    !isempty(t_all) || throw(ArgumentError("no selected evol data for $(variant.name)"))
    any(>(0), weight_all) || return nothing

    hint_default_A = (max=maximum(evol_all) - minimum(evol_all), min=0.0, init=(maximum(evol_all) - minimum(evol_all)) / 2)
    hint_default_C = (max=maximum(evol_all), min=minimum(evol_all), init=mean(evol_all))
    fit = fit_evol_oscillation_decay(
        t_all,
        evol_all;
        weight=weight_all,
        A_hint=get_hint_value(hints_user, :A_hint, hint_default_A),
        C_hint=get_hint_value(hints_user, :C_hint, hint_default_C),
        λ_hint=get_hint_value(hints_user, :λ_hint, (max=1000, min=20, init=50)),
        ν_hint=get_hint_value(hints_user, :ν_hint, (max=85, min=10, init=30)),
        φ_hint=get_hint_value(hints_user, :φ_hint, (max=2π, min=-2π, init=-0)),
    )
    return merge(
        fit,
        (;
            property=spec.name,
            variant=variant.name,
            selection_key=spec.selection_key,
            t_fit=(minimum(t_all), maximum(t_all)),
            n_sample=length(t_all),
            n_weight_nonzero=count(>(0), weight_all),
        ),
    )
end

function fit_evol_properties_from_trends(
    trend_sidepeak_nvlp::AbstractArray,
    property_specs::AbstractVector,
)
    return [
        begin
            fits = Dict{String,Any}()
            for spec in property_specs
                (hasproperty(spec, :fit_evol) && !isnothing(spec.fit_evol)) || continue
                for variant in spec.variants
                    fits[variant.name] = fit_evol_property_variant(
                        vec(trend_sidepeak_nvlp[c, :, i]),
                        spec,
                        variant;
                        fit_evol=spec.fit_evol,
                    )
                end
            end
            fits
        end
        for c in axes(trend_sidepeak_nvlp, 1), i in axes(trend_sidepeak_nvlp, 3)
    ]
end

function default_selector_t_spectrum(;
    selector_t_sidepeak::Function=t -> trues(length(t)),
    selector_t_envelope::Function=t -> trues(length(t)),
)
    return (;
        number=selector_t_sidepeak,
        sp_weight=selector_t_sidepeak,
        sp_height=selector_t_sidepeak,
        sp_width=selector_t_sidepeak,
        sp_wavenum=selector_t_sidepeak,
        nvlp=selector_t_envelope,
    )
end

function normalize_selector_t_spectrum(
    selector_t_spectrum::NamedTuple;
    selector_t_sidepeak::Union{Nothing,Function}=nothing,
    selector_t_envelope::Union{Nothing,Function}=nothing,
)
    selector_default = default_selector_t_spectrum(
        ;
        selector_t_sidepeak=isnothing(selector_t_sidepeak) ? (t -> trues(length(t))) : selector_t_sidepeak,
        selector_t_envelope=isnothing(selector_t_envelope) ? (t -> trues(length(t))) : selector_t_envelope,
    )
    return merge(selector_default, selector_t_spectrum)
end

function calc_selector_mask(selector::Function, t_vec::AbstractVector, name_selector)
    mask = selector(t_vec)
    mask isa AbstractVector{Bool} ||
        throw(ArgumentError("selector_t_spectrum.$name_selector must return a boolean vector for t_vec input."))
    length(mask) == length(t_vec) ||
        throw(DimensionMismatch("selector_t_spectrum.$name_selector returned length $(length(mask)), expected $(length(t_vec))."))
    any(mask) ||
        throw(ArgumentError("selector_t_spectrum.$name_selector selected no time points."))
    return mask
end

function anlz_trend_from_extr(
    t_vec::AbstractVector{<:Real},
    extr::AbstractVector{SoloExtract},
    freq_query::AbstractVector{<:Real};
    selector_t_spectrum::NamedTuple=NamedTuple(),
    selector_t_sidepeak::Union{Nothing,Function}=nothing,
    selector_t_envelope::Union{Nothing,Function}=nothing,
    query_weight_kwargs::NamedTuple=NamedTuple(),
)
    selector_t_spectrum = normalize_selector_t_spectrum(
        selector_t_spectrum;
        selector_t_sidepeak,
        selector_t_envelope,
    )
    mask_sel = Dict(
        key => calc_selector_mask(selector_t_spectrum[key], t_vec, key)
        for key in propertynames(selector_t_spectrum)
    )
    t_vec_sel = Dict(key => t_vec[mask] for (key, mask) in mask_sel)

    evol_fit_sp_fidl = extr |> e -> map(t -> t.sidepeak.params_tailess.rel_residue > 0.6 ? 0 : 1, e)
    evol_moment_sp_fidl = extr |> e -> map(t -> t.sidepeak.moments.weight < 0.05 ? 0 : 1, e)
    query_weight_sel = key -> evol -> query_weight(evol, mask_sel[key], t_vec, freq_query; query_weight_kwargs...)
    query_weight_sel_fidl = (key, fidl) -> evol -> query_weight(evol, mask_sel[key], t_vec, freq_query; weight=fidl, query_weight_kwargs...)
    evol_dens_sum = extr |> e -> map(t -> t.essentials.sum_dens_full, e)
    evol_fit_weight = extr |> e -> map(t -> t.sidepeak.params_tailess.weight, e)
    evol_fit_height = extr |> e -> map(t -> t.sidepeak.params_tailess.height, e)
    evol_fit_wavenum = extr |> e -> map(t -> t.sidepeak.params_tailess.wavenum, e)
    evol_fit_width = extr |> e -> map(t -> t.sidepeak.params_tailess.width, e)
    evol_fit_size_x = extr |> e -> map(t -> t.envelope.params_asymm.size[1], e)
    evol_fit_size_y = extr |> e -> map(t -> t.envelope.params_asymm.size[2], e)
    evol_fit_cent_x = extr |> e -> map(t -> t.envelope.params_asymm.cent[1], e)
    evol_fit_cent_y = extr |> e -> map(t -> t.envelope.params_asymm.cent[2], e)
    evol_moment_weight = extr |> e -> map(t -> t.sidepeak.moments.weight, e)
    evol_moment_height = extr |> e -> map(t -> t.sidepeak.moments.height, e)
    evol_moment_wavenum = extr |> e -> map(t -> t.sidepeak.moments.wavenum, e)
    evol_moment_width = extr |> e -> map(t -> t.sidepeak.moments.width, e)

    spct_dens_sum = evol_dens_sum |> query_weight_sel(:number)
    spct_fit_weight = evol_fit_weight |> query_weight_sel(:sp_weight)
    spct_fit_height = evol_fit_height |> query_weight_sel(:sp_height)
    spct_fit_wavenum = evol_fit_wavenum |> query_weight_sel_fidl(:sp_wavenum, evol_fit_sp_fidl)
    spct_fit_width = evol_fit_width |> query_weight_sel_fidl(:sp_width, evol_fit_sp_fidl)
    spct_moment_weight = evol_moment_weight |> query_weight_sel(:sp_weight)
    spct_moment_height = evol_moment_height |> query_weight_sel(:sp_height)
    spct_moment_wavenum = evol_moment_wavenum |> query_weight_sel_fidl(:sp_wavenum, evol_moment_sp_fidl)
    spct_moment_width = evol_moment_width |> query_weight_sel_fidl(:sp_width, evol_moment_sp_fidl)
    spct_fit_size_x = evol_fit_size_x |> query_weight_sel(:nvlp)
    spct_fit_size_y = evol_fit_size_y |> query_weight_sel(:nvlp)
    spct_fit_cent_x = evol_fit_cent_x |> query_weight_sel(:nvlp)
    spct_fit_cent_y = evol_fit_cent_y |> query_weight_sel(:nvlp)
    return Dict(
        "t_vec" => t_vec,
        "t_vec_sel_number" => t_vec_sel[:number],
        "t_vec_sel_sp_weight" => t_vec_sel[:sp_weight],
        "t_vec_sel_sp_height" => t_vec_sel[:sp_height],
        "t_vec_sel_sp_width" => t_vec_sel[:sp_width],
        "t_vec_sel_sp_wavenum" => t_vec_sel[:sp_wavenum],
        "t_vec_sel_nvlp_size" => t_vec_sel[:nvlp],
        "t_vec_sel_nvlp_cent" => t_vec_sel[:nvlp],
        "t_vec_sel_sp" => t_vec_sel[:sp_weight],
        "t_vec_sel_nvlp" => t_vec_sel[:nvlp],
        "mask_sel" => mask_sel[:sp_weight],
        "freq_query" => freq_query,
        "evol-all-dens-sum" => evol_dens_sum,
        "evol-all-fit-sp-fidl" => evol_fit_sp_fidl,
        "evol-all-moment-sp-fidl" => evol_moment_sp_fidl,
        "evol-all-fit-weight" => evol_fit_weight,
        "evol-all-fit-height" => evol_fit_height,
        "evol-all-fit-wavenum" => evol_fit_wavenum,
        "evol-all-fit-width" => evol_fit_width,
        "evol-all-moment-weight" => evol_moment_weight,
        "evol-all-moment-height" => evol_moment_height,
        "evol-all-moment-wavenum" => evol_moment_wavenum,
        "evol-all-moment-width" => evol_moment_width,
        "evol-all-fit-size-x" => evol_fit_size_x,
        "evol-all-fit-size-y" => evol_fit_size_y,
        "evol-all-fit-cent-x" => evol_fit_cent_x,
        "evol-all-fit-cent-y" => evol_fit_cent_y,
        "spct-sel-dens-sum" => spct_dens_sum,
        "spct-sel-fit-weight" => spct_fit_weight,
        "spct-sel-fit-height" => spct_fit_height,
        "spct-sel-fit-wavenum" => spct_fit_wavenum,
        "spct-sel-fit-width" => spct_fit_width,
        "spct-sel-moment-weight" => spct_moment_weight,
        "spct-sel-moment-height" => spct_moment_height,
        "spct-sel-moment-wavenum" => spct_moment_wavenum,
        "spct-sel-moment-width" => spct_moment_width,
        "spct-sel-fit-size-x" => spct_fit_size_x,
        "spct-sel-fit-size-y" => spct_fit_size_y,
        "spct-sel-fit-cent-x" => spct_fit_cent_x,
        "spct-sel-fit-cent-y" => spct_fit_cent_y,
    )
end

function matches_pattern(pattern, key)
    key_str = string(key)
    if pattern isa AbstractVector || pattern isa Tuple
        return any(p -> matches_pattern(p, key_str), pattern)
    end
    return occursin(Regex(pattern), key_str)
end

function mean_dict(dicts::AbstractArray{<:AbstractDict}; pattern_incl=".*", pattern_excl="(?!)")
    !isempty(dicts) || throw(ArgumentError("dicts must contain at least one dictionary."))

    keys_sel = dicts |>
               ds -> map(keys, ds) |>
                     Iterators.flatten |>
                     unique |>
                     ks -> filter(k -> matches_pattern(pattern_incl, k) && !matches_pattern(pattern_excl, k), ks)

    return Dict(
        key => begin
            all(d -> haskey(d, key), dicts) || throw(ArgumentError("selected key $key is not present in every dictionary."))
            mean([d[key] for d in dicts])
        end
        for key in keys_sel
    )
end
