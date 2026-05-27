using MultivariateStats: PCA, fit, predict, projection

struct ModeWeight{TProfile<:AbstractArray,TWeight<:AbstractArray}
    profile::TProfile
    weight::TWeight
end

function build_pca_matrix(samples::AbstractArray{<:AbstractArray{<:Real}})

    n_sample = length(samples)
    n_sample > 0 || throw(ArgumentError("samples must contain at least one sample array."))

    sample_first = first(samples)
    sz_sample = size(sample_first)
    n_feature = length(sample_first)
    n_feature > 0 || throw(ArgumentError("sample arrays must contain at least one value."))

    mat_sample = Matrix{Float64}(undef, n_feature, n_sample)
    for (idx_sample, sample) in pairs(samples)
        size(sample) == sz_sample || throw(DimensionMismatch("sample at index $idx_sample has size $(size(sample)); expected $sz_sample."))
        mat_sample[:, LinearIndices(samples)[idx_sample]] .= vec(sample)
    end

    return mat_sample, sz_sample
end

function fit_pca_modes(n_mode::Integer, samples::AbstractArray{<:AbstractArray{<:Real}})
    n_mode > 0 || throw(ArgumentError("n_mode=$n_mode must be positive."))
    n_mode_int = Int(n_mode)

    mat_sample, sample_shape = build_pca_matrix(samples)
    n_feature, n_sample = size(mat_sample)
    n_sample > 1 || throw(ArgumentError("PCA requires at least two samples; got n_sample=$n_sample."))
    n_mode_int <= min(n_feature, n_sample - 1) || throw(ArgumentError("n_mode=$n_mode exceeds min(n_feature=$n_feature, n_sample - 1=$(n_sample - 1))."))

    pca_fit = fit(PCA, mat_sample; maxoutdim=n_mode_int, pratio=1.0)
    mat_profile = projection(pca_fit)
    mat_weight = predict(pca_fit, mat_sample)
    size(mat_profile, 2) >= n_mode_int || throw(ArgumentError("PCA returned only $(size(mat_profile, 2)) modes; requested n_mode=$n_mode. Check for low-rank or constant sample data."))

    return [
        ModeWeight(
            reshape(copy(@view mat_profile[:, idx_mode]), sample_shape),
            reshape(copy(@view mat_weight[idx_mode, :]), size(samples)),
        )
        for idx_mode in 1:n_mode_int
    ]
end

function query_weight(evo, mask, t_vec, freq_query; scaling::Real=1000.0)
    weight = evo[mask] |> e -> e .- mean(e) |> e -> [
        sum(@. e * exp(-2im * pi * freq_query[f] * t_vec[mask] / scaling))
        for f in freq_query] |> e -> abs.(e) .^ 2
    return weight / sum(weight)
end

function anlz_trend_from_extr(
    t_vec::AbstractVector{<:Real},
    extr::AbstractVector{SoloExtract},
    freq_query::AbstractVector{<:Real};
    selector_t_sidepeak::Function,
    selector_t_envelope::Function,
    query_weight_kwargs::NamedTuple=NamedTuple(),
)
    mask_sel_sp = selector_t_sidepeak(t_vec)
    t_vec_sel_sp = t_vec[mask_sel_sp]
    mask_sel_nvlp = selector_t_envelope(t_vec)
    t_vec_sel_nvlp = t_vec[mask_sel_nvlp]

    query_weight_sel_sp = evo -> query_weight(evo, mask_sel_sp, t_vec, freq_query; query_weight_kwargs...)
    query_weight_sel_nvlp = evo -> query_weight(evo, mask_sel_nvlp, t_vec, freq_query; query_weight_kwargs...)
    evo_dens_sum = extr |> e -> map(t -> t.essentials.sum_dens_full, e)
    evo_fit_weight = extr |> e -> map(t -> t.sidepeak.params_tailess.weight, e)
    evo_fit_height = extr |> e -> map(t -> t.sidepeak.params_tailess.height, e)
    evo_fit_wavenum = extr |> e -> map(t -> t.sidepeak.params_tailess.wavenum, e)
    evo_fit_width = extr |> e -> map(t -> t.sidepeak.params_tailess.width, e)
    evo_fit_size_x = extr |> e -> map(t -> t.envelope.params_asymm.size[1], e)
    evo_fit_size_y = extr |> e -> map(t -> t.envelope.params_asymm.size[2], e)
    evo_moment_weight = extr |> e -> map(t -> t.sidepeak.moments.weight, e)
    evo_moment_height = extr |> e -> map(t -> t.sidepeak.moments.height, e)
    evo_moment_wavenum = extr |> e -> map(t -> t.sidepeak.moments.wavenum, e)
    evo_moment_width = extr |> e -> map(t -> t.sidepeak.moments.width, e)

    ft_fit_weight = evo_fit_weight |> query_weight_sel_sp
    ft_fit_height = evo_fit_height |> query_weight_sel_sp
    ft_fit_wavenum = evo_fit_wavenum |> query_weight_sel_sp
    ft_fit_width = evo_fit_width |> query_weight_sel_sp
    ft_moment_weight = evo_moment_weight |> query_weight_sel_sp
    ft_moment_height = evo_moment_height |> query_weight_sel_sp
    ft_moment_wavenum = evo_moment_wavenum |> query_weight_sel_sp
    ft_moment_width = evo_moment_width |> query_weight_sel_sp
    ft_fit_size_x = evo_fit_size_x |> query_weight_sel_nvlp
    ft_fit_size_y = evo_fit_size_y |> query_weight_sel_nvlp
    return Dict(
        "t_vec" => t_vec,
        "t_vec_sel_sp" => t_vec_sel_sp,
        "t_vec_sel_nvlp" => t_vec_sel_nvlp,
        "mask_sel" => mask_sel_sp,
        "freq_query" => freq_query,
        "evol-all-dens-sum" => evo_dens_sum,
        "evol-all-fit-weight" => evo_fit_weight,
        "evol-all-fit-height" => evo_fit_height,
        "evol-all-fit-wavenum" => evo_fit_wavenum,
        "evol-all-fit-width" => evo_fit_width,
        "evol-all-moment-weight" => evo_moment_weight,
        "evol-all-moment-height" => evo_moment_height,
        "evol-all-moment-wavenum" => evo_moment_wavenum,
        "evol-all-moment-width" => evo_moment_width,
        "evol-all-fit-size-x" => evo_fit_size_x,
        "evol-all-fit-size-y" => evo_fit_size_y,
        "freq-sel-fit-weight" => ft_fit_weight,
        "freq-sel-fit-height" => ft_fit_height,
        "freq-sel-fit-wavenum" => ft_fit_wavenum,
        "freq-sel-fit-width" => ft_fit_width,
        "freq-sel-moment-weight" => ft_moment_weight,
        "freq-sel-moment-height" => ft_moment_height,
        "freq-sel-moment-wavenum" => ft_moment_wavenum,
        "freq-sel-moment-width" => ft_moment_width,
        "freq-sel-fit-size-x" => ft_fit_size_x,
        "freq-sel-fit-size-y" => ft_fit_size_y,
    )
end
