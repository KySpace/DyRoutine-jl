using MultivariateStats: PCA, fit, predict, projection
using CairoMakie: Axis
using Colors: Oklch

struct ModeWeight{TProfile<:AbstractArray,TWeight<:AbstractArray}
    profile::TProfile
    weight::TWeight
end

function gen_clrmap_posneg(hue_pos, hue_neg)
    return [Oklch(1 - abs(t), 0.24 * abs(t), t > 0 ? hue_pos : hue_neg) |> c -> RGBAf(c) for t in range(-1, 1; length=256)]
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

function plot_mode_evol_freq_duet!(axs::Dict{String,Axis}, mode::ModeWeight, val_t::AbstractVector)
    ndims(mode.profile) == 3 && size(mode.profile, 1) == 2 || throw(ArgumentError("mode.profile must be a 3D array with size[1]==2."))
    clrmap = gen_clrmap_posneg(0.60 * 360, 0.96 * 360)
    c = maximum(abs, mode.profile)
    heatmap!(axs["l"], mode.profile[1, :, :]'; colormap=clrmap, colorrange=(-c, c))
    heatmap!(axs["r"], mode.profile[2, :, :]'; colormap=clrmap, colorrange=(-c, c))
    axs["l"].aspect = DataAspect()
    axs["r"].aspect = DataAspect()
    axs["l"] |> hidedecorations!
    axs["r"] |> hidedecorations!
    for rep = 1:size(mode.weight, 1)
        lines!(axs["evol"], val_t, mode.weight[rep, :]; color=(:black, 0.2))
    end
end

function plot_mode_evol_freq_solo!(axs::Dict{String,Axis}, mode::ModeWeight, val_t::AbstractVector)
    ndims(mode.profile) == 2 || throw(ArgumentError("mode.profile must be a 2D array. "))
    clrmap = gen_clrmap_posneg(0.60 * 360, 0.96 * 360)
    c = maximum(abs, mode.profile)
    heatmap!(axs["mode"], mode.profile[:, :]; colormap=clrmap, colorrange=(-c, c))
    axs["mode"].aspect = DataAspect()
    axs["mode"] |> hidedecorations!
    for rep = 1:size(mode.weight, 1)
        lines!(axs["evol"], val_t, mode.weight[rep, :]; color=(:black, 0.2))
    end
end

function anlz_trend_from_extr(t_vec::AbstractVector{<:Real}, extr::AbstractVector{SoloExtract}, selector_t::Function, freq_query::AbstractVector{<:Real})
    mask_sel = selector_t(t_vec)
    t_vec_sel = t_vec[mask_sel]
    query_weight = (evo, mask) -> evo[mask_sel] |> e -> e .- mean(e) |> e -> [
        sum(@. e * exp(-2im * pi * freq_query[f] * t_vec[mask_sel] / 1000.0))
        for f in freq_query] |> e -> abs.(e) .^ 2
    query_weight_sel = evo -> query_weight(evo, mask_sel)
    evo_fit_weight = extr |> e -> map(t -> t.sidepeak["weight"], e)
    evo_fit_height = extr |> e -> map(t -> t.sidepeak["height"], e)
    evo_fit_wavenum = extr |> e -> map(t -> t.sidepeak["wavenum"], e)
    evo_fit_width = extr |> e -> map(t -> t.sidepeak["width"], e)
    evo_moment_weight = extr |> e -> map(t -> t.moments["weight"], e)
    evo_moment_height = extr |> e -> map(t -> t.moments["height"], e)
    evo_moment_wavenum = extr |> e -> map(t -> t.moments["wavenum"], e)
    evo_moment_width = extr |> e -> map(t -> t.moments["width"], e)

    ft_fit_weight = evo_fit_weight |> query_weight_sel
    ft_fit_height = evo_fit_height |> query_weight_sel
    ft_fit_wavenum = evo_fit_wavenum |> query_weight_sel
    ft_fit_width = evo_fit_width |> query_weight_sel
    ft_moment_weight = evo_moment_weight |> query_weight_sel
    ft_moment_height = evo_moment_height |> query_weight_sel
    ft_moment_wavenum = evo_moment_wavenum |> query_weight_sel
    ft_moment_width = evo_moment_width |> query_weight_sel
    return Dict(
        "t_vec" => t_vec,
        "t_vec_sel" => t_vec_sel,
        "mask_sel" => mask_sel,
        "freq_query" => freq_query,
        "evol-all-fit-weight" => evo_fit_weight,
        "evol-all-fit-height" => evo_fit_height,
        "evol-all-fit-wavenum" => evo_fit_wavenum,
        "evol-all-fit-width" => evo_fit_width,
        "evol-all-moment-weight" => evo_moment_weight,
        "evol-all-moment-height" => evo_moment_height,
        "evol-all-moment-wavenum" => evo_moment_wavenum,
        "evol-all-moment-width" => evo_moment_width,
        "freq-sel-fit-weight" => ft_fit_weight,
        "freq-sel-fit-height" => ft_fit_height,
        "freq-sel-fit-wavenum" => ft_fit_wavenum,
        "freq-sel-fit-width" => ft_fit_width,
        "freq-sel-moment-weight" => ft_moment_weight,
        "freq-sel-moment-height" => ft_moment_height,
        "freq-sel-moment-wavenum" => ft_moment_wavenum,
        "freq-sel-moment-width" => ft_moment_width,
    )
end

function plot_trend_sidepeak!(axs_trend::Dict, trend_sidepeak::AbstractVector, istp)
    hue = hue_theme_istp[istp]
    clr_fit = (:black, 1.0)
    clr_mmt = (:seagreen, 1.0)
    for r = axes(trend_sidepeak, 1)
        trend = trend_sidepeak[r]
        axs = axs_trend["repeats"][r]
        for (k, obj) in axs
            obj isa Axis && empty!(obj)
        end
        vspan!(axs["evol-weight"], trend["t_vec_sel"][1], trend["t_vec_sel"][end]; color=RGBAf(Oklch(0.98, 0.1, hue), 0.2))
        vspan!(axs["evol-height"], trend["t_vec_sel"][1], trend["t_vec_sel"][end]; color=RGBAf(Oklch(0.98, 0.1, hue), 0.2))
        vspan!(axs["evol-width"], trend["t_vec_sel"][1], trend["t_vec_sel"][end]; color=RGBAf(Oklch(0.98, 0.1, hue), 0.2))
        vspan!(axs["evol-wavenum"], trend["t_vec_sel"][1], trend["t_vec_sel"][end]; color=RGBAf(Oklch(0.98, 0.1, hue), 0.2))
        lines!(axs["evol-weight"], trend["t_vec"], trend["evol-all-fit-weight"]; color=clr_fit)
        lines!(axs["evol-height"], trend["t_vec"], trend["evol-all-fit-height"]; color=clr_fit)
        lines!(axs["evol-width"], trend["t_vec"], trend["evol-all-fit-width"]; color=clr_fit)
        lines!(axs["evol-wavenum"], trend["t_vec"], trend["evol-all-fit-wavenum"]; color=clr_fit)
        lines!(axs["evol-weight"], trend["t_vec"], trend["evol-all-moment-weight"]; color=clr_mmt)
        lines!(axs["evol-height"], trend["t_vec"], trend["evol-all-moment-height"]; color=clr_mmt)
        lines!(axs["evol-width"], trend["t_vec"], trend["evol-all-moment-width"]; color=clr_mmt)
        lines!(axs["evol-wavenum"], trend["t_vec"], trend["evol-all-moment-wavenum"]; color=clr_mmt)
        lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-fit-weight"]; color=clr_fit)
        lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-fit-height"]; color=clr_fit)
        lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-fit-width"]; color=clr_fit)
        lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-fit-wavenum"]; color=clr_fit)
        lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-moment-weight"]; color=clr_mmt)
        lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-moment-height"]; color=clr_mmt)
        lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-moment-width"]; color=clr_mmt)
        lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-moment-wavenum"]; color=clr_mmt)
    end
end
