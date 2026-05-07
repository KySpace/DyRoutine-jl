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

function plot_mode_evol_freq(axs::Dict{String,Axis}, mode::ModeWeight, val_t::AbstractVector)
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
