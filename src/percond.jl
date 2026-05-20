using Statistics

function reshape_variation(dens::AbstractArray{<:Real,3}, val)
    variation, height, width = size(dens)
    variation_expected = length(val[1]) * length(val[2]) * length(val[3])
    variation == variation_expected || throw(
        ArgumentError("Expected variation axis to have length $variation_expected, got $variation."),
    )

    return dens |>
           x -> reshape(x, length(val[3]), length(val[2]), length(val[1]), height, width) |>
                x -> permutedims(x, (3, 2, 1, 4, 5))
end

function calc_mean_std(q::AbstractArray{<:Real})
    mean_q = Statistics.mean(q)
    std_q = Statistics.std(q)
    return mean_q, std_q
end

function calc_stacked_essn(essns::AbstractVector{SoloEssentials})::SoloEssentials
    n_essn = length(essns)
    n_essn > 0 || throw(ArgumentError("essns must contain at least one SoloEssentials."))
    essn_ref = first(essns)
    return SoloEssentials(
        mean(map(essn -> essn.dens2d, essns)),
        mean(map(essn -> essn.modl2d, essns)),
        mean(map(essn -> essn.prfl_strip, essns)),
        mean(map(essn -> essn.prfl_modl, essns)),
        mean(map(essn -> essn.prfl_modl_norm_px, essns)),
        essn_ref.smwh,
        essn_ref.smwh_strip,
        essn_ref.smw_modl,
        essn_ref.step_posi,
        essn_ref.step_modl,
        mean(map(essn -> essn.sum_dens_full, essns)),
    )
end
