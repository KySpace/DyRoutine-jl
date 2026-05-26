using Statistics

function calc_mean_std(q::AbstractArray{<:Real})
    mean_q = Statistics.mean(q)
    std_q = Statistics.std(q)
    return mean_q, std_q
end

function calc_stacked_essn(essns::AbstractVector{SoloEssentials})::SoloEssentials
    n_essn = length(essns)
    n_essn > 0 || throw(ArgumentError("essns must contain at least one SoloEssentials."))
    essn_ref = first(essns)
    mean_tuple(ts) = begin
        t_ref = first(ts)
        ntuple(idx -> mean(t -> t[idx], ts), length(t_ref))
    end
    return SoloEssentials(
        mean(map(essn -> essn.dens2d, essns)),
        mean(map(essn -> essn.modl2d, essns)),
        mean(map(essn -> essn.dens2d_core, essns)),
        mean_tuple(map(essn -> essn.offset_cent_core, essns)),
        mean_tuple(map(essn -> essn.smwh_core, essns)),
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
