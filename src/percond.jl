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
        mean_tuple(map(essn -> essn.smwh_core, essns)) |> wh -> map(Int, wh),
        mean(map(essn -> essn.prfl_strip, essns)),
        (;
            main=(;
                norm=mean(map(essn -> essn.prfl_modl.main.norm, essns)),
                raw=mean(map(essn -> essn.prfl_modl.main.raw, essns)),
                normed_px=mean(map(essn -> essn.prfl_modl.main.normed_px, essns)),
            ),
            side=(;
                raw=mean(map(essn -> essn.prfl_modl.side.raw, essns)),
                normed_px=mean(map(essn -> essn.prfl_modl.side.normed_px, essns)),
            ),
            mask=essn_ref.prfl_modl.mask,
        ),
        essn_ref.smwh,
        essn_ref.smwh_strip,
        essn_ref.step_posi,
        essn_ref.step_modl,
        mean(map(essn -> essn.sum_dens_full, essns)),
    )
end
