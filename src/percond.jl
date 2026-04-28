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
    std_q =  Statistics.std(q)
    return mean_q, std_q
end
