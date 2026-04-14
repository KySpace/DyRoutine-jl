function reshape_variation(dens::AbstractArray{<:Real, 3}, val)
    variation, height, width = size(dens)
    variation_expected = length(val[1]) * length(val[2]) * length(val[3])
    variation == variation_expected || throw(
        ArgumentError("Expected variation axis to have length $variation_expected, got $variation."),
    )

    return dens |>
           x -> reshape(x, length(val[3]), length(val[2]), length(val[1]), height, width) |>
           x -> permutedims(x, (3, 2, 1, 4, 5))
end

function summarize_repeat_number(dens::AbstractArray{<:Real, 3}, val)
    dens_by_variation = dens |> x -> reshape_variation(x, val)
    number_by_repeat = dens_by_variation |>
                       x -> stack(map(calc_number, eachslice(x; dims=(1, 2, 3))); dims=1) |>
                       x -> reshape(x, size(dens_by_variation, 1), size(dens_by_variation, 2), size(dens_by_variation, 3))

    repeat_count = size(number_by_repeat, 1)
    val_number = number_by_repeat |>
                 x -> sum(x; dims=1) ./ repeat_count |>
                 x -> dropdims(x; dims=1)

    err_number = if repeat_count > 1
        number_by_repeat |>
        x -> x .- reshape(val_number, 1, size(val_number, 1), size(val_number, 2)) |>
        x -> sum(x .^ 2; dims=1) ./ (repeat_count - 1) |>
        x -> sqrt.(x) |>
        x -> dropdims(x; dims=1)
    else
        zero(val_number)
    end

    return (; dens_by_variation, number_by_repeat, val_number, err_number)
end
