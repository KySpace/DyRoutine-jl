using Printf

function gen_date_path(mmdd::AbstractString, year::Integer)
    @assert length(mmdd) == 4 "mmdd should be a 4-character string like \"0325\""
    month = mmdd[1:2]
    return "$(year)-$(month)\\$(mmdd)"
end

function gen_h5name(mmdd::AbstractString, run::Integer)
    @assert length(mmdd) == 4 "mmdd should be a 4-character string like \"0325\""
    month = mmdd[1:2]
    return "d$(mmdd)r$(@sprintf("%02d", run)).h5"
end
"""
    flatten_with_rebuilder(x)

Flatten either

1. an N-dimensional array of real numbers, e.g. `Array{Float64, N}`
2. an N-dimensional array whose elements are themselves arrays of real numbers

Return `(v, rebuild)`, where `v` is a flat vector and `rebuild(v2)` restores
the original structure.
"""
function flatten_with_rebuilder(x::AbstractArray{<:Real})
    sz = size(x)
    T = eltype(x)
    v = collect(vec(x))
    function rebuild(v2::AbstractVector)
        length(v2) == length(v) || throw(DimensionMismatch("wrong vector length"))
        return reshape(collect(T, v2), sz)
    end
    return v, rebuild
end

function flatten_with_rebuilder(x::AbstractArray{<:AbstractArray{<:Real}})
    outer_sz = size(x)
    inner_sizes = map(size, x)
    inner_lengths = map(length, x)
    inner_types = map(eltype, x)
    offsets = cumsum([0; inner_lengths...])
    total_len = offsets[end]
    T = promote_type(inner_types...)
    v = Vector{T}(undef, total_len)
    for (i, a) in enumerate(x)
        r = offsets[i] + 1 : offsets[i + 1]
        v[r] .= vec(a)
    end
    function rebuild(v2::AbstractVector)
        length(v2) == total_len || throw(DimensionMismatch("wrong vector length"))
        ys = Vector{Any}(undef, length(x))
        for i in eachindex(x)
            r = offsets[i] + 1 : offsets[i + 1]
            Ti = inner_types[i]
            ys[i] = reshape(collect(Ti, v2[r]), inner_sizes[i])
        end
        return reshape(ys, outer_sz)
    end
    return v, rebuild
end
