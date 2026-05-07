using Printf

function gen_date_path(mmdd::AbstractString, year::Integer)
    @assert length(mmdd) == 4 "mmdd should be a 4-character string like \"0325\""
    month = mmdd[1:2]
    return "$(year)-$(month)/$(mmdd)"
end

function gen_h5name(mmdd::AbstractString, run::Integer)
    @assert length(mmdd) == 4 "mmdd should be a 4-character string like \"0325\""
    month = mmdd[1:2]
    return "d$(mmdd)r$(@sprintf("%02d", run)).h5"
end
