include(joinpath(@__DIR__, "dualmotcommons.jl"))

function max_number_with_error(nums::AbstractArray, stds::AbstractArray)
    size(nums) == size(stds) || throw(DimensionMismatch("number and deviation arrays must match"))
    values = vec(nums .+ ifelse.(isfinite.(stds), max.(stds, 0), 0.0))
    finite_values = filter(isfinite, values)
    isempty(finite_values) && throw(ArgumentError("no finite loading numbers available"))
    maximum(finite_values)
end


function loading_421_extrema(stats::NamedTuple, γ_active::Real)
    idx_axis = Dict(key => idx for (idx, key) in enumerate(stats.name_stat))
    idx_bias = idx_axis[:β_MOT]
    idx_loadcfg = idx_axis[:loadcfg]
    maxima_num = Dict(:t => Float64[], :n => Float64[])
    maxima_time = Dict(:t => Float64[], :n => Float64[])
    for (ibias, bias) in enumerate(stats.vars.β_MOT),
        (iloadcfg, loadcfg) in enumerate(stats.vars.loadcfg)
        # Preserve the n-balanced plotting rule: DCS is drawn but does not set limits.
        iszero(bias) || loadcfg != :DCS || continue
        balance = Symbol(balance_tag(bias))
        indices = Any[Colon() for _ in stats.name_stat]
        indices[idx_bias] = ibias
        indices[idx_loadcfg] = iloadcfg
        push!(maxima_num[balance], max_number_with_error(
            @view(stats.num_stat[indices...]), @view(stats.std_num_stat[indices...])))
        scale_time = iszero(bias) ? γ_active : 1.0
        push!(maxima_time[balance], maximum(stats.vars.t_load) * scale_time)
    end
    extrema = Dict{Symbol, Tuple{Float64, Float64}}()
    for balance in (:t, :n)
        isempty(maxima_num[balance]) ||
            (extrema[balance] = (maximum(maxima_time[balance]), maximum(maxima_num[balance])))
    end
    isempty(extrema) && throw(ArgumentError("no 421 loading curves eligible to set limits"))
    extrema
end


function loading_runinfo_stats(runinfo::NamedTuple;
    label::AbstractString,
    bounds_sigmax_num::Tuple{<:Real,<:Real},
    bounds_sigmay_num::Tuple{<:Real,<:Real},
    num_max_num::Real,
)
    blocks = [calc_num_evol_block(data;
        label="$label data[$idx]", bounds_sigmax_num, bounds_sigmay_num,
        num_max_num) for (idx, data) in enumerate(runinfo.data)]
    combined = combine_num_evol_blocks(blocks; label)
    idx_rep_axis = findfirst(==(:rep), keys(combined.vars))
    num_stat = dropdims(mapslices(combined.num_fmt; dims=idx_rep_axis) do values
        valid = collect(skipmissing(vec(values)))
        isempty(valid) ? NaN : mean(valid)
    end; dims=idx_rep_axis)
    std_num_stat = dropdims(mapslices(combined.num_fmt; dims=idx_rep_axis) do values
        valid = collect(skipmissing(vec(values)))
        length(valid) < 2 ? NaN : std(valid)
    end; dims=idx_rep_axis)
    name_stat = Tuple(key for key in keys(combined.vars) if key != :rep)
    (; vars=combined.vars, name_stat, num_stat, std_num_stat)
end

function preprocess_loading_limits()
    path_data = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data"
    pairs = ("162-164", "160-162", "161-162", "161-164", "163-164",
        "162-163", "161-163")
    bounds_sigmax_num = (4e-4, Inf)
    bounds_sigmay_num = (4e-4, Inf)
    num_max_num = 2e8
    γ_active = 0.43
    limits = Dict{Tuple{String, Symbol}, NamedTuple}()

    for pair in pairs
        runinfos_421 = read_num_evol_runinfos(joinpath(path_data, "MOT loading 421"), pair;
            var_specs=DUALMOT_LOADING_VAR_SPECS,
            validate_vars=validate_dualmot_vars)
        extrema_421 = [loading_421_extrema(loading_runinfo_stats(runinfo;
                label="preprocess 421 $pair/$(runinfo.tag)", bounds_sigmax_num,
                bounds_sigmay_num, num_max_num), γ_active)
            for runinfo in runinfos_421]
        runinfos_626 = read_num_evol_runinfos(
            joinpath(path_data, "MOT loading 626"), pair;
            var_specs=DUALMOT_LOADCFG_VAR_SPECS,
            validate_vars=validate_loadcfg_comparison_vars)
        stats_626 = [loading_runinfo_stats(runinfo;
            label="preprocess 626 $pair/$(runinfo.tag)", bounds_sigmax_num,
            bounds_sigmay_num, num_max_num) for runinfo in runinfos_626]
        max_num_626 = maximum(max_number_with_error(stats.num_stat,
            stats.std_num_stat) for stats in stats_626)
        max_time_626 = maximum(maximum(stats.vars.t_load) for stats in stats_626) * γ_active

        extrema_t = [(max_time_626, max_num_626)]
        append!(extrema_t, [extrema[:t] for extrema in extrema_421 if haskey(extrema, :t)])
        xmax_t = maximum(first, extrema_t)
        ymax_t = maximum(last, extrema_t) / 1e7
        limits[(pair, :t)] = (x=(0.0, 1.03 * xmax_t), y=(0.0, 1.05 * ymax_t))
        println("$pair shared t-balanced linear loading limits: $(limits[(pair, :t)])")

        extrema_n = [extrema[:n] for extrema in extrema_421 if haskey(extrema, :n)]
        if !isempty(extrema_n)
            xmax_n = maximum(first, extrema_n)
            ymax_n = maximum(last, extrema_n) / 1e7
            limits[(pair, :n)] = (x=(0.0, 1.03 * xmax_n), y=(0.0, 1.05 * ymax_n))
            println("$pair shared n-balanced linear loading limits: $(limits[(pair, :n)])")
        end
    end
    limits
end

limits_loading_linear = preprocess_loading_limits()
