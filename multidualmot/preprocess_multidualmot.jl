include(joinpath(@__DIR__, "dualmotcommons.jl"))

function max_number_with_error(nums::AbstractArray, stds::AbstractArray)
    size(nums) == size(stds) || throw(DimensionMismatch("number and deviation arrays must match"))
    values = vec(nums .+ ifelse.(isfinite.(stds), max.(stds, 0), 0.0))
    finite_values = filter(isfinite, values)
    isempty(finite_values) && throw(ArgumentError("no finite loading numbers available"))
    maximum(finite_values)
end

function preprocess_loading_limits()
    path_data = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data"
    pairs = ("162-164", "160-162", "161-162", "161-164", "163-164",
        "162-163", "161-163")
    bounds_sigmax_num = (4e-4, Inf)
    bounds_sigmay_num = (4e-4, Inf)
    num_max_num = 2e8
    γ_active = 0.43
    limits = Dict{String, NamedTuple}()

    for pair in pairs
        runinfos_421 = read_num_evol_runinfos(joinpath(path_data, "MOT loading 421"), pair;
            var_specs=DUALMOT_LOADING_VAR_SPECS,
            validate_vars=validate_dualmot_vars)
        runinfo_421 = only(filter(info -> length(info.data) == 1 &&
            0.0 in only(info.data).vars.β_MOT, runinfos_421))
        stats_421 = calc_num_evol_block(only(runinfo_421.data);
            label="preprocess 421 $pair", bounds_sigmax_num, bounds_sigmay_num,
            num_max_num)
        idx_bias = findfirst(==(:β_MOT), stats_421.name_stat)
        idx_zero = findfirst(==(0.0), stats_421.vars.β_MOT)
        nums_421 = selectdim(stats_421.num_stat, idx_bias, idx_zero)
        stds_421 = selectdim(stats_421.std_num_stat, idx_bias, idx_zero)
        max_num_421 = max_number_with_error(nums_421, stds_421)
        max_time_421 = maximum(stats_421.vars.t_load) * γ_active

        runinfo_626 = only(read_num_evol_runinfos(
            joinpath(path_data, "MOT loading 626"), pair;
            var_specs=DUALMOT_LOADCFG_VAR_SPECS,
            validate_vars=validate_loadcfg_comparison_vars))
        stats_626 = [calc_num_evol_block(data;
            label="preprocess 626 $pair data[$idx]", bounds_sigmax_num,
            bounds_sigmay_num, num_max_num)
            for (idx, data) in enumerate(runinfo_626.data)]
        max_num_626 = maximum(max_number_with_error(stats.num_stat,
            stats.std_num_stat) for stats in stats_626)
        max_time_626 = maximum(maximum(stats.vars.t_load) for stats in stats_626) * γ_active

        xmax = max(max_time_421, max_time_626)
        ymax = max(max_num_421, max_num_626) / 1e7
        limits[pair] = (x=(0.0, 1.03 * xmax), y=(0.0, 1.05 * ymax))
        println("$pair shared linear loading limits: $(limits[pair])")
    end
    limits
end

limits_loading_linear = preprocess_loading_limits()
