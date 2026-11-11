# Reuse production detection/fitting without running the spectrum exports.
if !isdefined(@__MODULE__, :Spectral626Calibration)
    @eval module Spectral626Calibration
        const ARGS = ["--definitions-only"]
        include(joinpath(@__DIR__, "spectral626.jl"))
    end
end
if !isdefined(@__MODULE__, :CalibrationBenchmark)
    include(joinpath(@__DIR__, "..", "src", "calibration_benchmark.jl"))
end
using .CalibrationBenchmark
using .Spectral626Calibration: load_calibration_references,
    fit_spline_with_offset_and_error, interpolate_natural_residual,
    calc_spline_calibrated_frequency, calc_spline_calibrated_slope,
    calc_spline_value_with_linear_tails, calc_spline_slope_with_linear_tails
using CSV, DataFrames, JLD2, Printf, YAML
import CairoMakie as CM

function fit_benchmark_calibration(entry::AbstractDict, references::NamedTuple)
    (; points_reference, ids_file) = references
    model = Symbol(entry["model"])
    if model == :spline_with_offset_and_error
        fit = fit_spline_with_offset_and_error(points_reference, ids_file,
            Float64.(entry["knots"]), Float64(entry["regularisation"]);
            fit_source_offsets=get(entry, "fit_source_offsets", true))
        return (; tag=String(entry["tag"]), model, points_reference, ids_file, fit...)
    end
    model == :naturalcubicspline ||
        throw(ArgumentError("Benchmark model must be spline_with_offset_and_error or naturalcubicspline; got $model"))
    splines_by_file = Dict{Int,Any}()
    bounds_reference_by_file = Dict{Int,Tuple{Float64,Float64}}()
    for id_file in ids_file
        points_source = sort(filter(point -> point.id_file == id_file, points_reference);
            by=point -> point.time_peak_s)
        times_s = [point.time_peak_s for point in points_source]
        length(times_s) >= 3 && all(>(0), diff(times_s)) ||
            throw(ArgumentError("Natural cubic source $id_file needs at least three distinct times"))
        splines_by_file[id_file] = interpolate_natural_residual(times_s,
            [point.freq_mhz for point in points_source])
        bounds_reference_by_file[id_file] = extrema(times_s)
    end
    return (; tag=String(entry["tag"]), model, points_reference, ids_file,
        splines_by_file, bounds_reference_by_file)
end

function calc_benchmark_frequency(calibration::NamedTuple, id_file::Integer, time_s::Real)::Float64
    calibration.model == :spline_with_offset_and_error &&
        return calc_spline_calibrated_frequency(calibration, id_file, time_s)
    return calc_spline_value_with_linear_tails(calibration.splines_by_file[id_file], time_s,
        calibration.bounds_reference_by_file[id_file])
end

function calc_benchmark_slope(calibration::NamedTuple, id_file::Integer, time_s::Real)::Float64
    calibration.model == :spline_with_offset_and_error &&
        return calc_spline_calibrated_slope(calibration, id_file, time_s)
    return calc_spline_slope_with_linear_tails(calibration.splines_by_file[id_file], time_s,
        calibration.bounds_reference_by_file[id_file])
end

function find_benchmark_reference_time(calibration::NamedTuple, id_file::Integer, frequency_mhz::Real)::Float64
    time_left_s, time_right_s = calibration.bounds_reference_by_file[id_file]
    freq_left_mhz = calc_benchmark_frequency(calibration, id_file, time_left_s)
    freq_right_mhz = calc_benchmark_frequency(calibration, id_file, time_right_s)
    freq_left_mhz > freq_right_mhz || throw(ArgumentError("Range-reference frequency must decrease with time"))
    if frequency_mhz > freq_left_mhz || frequency_mhz < freq_right_mhz
        time_endpoint_s = frequency_mhz > freq_left_mhz ? time_left_s : time_right_s
        freq_endpoint_mhz = calc_benchmark_frequency(calibration, id_file, time_endpoint_s)
        slope = calc_benchmark_slope(calibration, id_file, time_endpoint_s)
        slope < 0 || throw(ArgumentError("Range-reference tail must have negative slope"))
        return time_endpoint_s + (frequency_mhz - freq_endpoint_mhz) / slope
    end
    for _ in 1:60
        time_mid_s = (time_left_s + time_right_s) / 2
        if calc_benchmark_frequency(calibration, id_file, time_mid_s) > frequency_mhz
            time_left_s = time_mid_s
        else
            time_right_s = time_mid_s
        end
    end
    return (time_left_s + time_right_s) / 2
end

function run_spectral626_benchmark(
    path_configs::AbstractString=joinpath(@__DIR__, "spectrum626_benchmark");
    path_output::AbstractString=joinpath(path_configs, "results"),
    source_cache::AbstractDict=Dict{Any,Any}(),
)
    paths_config = sort(filter(isfile,
        [joinpath(path, "config.yaml") for path in readdir(path_configs; join=true) if isdir(path)]))
    length(paths_config) >= 2 || throw(ArgumentError("Need at least two variant/config.yaml files in $path_configs"))
    configs = YAML.load_file.(paths_config)
    config_benchmark = first(configs)["benchmark"]
    path_data = abspath(String(config_benchmark["path_data"]))
    frequency_limits_mhz = Float64.(config_benchmark["frequency_range_mhz"])
    length(frequency_limits_mhz) == 2 && all(isfinite, frequency_limits_mhz) &&
        first(frequency_limits_mhz) < last(frequency_limits_mhz) ||
        throw(ArgumentError("frequency_range_mhz must contain two increasing finite limits"))
    id_range_reference = Int(config_benchmark["range_reference_index"])
    n_plot_samples = Int(get(config_benchmark, "plot_samples", 2001))
    n_plot_samples >= 2 || throw(ArgumentError("plot_samples must be at least two"))
    tags = String[entry["tag"] for entry in first(configs)["calibration"]]
    length(unique(tags)) == length(tags) || throw(ArgumentError("Calibration tags must be unique"))
    names_variant = String[config["benchmark"]["name"] for config in configs]
    length(unique(names_variant)) == length(names_variant) || throw(ArgumentError("Variant names must be unique"))
    for config in configs
        config["benchmark"] == config_benchmark ||
            all(get(config["benchmark"], key, nothing) == get(config_benchmark, key, nothing)
                for key in ("path_data", "frequency_range_mhz", "range_reference_index", "plot_samples")) ||
            throw(ArgumentError("All variants must use the same data, range reference, limits and plot sample count"))
        String[entry["tag"] for entry in config["calibration"]] == tags ||
            throw(ArgumentError("All variants must have the same ordered calibration tags"))
    end

    summaries = NamedTuple[]
    curves = NamedTuple[]
    references_export = NamedTuple[]
    cases = Dict{Tuple{String,String},Any}()
    fig = CM.Figure(size=(1100, 430length(tags)), fontsize=17)
    for (idx_tag, tag) in enumerate(tags)
        entry_first = first(configs)["calibration"][idx_tag]
        signature(entry) = [(Int(source["index"]),
            [(Float64(peak["freq"]), Float64(peak["time"])) for peak in source["peaks"]])
            for source in entry["source"]]
        all(signature(config["calibration"][idx_tag]) == signature(entry_first) for config in configs) ||
            throw(ArgumentError("Variants for '$tag' must use identical ordered sources and reference peaks"))
        references = load_calibration_references(path_data, entry_first["source"]; source_cache)
        id_range_reference in references.ids_file || throw(ArgumentError("Range-reference source $id_range_reference absent from '$tag'"))
        append!(references_export, [merge((; tag), point) for point in references.points_reference])
        entry_natural = merge(entry_first, Dict("model" => "naturalcubicspline"))
        calibration_reference = fit_benchmark_calibration(entry_natural, references)
        bounds_time_s = extrema(find_benchmark_reference_time(calibration_reference,
            id_range_reference, freq) for freq in frequency_limits_mhz)
        times_plot_s = collect(range(bounds_time_s...; length=n_plot_samples))
        all(t -> calc_benchmark_slope(calibration_reference, id_range_reference, t) < 0, times_plot_s) ||
            throw(ArgumentError("Natural range-reference curve must be strictly decreasing in the selected span"))
        frequencies_plot_mhz = [calc_benchmark_frequency(calibration_reference, id_range_reference, t) for t in times_plot_s]
        ax = CM.Axis(fig[idx_tag, 1]; title="$tag: maximum source disagreement",
            xlabel="Fixed natural reference frequency (MHz), source $id_range_reference",
            ylabel="Maximum − minimum calibration (MHz)",
            xgridvisible=false, ygridvisible=false, xticksmirrored=true, yticksmirrored=true,
            xtickalign=1, ytickalign=1)
        for (idx_variant, config) in enumerate(configs)
            name_variant = names_variant[idx_variant]
            entry = config["calibration"][idx_tag]
            calibration = fit_benchmark_calibration(entry, references)
            calc_frequency = (id, t) -> calc_benchmark_frequency(calibration, id, t)
            calc_slope = (id, t) -> calc_benchmark_slope(calibration, id, t)
            times_break = [bounds_time_s...,
                (point.time_peak_s for point in references.points_reference)...]
            calibration.model == :spline_with_offset_and_error && append!(times_break, Float64.(entry["knots"]))
            times_break = sort!(unique(filter(t -> first(bounds_time_s) <= t <= last(bounds_time_s), times_break)))
            result = calc_max_curve_disagreement(calc_frequency, calc_slope, references.ids_file, times_break)
            mat_frequency_mhz = [calc_frequency(id, t) for t in times_plot_s, id in references.ids_file]
            disagreement_mhz = vec(maximum(mat_frequency_mhz; dims=2) - minimum(mat_frequency_mhz; dims=2))
            maximum(disagreement_mhz) <= result.disagreement_mhz + 1e-8 ||
                error("Continuous maximum is smaller than a sampled disagreement")
            reference_error_mhz = maximum(abs(calc_frequency(point.id_file, point.time_peak_s) - point.freq_mhz)
                for point in references.points_reference)
            reference_error_mhz < 1e-8 || error("Reference interpolation failed for $tag / $name_variant")
            offsets_enabled = calibration.model == :spline_with_offset_and_error && calibration.fit_source_offsets
            # This catches accidental changes to the no-offset solve.
            calibration.model == :spline_with_offset_and_error && !offsets_enabled &&
                !all(iszero, values(calibration.offset_by_file)) && error("No-offset fit returned a nonzero offset")
            frequency_at_max_mhz = calc_benchmark_frequency(calibration_reference, id_range_reference, result.time_s)
            span_extrapolated = any(id -> !(calibration.bounds_reference_by_file[id][1] <= result.time_s <=
                calibration.bounds_reference_by_file[id][2]), (result.id_upper, result.id_lower))
            push!(summaries, (; tag, variant=name_variant, model=String(calibration.model), offsets_enabled,
                n_sources=length(references.ids_file), n_references=length(references.points_reference),
                frequency_min_mhz=first(frequency_limits_mhz), frequency_max_mhz=last(frequency_limits_mhz),
                time_min_s=first(bounds_time_s), time_max_s=last(bounds_time_s),
                max_disagreement_mhz=result.disagreement_mhz, frequency_at_max_mhz,
                time_at_max_s=result.time_s, result.id_upper, result.id_lower,
                maximizing_pair_extrapolated=span_extrapolated, max_reference_error_mhz=reference_error_mhz))
            for (idx_time, time_s) in enumerate(times_plot_s)
                push!(curves, (; tag, variant=name_variant, time_s,
                    frequency_reference_mhz=frequencies_plot_mhz[idx_time], disagreement_mhz=disagreement_mhz[idx_time]))
            end
            cases[(tag, name_variant)] = (; calibration, result, times_break, bounds_time_s,
                calibration_reference, times_plot_s, mat_frequency_mhz, disagreement_mhz)
            label = @sprintf("%s: %.6f MHz", name_variant, result.disagreement_mhz)
            CM.lines!(ax, frequencies_plot_mhz, disagreement_mhz; label, linewidth=2)
            CM.scatter!(ax, [frequency_at_max_mhz], [result.disagreement_mhz]; markersize=8,
                color=CM.Makie.wong_colors()[mod1(idx_variant, 7)])
            @printf("%s / %s: maximum %.9f MHz at %.6f MHz, t=%.9f s; sources %d/%d; reference error %.3g MHz\n",
                tag, name_variant, result.disagreement_mhz, frequency_at_max_mhz, result.time_s,
                result.id_upper, result.id_lower, reference_error_mhz)
        end
        CM.xlims!(ax, frequency_limits_mhz...)
        CM.axislegend(ax; position=:lt, framevisible=false, backgroundcolor=:transparent, labelsize=13)
    end
    mkpath(path_output)
    CSV.write(joinpath(path_output, "summary.csv"), DataFrame(summaries))
    CSV.write(joinpath(path_output, "disagreement_curves.csv"), DataFrame(curves))
    CSV.write(joinpath(path_output, "detected_references.csv"), DataFrame(references_export))
    # Plain arrays make the saved benchmark independent of Julia spline serialization.
    frequencies_export = [(; tag, variant, time_s=case.times_plot_s,
        ids_source=case.calibration.ids_file, frequency_mhz=case.mat_frequency_mhz)
        for ((tag, variant), case) in cases]
    JLD2.jldsave(joinpath(path_output, "calibration_curves.jld2"); frequencies_export, summaries)
    CM.save(joinpath(path_output, "disagreement.png"), fig; px_per_unit=2)
    for (path_config, name_variant) in zip(paths_config, names_variant)
        path_snapshot = joinpath(path_output, "config_snapshots", basename(dirname(path_config)))
        mkpath(path_snapshot)
        cp(path_config, joinpath(path_snapshot, "config.yaml"); force=true)
    end
    println("Saved benchmark to $(abspath(path_output)); source cache contains $(length(source_cache)) captures")
    return (; summaries, cases, fig, source_cache, path_output=abspath(path_output))
end

if !("--definitions-only" in ARGS)
length(ARGS) <= 2 || throw(ArgumentError("Usage: benchmark_spectral626.jl [config-root] [output-folder]"))
path_benchmark_configs = isempty(ARGS) ? joinpath(@__DIR__, "spectrum626_benchmark") : abspath(ARGS[1])
path_benchmark_output = length(ARGS) == 2 ? abspath(ARGS[2]) : joinpath(path_benchmark_configs, "results")
if !isdefined(@__MODULE__, :benchmark_source_cache)
    benchmark_source_cache = Dict{Any,Any}()
end
benchmark_data = run_spectral626_benchmark(path_benchmark_configs;
    path_output=path_benchmark_output, source_cache=benchmark_source_cache)
end
nothing
