using CSV
using GLMakie
using Colors
using DataFrames
using ImageFiltering
using JLD2
using LsqFit
using MAT
using Peaks
using Statistics
using YAML

GLMakie.activate!()
const path_data_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\spectrum 626"
const n_horizontal_divisions = 14
const val_smoothing_sigma_points = 2.5
const n_smoothing_window_points = 11 # symmetric kernel spanning ten sample intervals
const bounds_spectrum_frequency_mhz = (-800.0, 800.0)
const bounds_result_frequency_mhz = (-720.0, 720.0)

function smooth_scope_channel(values::AbstractVector{<:Real})::Vector{Float64}
    kernel_gaussian = Kernel.gaussian(
        (val_smoothing_sigma_points,), (n_smoothing_window_points,))
    return Float64.(imfilter(values, kernel_gaussian, Pad(:replicate)))
end

function read_scope_metadata(path_csv::AbstractString, n_lines::Int=10)::Vector{String}
    open(path_csv, "r") do io
        return [readline(io) for _ in 1:n_lines]
    end
end

function parse_record_length(metadata_lines::AbstractVector{<:AbstractString})::Union{Nothing,Int}
    idx_record = findfirst(line -> startswith(line, "Record Length,"), metadata_lines)
    isnothing(idx_record) && return nothing
    line_record = metadata_lines[idx_record]
    match_record = match(r"Record Length,Analog:(\d+)", line_record)
    isnothing(match_record) && throw(ArgumentError("Could not parse record length from metadata: $line_record"))
    return parse(Int, match_record.captures[1])
end

function convert_scope_csv_to_jld2(path_csv::AbstractString, path_jld2::AbstractString)::Nothing
    isfile(path_csv) || throw(ArgumentError("Input CSV does not exist: $path_csv"))
    metadata_lines = read_scope_metadata(path_csv)

    # Lines 1–10 are scope metadata; lines 11–12 label the three numeric columns.
    table = CSV.read(path_csv, DataFrame;
        header=false,
        skipto=13,
        types=Float64,
        ntasks=1,
        strict=true,
    )
    size(table, 2) == 3 || throw(DimensionMismatch("Expected 3 numeric columns (time, CH1, CH2), got $(size(table, 2))"))
    n_expected = parse_record_length(metadata_lines)
    isnothing(n_expected) || nrow(table) == n_expected ||
        throw(DimensionMismatch("Record Length metadata says $n_expected samples, CSV contains $(nrow(table))"))

    mkpath(dirname(abspath(path_jld2)))
    JLD2.jldsave(path_jld2;
        time_s=table[!, 1],
        ch1_v=table[!, 2],
        ch2_v=table[!, 3],
        column_names=["time", "CH1", "CH2"],
        column_units=["s", "V", "V"],
        metadata_lines,
        source_file=basename(path_csv),
    )

    println("Saved $(nrow(table)) samples to $path_jld2")
    return nothing
end

function convert_missing_scope_csvs(path_data::AbstractString)::Nothing
    isdir(path_data) || throw(ArgumentError("Data folder does not exist: $path_data"))
    paths_csv = filter(path -> occursin(r"^SDS\d{5}\.csv$", basename(path)), readdir(path_data; join=true))
    sort!(paths_csv)
    println("Found $(length(paths_csv)) scope CSV files in $path_data")

    n_skipped = 0
    n_converted = 0
    for path_csv in paths_csv
        path_jld2 = splitext(path_csv)[1] * ".jld2"
        if isfile(path_jld2)
            println("Skipping existing cache: $path_jld2")
            n_skipped += 1
            continue
        end

        convert_scope_csv_to_jld2(path_csv, path_jld2)
        n_converted += 1
        GC.gc()
    end

    println("Scope conversion complete: converted $n_converted; skipped $n_skipped existing caches")
    return nothing
end

function find_plateau_extrema(y::AbstractVector{<:Real}, find_minima::Bool)::Vector{Int}
    idx_extrema = Int[]
    n_samples = length(y)
    idx_run_start = 1

    while idx_run_start <= n_samples
        idx_run_end = idx_run_start
        while idx_run_end < n_samples && y[idx_run_end + 1] == y[idx_run_start]
            idx_run_end += 1
        end

        if idx_run_start > 1 && idx_run_end < n_samples
            val = y[idx_run_start]
            val_before = y[idx_run_start - 1]
            val_after = y[idx_run_end + 1]
            is_extremum = find_minima ? val < val_before && val < val_after : val > val_before && val > val_after
            is_extremum && push!(idx_extrema, idx_run_start)
        end
        idx_run_start = idx_run_end + 1
    end

    return idx_extrema
end

function select_deadtime_extrema(
    idx_candidates::AbstractVector{<:Integer},
    y::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    find_minima::Bool,
    dead_time_s::Real,
)::Vector{Int}
    length(y) == length(time_s) || throw(DimensionMismatch("Signal and time arrays must have equal lengths"))
    dead_time_s > 0 || throw(ArgumentError("dead_time_s must be positive; got $dead_time_s"))

    idx_selected = Int[]
    pos_candidate = 1
    while pos_candidate <= length(idx_candidates)
        idx_anchor = idx_candidates[pos_candidate]
        time_limit = time_s[idx_anchor] + dead_time_s
        pos_group_end = pos_candidate
        while pos_group_end < length(idx_candidates) && time_s[idx_candidates[pos_group_end + 1]] <= time_limit
            pos_group_end += 1
        end

        # Keep the strongest candidate in the dead-time window; strict comparison keeps the first tie.
        idx_best = idx_candidates[pos_candidate]
        for pos in (pos_candidate + 1):pos_group_end
            idx_candidate = idx_candidates[pos]
            is_better = find_minima ? y[idx_candidate] < y[idx_best] : y[idx_candidate] > y[idx_best]
            is_better && (idx_best = idx_candidate)
        end
        push!(idx_selected, idx_best)

        pos_candidate = pos_group_end + 1
        while pos_candidate <= length(idx_candidates) && time_s[idx_candidates[pos_candidate]] - time_s[idx_best] < dead_time_s
            pos_candidate += 1
        end
    end

    return idx_selected
end

function find_windowed_extrema(
    y::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    find_minima::Bool,
    dead_time_s::Real,
)::Vector{Int}
    length(y) == length(time_s) || throw(DimensionMismatch("Signal and time arrays must have equal lengths"))
    dead_time_s > 0 || throw(ArgumentError("dead_time_s must be positive; got $dead_time_s"))

    bin_width_s = dead_time_s / 4
    n_bins = ceil(Int, (time_s[end] - time_s[1]) / bin_width_s)
    val_bin_mean = Float64[]
    idx_bin_first = Int[]
    idx_bin_last = Int[]
    sizehint!(val_bin_mean, n_bins)
    sizehint!(idx_bin_first, n_bins)
    sizehint!(idx_bin_last, n_bins)

    for idx_bin in 1:n_bins
        time_first = time_s[1] + (idx_bin - 1) * bin_width_s
        time_last = min(time_first + bin_width_s, time_s[end])
        idx_first = searchsortedfirst(time_s, time_first)
        idx_last = searchsortedlast(time_s, time_last)
        idx_first <= idx_last || continue
        push!(idx_bin_first, idx_first)
        push!(idx_bin_last, idx_last)
        push!(val_bin_mean, sum(@view y[idx_first:idx_last]) / (idx_last - idx_first + 1))
    end

    idx_bin_candidates = find_plateau_extrema(val_bin_mean, find_minima)
    idx_candidates = Int[]
    for idx_bin in idx_bin_candidates
        idx_first = idx_bin_first[idx_bin]
        idx_last = idx_bin_last[idx_bin]
        val_window = @view y[idx_first:idx_last]
        _, pos_extremum = find_minima ? findmin(val_window) : findmax(val_window)
        push!(idx_candidates, idx_first + pos_extremum - 1)
    end

    return select_deadtime_extrema(idx_candidates, y, time_s, find_minima, dead_time_s)
end

function find_local_peaks_with_deadtime(
    y::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    dead_time_s::Real,
    baseline_v::Real,
    min_prominence_v::Real,
)::Vector{Int}
    length(y) == length(time_s) || throw(DimensionMismatch("Signal and time arrays must have equal lengths"))
    dead_time_s > 0 || throw(ArgumentError("dead_time_s must be positive; got $dead_time_s"))
    min_prominence_v >= 0 || throw(ArgumentError("min_prominence_v must be nonnegative; got $min_prominence_v"))
    idx_candidates = filter!(
        idx_peak -> y[idx_peak] - baseline_v > min_prominence_v,
        findmaxima(y).indices,
    )
    idx_selected = select_deadtime_extrema(idx_candidates, y, time_s, false, dead_time_s)
    return idx_selected
end

function calc_mean_excluding_top_fraction(y::AbstractVector{<:Real}, fraction_top::Real)::Float64
    0 <= fraction_top < 1 || throw(ArgumentError("fraction_top must be in [0, 1); got $fraction_top"))
    n_keep = floor(Int, (1 - fraction_top) * length(y))
    n_keep > 0 || throw(ArgumentError("Not enough samples to calculate a trimmed mean"))
    val_sorted = Float64.(y)
    partialsort!(val_sorted, n_keep)
    return sum(@view val_sorted[1:n_keep]) / n_keep
end

function write_ch1_periods_csv(
    path_csv::AbstractString,
    time_s::AbstractVector{<:Real},
    ch1_v::AbstractVector{<:Real},
    idx_minima::AbstractVector{<:Integer},
    idx_maxima::AbstractVector{<:Integer},
)::Nothing
    period_rows = NamedTuple[]
    pos_max = 1

    for idx_period in 1:(length(idx_minima) - 1)
        idx_min_start = idx_minima[idx_period]
        idx_min_end = idx_minima[idx_period + 1]
        while pos_max <= length(idx_maxima) && idx_maxima[pos_max] <= idx_min_start
            pos_max += 1
        end
        pos_max <= length(idx_maxima) && idx_maxima[pos_max] < idx_min_end ||
            throw(ArgumentError("No maximum found between minima at sample indices $idx_min_start and $idx_min_end"))
        idx_max = idx_maxima[pos_max]
        pos_max += 1

        push!(period_rows, (
            period=idx_period,
            start_min_idx=idx_min_start,
            start_min_time_s=time_s[idx_min_start],
            start_min_ch1_v=ch1_v[idx_min_start],
            peak_idx=idx_max,
            peak_time_s=time_s[idx_max],
            peak_ch1_v=ch1_v[idx_max],
            end_min_idx=idx_min_end,
            end_min_time_s=time_s[idx_min_end],
            end_min_ch1_v=ch1_v[idx_min_end],
            period_s=time_s[idx_min_end] - time_s[idx_min_start],
            period_n_datapoints=idx_min_end - idx_min_start + 1,
            uphill_s=time_s[idx_max] - time_s[idx_min_start],
            uphill_n_datapoints=idx_max - idx_min_start + 1,
            downhill_s=time_s[idx_min_end] - time_s[idx_max],
            downhill_n_datapoints=idx_min_end - idx_max + 1,
        ))
    end

    if isempty(period_rows)
        println("No completed CH1 periods found; skipping period CSV $path_csv")
        return nothing
    end
    mkpath(dirname(abspath(path_csv)))
    CSV.write(path_csv, DataFrame(period_rows))
    println("Saved $(length(period_rows)) completed periods to $path_csv")
    return nothing
end

function find_first_uphill_zero_crossing(
    y::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    idx_start::Int,
    idx_end::Int,
)::Union{Nothing,Float64}
    idx_start < idx_end || throw(ArgumentError("A rising ramp must contain at least two samples"))
    for idx in idx_start:(idx_end - 1)
        val = y[idx]
        val_next = y[idx + 1]
        iszero(val) && return Float64(time_s[idx])
        if val < 0 && val_next >= 0
            frac = -val / (val_next - val)
            return Float64(time_s[idx] + frac * (time_s[idx + 1] - time_s[idx]))
        end
    end
    return nothing
end

function find_ch1_rising_ramp(
    ch1_v::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    idx_minima::AbstractVector{<:Integer},
    idx_maxima::AbstractVector{<:Integer},
)
    length(ch1_v) == length(time_s) || throw(DimensionMismatch("CH1 and time arrays must have equal lengths"))
    best_ramp = nothing
    for (pos_min, idx_start) in enumerate(idx_minima)
        pos_max = searchsortedfirst(idx_maxima, idx_start + 1)
        pos_max <= length(idx_maxima) || continue
        idx_end = idx_maxima[pos_max]
        pos_min == length(idx_minima) || idx_end < idx_minima[pos_min + 1] || continue
        time_zero_s = find_first_uphill_zero_crossing(ch1_v, time_s, idx_start, idx_end)
        isnothing(time_zero_s) && continue
        if isnothing(best_ramp) || abs(time_zero_s) < abs(best_ramp.time_zero_s)
            best_ramp = (; idx_start, idx_end, time_zero_s)
        end
    end
    isnothing(best_ramp) && throw(ArgumentError(
        "No complete CH1 rising ramp with a zero crossing was found between a minimum and a later maximum"))
    return best_ramp
end

function find_ch1_alignment_zero_crossing(
    ch1_v::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    dead_time_s::Real=2.0,
)::Float64
    length(ch1_v) == length(time_s) || throw(DimensionMismatch("CH1 and time arrays must have equal lengths"))
    idx_minima = find_windowed_extrema(ch1_v, time_s, true, dead_time_s)
    idx_maxima = find_windowed_extrema(ch1_v, time_s, false, dead_time_s)
    return find_ch1_rising_ramp(ch1_v, time_s, idx_minima, idx_maxima).time_zero_s
end

function oklch_to_rgb(h_deg::Real, l::Real, c::Real)::RGB{Float64}
    angle = deg2rad(h_deg)
    a = c * cos(angle)
    b = c * sin(angle)
    l_ = l + 0.3963377774 * a + 0.2158037573 * b
    m_ = l - 0.1055613458 * a - 0.0638541728 * b
    s_ = l - 0.0894841775 * a - 1.2914855480 * b
    l_cubed, m_cubed, s_cubed = l_^3, m_^3, s_^3

    red_linear = 4.0767416621 * l_cubed - 3.3077115913 * m_cubed + 0.2309699292 * s_cubed
    green_linear = -1.2684380046 * l_cubed + 2.6097574011 * m_cubed - 0.3413193965 * s_cubed
    blue_linear = -0.0041960863 * l_cubed - 0.7034186147 * m_cubed + 1.7076147010 * s_cubed
    to_srgb(channel) = clamp(channel <= 0.0031308 ? 12.92 * channel : 1.055 * channel^(1 / 2.4) - 0.055, 0, 1)
    return RGB(to_srgb(red_linear), to_srgb(green_linear), to_srgb(blue_linear))
end

function parse_scope_metadata_value(metadata_lines::AbstractVector{<:AbstractString}, field::AbstractString)::String
    idx_field = findfirst(line -> startswith(line, "$field,"), metadata_lines)
    isnothing(idx_field) && throw(ArgumentError("CSV metadata is missing '$field'"))
    fields = split(metadata_lines[idx_field], ',')
    length(fields) >= 2 || throw(ArgumentError("Could not parse '$field' metadata: $(metadata_lines[idx_field])"))
    return strip(fields[2])
end

function calc_time_axis(metadata_lines::AbstractVector{<:AbstractString}, n_samples::Int)::AbstractRange{Float64}
    units = parse_scope_metadata_value(metadata_lines, "Horizontal Units")
    unit_scale = units == "us" ? 1e-6 : units == "ms" ? 1e-3 : units == "s" ? 1.0 :
        throw(ArgumentError("Unsupported Horizontal Units '$units'; expected us, ms, or s"))
    scale_str = parse_scope_metadata_value(metadata_lines, "Horizontal Scale")
    scale_per_div = parse(Float64, scale_str) * unit_scale
    half_span_s = (n_horizontal_divisions / 2) * scale_per_div
    return range(-half_span_s, stop=half_span_s, length=n_samples)
end

function load_scope_period_data(path_jld2::AbstractString, dead_time_s::Real)
    payload = JLD2.load(path_jld2)
    ch1_v = smooth_scope_channel(payload["ch1_v"])
    ch2_v = smooth_scope_channel(payload["ch2_v"])
    metadata_lines = payload["metadata_lines"]
    length(ch1_v) == length(ch2_v) ||
        throw(DimensionMismatch("CH1 and CH2 arrays must have equal lengths in $path_jld2"))
    time_s = calc_time_axis(metadata_lines, length(ch1_v))
    idx_minima = find_windowed_extrema(ch1_v, time_s, true, dead_time_s)
    idx_maxima = find_windowed_extrema(ch1_v, time_s, false, dead_time_s)
    return (; time_s, ch1_v, ch2_v, idx_minima, idx_maxima)
end

function plot_ch2_rising_ramps(
    path_data::AbstractString,
    ids_file::AbstractVector{<:Integer},
    path_png::AbstractString;
    val_min_peak_prominence_v::Float64 = 0.005
)
    isempty(ids_file) && throw(ArgumentError("At least one spectrum file is required"))

    val_offset_per_file_v = 0.5
    val_fraction_top_excluded = 0.01
    val_peak_dead_time_s = 0.05
    val_hue_deg_file = [mod(24 + (idx - 1) * 360 / length(ids_file), 360)
        for idx in eachindex(ids_file)]
    val_lightness_ch2 = 0.60
    val_lightness_ch1 = 0.88
    val_chroma = 0.15
    colors_ch2 = [oklch_to_rgb(val_hue_deg_file[idx], val_lightness_ch2, val_chroma) for idx in eachindex(ids_file)]
    colors_ch1 = [oklch_to_rgb(val_hue_deg_file[idx], val_lightness_ch1, val_chroma) for idx in eachindex(ids_file)]

    px_width = 3840
    px_height = 2160
    fig = Figure(size=(px_width, px_height), fontsize=28)
    ax = Axis(fig[1, 1];
        title="CH1 and CH2 rising ramps, aligned to CH1 zero-crossings",
        xlabel="Time from first rising CH1 zero-crossing (s)",
        ylabel="CH1 and CH2 voltage (V), offset by 0.5 V per file",
        titlesize=36,
        xlabelsize=32,
        ylabelsize=32,
        xticklabelsize=26,
        yticklabelsize=26,
    )

    for (idx_file, id_file) in enumerate(ids_file)
        name_file = "SDS$(lpad(string(id_file), 5, '0'))"
        path_jld2 = joinpath(path_data, "$name_file.jld2")
        isfile(path_jld2) || throw(ArgumentError("JLD2 cache does not exist: $path_jld2"))
        (; time_s, ch1_v, ch2_v, idx_minima, idx_maxima) = load_scope_period_data(path_jld2, 2.0)

        (; idx_start, idx_end, time_zero_s) = find_ch1_rising_ramp(
            ch1_v, time_s, idx_minima, idx_maxima)
        val_offset_v = (length(ids_file) - idx_file) * val_offset_per_file_v
        time_ramp_offset = @view(time_s[idx_start:idx_end]) .- time_zero_s
        ch1_ramp_offset = (@view ch1_v[idx_start:idx_end]) ./ 8 .+ val_offset_v
        time_ramp = @view time_s[idx_start:idx_end]
        ch2_ramp = @view ch2_v[idx_start:idx_end]
        ch2_ramp_offset = ch2_ramp .+ val_offset_v
        lines!(ax, time_ramp_offset, ch1_ramp_offset;
            color=colors_ch1[idx_file], linewidth=1.5)
        lines!(ax, time_ramp_offset, ch2_ramp_offset;
            color=colors_ch2[idx_file], linewidth=1.5, label=name_file)
        val_ch2_baseline_v = calc_mean_excluding_top_fraction(ch2_ramp, val_fraction_top_excluded)
        idx_peaks_ch2 = find_local_peaks_with_deadtime(
            ch2_ramp, time_ramp, val_peak_dead_time_s, val_ch2_baseline_v, val_min_peak_prominence_v)
        for idx_peak in idx_peaks_ch2
            time_peak = time_ramp[idx_peak] - time_zero_s
            value_peak = ch2_ramp[idx_peak] + val_offset_v
            scatter!(ax, [time_peak], [value_peak];
                color=colors_ch2[idx_file], markersize=10, strokecolor=:black, strokewidth=1)
            text!(ax, time_peak, value_peak + 0.025;
                text=string(round(time_peak; digits=3)),
                color=colors_ch2[idx_file], fontsize=16,
                rotation=π / 2, align=(:center, :bottom),
                strokecolor=:white, strokewidth=0.4)
        end
        println("$name_file: selected rising ramp $(round(time_s[idx_start]; digits=6)) to $(round(time_s[idx_end]; digits=6)) s; rising zero-crossing $(round(time_zero_s; digits=6)) s; CH2 offset $val_offset_v V; trimmed CH2 baseline $(round(val_ch2_baseline_v; digits=4)) V; annotated $(length(idx_peaks_ch2)) peaks (>$val_min_peak_prominence_v V above baseline, $(val_peak_dead_time_s) s dead-time)")
        GC.gc()
    end

    axislegend(ax; position=:rt, framevisible=false, labelsize=24)
    save(path_png, fig; px_per_unit=1)
    println("Saved $(px_width)×$(px_height) multi-file CH1/CH2 PNG to $path_png")
    return fig
end

function find_local_peak_near_time(
    y::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    target_time_s::Real,
    time_uncertainty_s::Real,
)::Int
    length(y) == length(time_s) || throw(DimensionMismatch("Signal and time arrays must have equal lengths"))
    time_uncertainty_s > 0 || throw(ArgumentError("time_uncertainty_s must be positive"))
    idx_first = searchsortedfirst(time_s, target_time_s - time_uncertainty_s)
    idx_last = searchsortedlast(time_s, target_time_s + time_uncertainty_s)
    idx_first <= idx_last || throw(ArgumentError("No samples within ±$time_uncertainty_s s of $target_time_s s"))

    y_window = @view y[idx_first:idx_last]
    idx_local_peaks = findmaxima(y_window).indices
    isempty(idx_local_peaks) && throw(ArgumentError("No local CH2 peak within ±$time_uncertainty_s s of $target_time_s s"))
    idx_peaks = idx_first .+ idx_local_peaks .- 1
    return idx_peaks[argmax(y[idx_peaks])]
end

function is_missing_config_value(value)
    return ismissing(value) || isnothing(value) ||
        (value isa AbstractString && lowercase(strip(value)) == "missing")
end

function find_local_peak_near_frequency(
    y::AbstractVector{<:Real},
    frequency_mhz::AbstractVector{<:Real},
    target_frequency_mhz::Real,
    frequency_tolerance_mhz::Real,
)::Int
    length(y) == length(frequency_mhz) ||
        throw(DimensionMismatch("Signal and frequency arrays must have equal lengths"))
    frequency_tolerance_mhz > 0 || throw(ArgumentError("frequency_tolerance_mhz must be positive"))
    idx_window = findall(abs.(frequency_mhz .- target_frequency_mhz) .<= frequency_tolerance_mhz)
    isempty(idx_window) && throw(ArgumentError(
        "No samples within ±$frequency_tolerance_mhz MHz of $target_frequency_mhz MHz"))

    idx_first, idx_last = first(idx_window), last(idx_window)
    idx_peaks_local = findmaxima(@view y[idx_first:idx_last]).indices
    idx_peaks = filter(idx -> abs(frequency_mhz[idx] - target_frequency_mhz) <= frequency_tolerance_mhz,
        idx_first .+ idx_peaks_local .- 1)
    isempty(idx_peaks) && (idx_peaks = idx_window)
    return idx_peaks[argmax(y[idx_peaks])]
end

function calc_frequency_from_time(
    alpha_mhz_per_s::Real,
    val_p_shape::AbstractVector{<:Real},
    u_s::Real,
)::Float64
    return alpha_mhz_per_s * u_s *
        (1 + val_p_shape[1] * u_s + val_p_shape[2] * u_s^3 + val_p_shape[3] * u_s^4)
end

function calc_frequency_shape_slope(
    val_p_shape::AbstractVector{<:Real},
    u_s::Real,
)::Float64
    return 1 + 2val_p_shape[1] * u_s +
        4val_p_shape[2] * u_s^3 + 5val_p_shape[3] * u_s^4
end

function find_monotonic_model_limits(
    alpha_mhz_per_s::Real,
    val_p_shape::AbstractVector{<:Real},
    freq_min_mhz::Real,
    freq_max_mhz::Real,
)::Tuple{Float64,Float64}
    length(val_p_shape) == 3 || throw(DimensionMismatch("Expected three frequency-shape coefficients"))
    freq_min_mhz < 0 < freq_max_mhz ||
        throw(ArgumentError("Frequency limits must span zero; got $freq_min_mhz to $freq_max_mhz MHz"))
    iszero(alpha_mhz_per_s) && throw(ArgumentError("Calibration scale α cannot be zero"))

    step_s = 0.001
    max_steps = 100_000
    limits_u = Float64[]
    for direction in (-1, 1)
        u_s = 0.0
        for _ in 1:max_steps
            u_next_s = u_s + direction * step_s
            frequency_next_mhz = calc_frequency_from_time(alpha_mhz_per_s, val_p_shape, u_next_s)
            if calc_frequency_shape_slope(val_p_shape, u_next_s) <= 0 ||
                !(freq_min_mhz <= frequency_next_mhz <= freq_max_mhz)
                break
            end
            u_s = u_next_s
        end
        abs(u_s) < max_steps * step_s ||
            throw(ArgumentError("Calibration model did not reach a frequency limit or turning point within $(max_steps * step_s) s"))
        push!(limits_u, u_s)
    end
    return (limits_u[1], limits_u[2])
end

function find_monotonic_frequency_indices(
    time_relative_s::AbstractVector{<:Real},
    frequency_mhz::AbstractVector{<:Real},
    val_p_shape::AbstractVector{<:Real},
    freq_min_mhz::Real,
    freq_max_mhz::Real,
)::UnitRange{Int}
    length(time_relative_s) == length(frequency_mhz) ||
        throw(DimensionMismatch("Relative time and frequency arrays must have equal lengths"))
    length(val_p_shape) == 3 || throw(DimensionMismatch("Expected three frequency-shape coefficients"))
    isempty(time_relative_s) && throw(ArgumentError("Relative time array cannot be empty"))
    time_relative_s[1] <= 0 <= time_relative_s[end] ||
        throw(ArgumentError("Relative time array must span zero"))
    freq_min_mhz < 0 < freq_max_mhz ||
        throw(ArgumentError("Frequency limits must span zero; got $freq_min_mhz to $freq_max_mhz MHz"))

    idx_zero = searchsortedfirst(time_relative_s, 0)
    idx_first = idx_zero
    while idx_first > 1
        idx_next = idx_first - 1
        calc_frequency_shape_slope(val_p_shape, time_relative_s[idx_next]) > 0 || break
        freq_min_mhz <= frequency_mhz[idx_next] <= freq_max_mhz || break
        idx_first = idx_next
    end
    idx_last = idx_zero
    while idx_last < length(time_relative_s)
        idx_next = idx_last + 1
        calc_frequency_shape_slope(val_p_shape, time_relative_s[idx_next]) > 0 || break
        freq_min_mhz <= frequency_mhz[idx_next] <= freq_max_mhz || break
        idx_last = idx_next
    end
    return idx_first:idx_last
end

function fit_time_zero_and_scale(
    time_center_s::Real,
    time_peak_s::AbstractVector{<:Real},
    freq_mhz::AbstractVector{<:Real},
    val_p_shape::AbstractVector{<:Real},
)
    length(time_peak_s) == length(freq_mhz) ||
        throw(DimensionMismatch("Reference times and frequencies must have equal lengths"))
    length(time_peak_s) >= 2 || throw(ArgumentError("At least two references are needed to fit t₀ and α"))
    length(val_p_shape) == 3 || throw(DimensionMismatch("Expected shared p₂–p₄ coefficients"))
    val_time_centered_s = Float64.(time_peak_s .- time_center_s)

    function model_result(x, theta)
        return [calc_frequency_from_time(theta[2], val_p_shape, x[idx] - theta[1]) for idx in eachindex(x)]
    end
    function jacobian_result(x, theta)
        jac = zeros(Float64, length(x), 2)
        alpha_mhz_per_s = theta[2]
        for idx in eachindex(x)
            u_s = x[idx] - theta[1]
            shape = 1 + val_p_shape[1] * u_s + val_p_shape[2] * u_s^3 + val_p_shape[3] * u_s^4
            slope_shape = 1 + 2val_p_shape[1] * u_s +
                4val_p_shape[2] * u_s^3 + 5val_p_shape[3] * u_s^4
            jac[idx, 1] = -alpha_mhz_per_s * slope_shape
            jac[idx, 2] = u_s * shape
        end
        return jac
    end

    fit = LsqFit.curve_fit(
        model_result, jacobian_result, val_time_centered_s, Float64.(freq_mhz), [0.0, -472.0])
    fit.converged || @warn "Per-index t₀/α fit did not report convergence"
    return (time_zero_s=Float64(time_center_s + fit.param[1]),
        alpha_mhz_per_s=Float64(fit.param[2]), fit)
end

function fit_scale_with_fixed_time_zero(
    time_zero_s::Real,
    time_peak_s::AbstractVector{<:Real},
    freq_mhz::AbstractVector{<:Real},
    val_p_shape::AbstractVector{<:Real},
)
    length(time_peak_s) == length(freq_mhz) ||
        throw(DimensionMismatch("Reference times and frequencies must have equal lengths"))
    length(time_peak_s) >= 2 || throw(ArgumentError("At least two references are needed to fit α"))
    val_u_s = Float64.(time_peak_s .- time_zero_s)
    function model_scale(u_s, p)
        return [calc_frequency_from_time(p[1], val_p_shape, u) for u in u_s]
    end
    function jacobian_scale(u_s, p)
        return reshape([
            u * (1 + val_p_shape[1] * u + val_p_shape[2] * u^3 + val_p_shape[3] * u^4)
            for u in u_s
        ], :, 1)
    end
    fit = LsqFit.curve_fit(
        model_scale, jacobian_scale, val_u_s, Float64.(freq_mhz), [-472.0])
    fit.converged || @warn "Fixed-t₀ α fit did not report convergence"
    return (alpha_mhz_per_s=Float64(only(fit.param)), fit)
end

function fit_monotone_cubic(
    time_s::AbstractVector{<:Real}, freq_mhz::AbstractVector{<:Real};
    extrapolation::Symbol=:linear,
)
    extrapolation in (:linear, :quadratic) ||
        throw(ArgumentError("Spline extrapolation must be linear or quadratic; got $extrapolation"))
    length(time_s) == length(freq_mhz) ||
        throw(DimensionMismatch("Spline reference times and frequencies must have equal lengths"))
    n_required = extrapolation == :quadratic ? 3 : 2
    length(time_s) >= n_required ||
        throw(ArgumentError("$extrapolation spline extrapolation needs at least $n_required references"))
    ids_sorted = sortperm(time_s)
    val_time_s = Float64.(time_s[ids_sorted])
    val_freq_mhz = Float64.(freq_mhz[ids_sorted])
    val_h_s = diff(val_time_s)
    all(>(0), val_h_s) || throw(ArgumentError("Spline reference times must be distinct"))
    val_secant = diff(val_freq_mhz) ./ val_h_s
    all(>(0), val_secant) || all(<(0), val_secant) ||
        throw(ArgumentError("Spline reference frequencies must be strictly monotonic in time"))
    n_points = length(val_time_s)
    val_tangent = zeros(Float64, n_points)
    if n_points == 2
        val_tangent .= val_secant[1]
    else
        for idx in 2:(n_points - 1)
            h_left, h_right = val_h_s[idx - 1], val_h_s[idx]
            slope_left, slope_right = val_secant[idx - 1], val_secant[idx]
            weight_left = 2h_right + h_left
            weight_right = h_right + 2h_left
            val_tangent[idx] = (weight_left + weight_right) /
                (weight_left / slope_left + weight_right / slope_right)
        end
        val_tangent[1] = ((2val_h_s[1] + val_h_s[2]) * val_secant[1] -
            val_h_s[1] * val_secant[2]) / (val_h_s[1] + val_h_s[2])
        val_tangent[end] = ((2val_h_s[end] + val_h_s[end - 1]) * val_secant[end] -
            val_h_s[end] * val_secant[end - 1]) / (val_h_s[end] + val_h_s[end - 1])
        for idx in (1, n_points)
            idx_secant = idx == 1 ? 1 : length(val_secant)
            if sign(val_tangent[idx]) != sign(val_secant[idx_secant])
                val_tangent[idx] = 0.0
            elseif abs(val_tangent[idx]) > 3abs(val_secant[idx_secant])
                val_tangent[idx] = 3val_secant[idx_secant]
            end
        end
    end
    curvature_left = n_points >= 3 ?
        (val_secant[2] - val_secant[1]) / (val_time_s[3] - val_time_s[1]) : 0.0
    curvature_right = n_points >= 3 ?
        (val_secant[end] - val_secant[end - 1]) /
        (val_time_s[end] - val_time_s[end - 2]) : 0.0
    return (; time_s=val_time_s, freq_mhz=val_freq_mhz,
        tangent_mhz_per_s=val_tangent, extrapolation,
        curvature_tails=(curvature_left, curvature_right))
end

function calc_monotone_cubic_frequency(spline, time_s::Real)::Float64
    val_time_s = spline.time_s
    val_freq_mhz = spline.freq_mhz
    val_tangent = spline.tangent_mhz_per_s
    if time_s <= first(val_time_s)
        dt_s = time_s - first(val_time_s)
        curvature = spline.extrapolation == :quadratic ? spline.curvature_tails[1] : 0.0
        return first(val_freq_mhz) + first(val_tangent) * dt_s + curvature * dt_s^2
    elseif time_s >= last(val_time_s)
        dt_s = time_s - last(val_time_s)
        curvature = spline.extrapolation == :quadratic ? spline.curvature_tails[2] : 0.0
        return last(val_freq_mhz) + last(val_tangent) * dt_s + curvature * dt_s^2
    end
    idx = searchsortedlast(val_time_s, time_s)
    h_s = val_time_s[idx + 1] - val_time_s[idx]
    z = (time_s - val_time_s[idx]) / h_s
    return (2z^3 - 3z^2 + 1) * val_freq_mhz[idx] +
        (z^3 - 2z^2 + z) * h_s * val_tangent[idx] +
        (-2z^3 + 3z^2) * val_freq_mhz[idx + 1] +
        (z^3 - z^2) * h_s * val_tangent[idx + 1]
end

function find_monotone_cubic_tail_limit(spline, target_freq_mhz::Real, side::Symbol)::Float64
    side in (:left, :right) || throw(ArgumentError("Spline tail side must be left or right"))
    idx_endpoint = side == :left ? 1 : length(spline.time_s)
    direction = side == :left ? -1 : 1
    time_endpoint_s = spline.time_s[idx_endpoint]
    freq_endpoint_mhz = spline.freq_mhz[idx_endpoint]
    slope_endpoint = spline.tangent_mhz_per_s[idx_endpoint]
    curvature = spline.extrapolation == :quadratic ?
        spline.curvature_tails[side == :left ? 1 : 2] : 0.0
    movement_sign = sign(slope_endpoint * direction)
    iszero(movement_sign) && return time_endpoint_s
    distance_target_mhz = (target_freq_mhz - freq_endpoint_mhz) * movement_sign
    distance_target_mhz >= 0 ||
        throw(ArgumentError("$target_freq_mhz MHz lies inside the $side spline references"))
    iszero(distance_target_mhz) && return time_endpoint_s

    time_turn_s = iszero(curvature) ? time_endpoint_s :
        time_endpoint_s - slope_endpoint / (2curvature)
    if !iszero(curvature) && direction * (time_turn_s - time_endpoint_s) > 0
        freq_turn_mhz = calc_monotone_cubic_frequency(spline, time_turn_s)
        if distance_target_mhz >= (freq_turn_mhz - freq_endpoint_mhz) * movement_sign
            return time_turn_s
        end
        time_far_s = time_turn_s
    else
        span_s = max(1.0, 2distance_target_mhz / abs(slope_endpoint))
        time_far_s = time_endpoint_s + direction * span_s
        while (calc_monotone_cubic_frequency(spline, time_far_s) - target_freq_mhz) *
            movement_sign < 0
            span_s *= 2
            span_s <= 100_000 ||
                throw(ArgumentError("Spline tail did not reach $target_freq_mhz MHz"))
            time_far_s = time_endpoint_s + direction * span_s
        end
    end
    time_near_s = time_endpoint_s
    for _ in 1:60
        time_mid_s = (time_near_s + time_far_s) / 2
        if (calc_monotone_cubic_frequency(spline, time_mid_s) - target_freq_mhz) *
            movement_sign < 0
            time_near_s = time_mid_s
        else
            time_far_s = time_mid_s
        end
    end
    return (time_near_s + time_far_s) / 2
end

function find_monotone_cubic_model_limits(
    spline, freq_min_mhz::Real, freq_max_mhz::Real,
)::Tuple{Float64,Float64}
    freq_min_mhz < freq_max_mhz ||
        throw(ArgumentError("Frequency limits must increase; got $freq_min_mhz to $freq_max_mhz MHz"))
    is_descending = first(spline.freq_mhz) > last(spline.freq_mhz)
    freq_left_mhz = is_descending ? freq_max_mhz : freq_min_mhz
    freq_right_mhz = is_descending ? freq_min_mhz : freq_max_mhz
    freq_ref_min_mhz, freq_ref_max_mhz = extrema(spline.freq_mhz)
    time_left_s = freq_ref_min_mhz <= freq_left_mhz <= freq_ref_max_mhz ?
        find_monotone_cubic_time(spline, freq_left_mhz) :
        find_monotone_cubic_tail_limit(spline, freq_left_mhz, :left)
    time_right_s = freq_ref_min_mhz <= freq_right_mhz <= freq_ref_max_mhz ?
        find_monotone_cubic_time(spline, freq_right_mhz) :
        find_monotone_cubic_tail_limit(spline, freq_right_mhz, :right)
    return (time_left_s, time_right_s)
end

function find_monotone_cubic_time(spline, target_freq_mhz::Real)::Float64
    val_freq_mhz = spline.freq_mhz
    val_time_s = spline.time_s
    idx_knot = findfirst(==(target_freq_mhz), val_freq_mhz)
    isnothing(idx_knot) || return val_time_s[idx_knot]
    if (target_freq_mhz - val_freq_mhz[1]) * sign(val_freq_mhz[end] - val_freq_mhz[1]) <= 0
        time_target_s = find_monotone_cubic_tail_limit(spline, target_freq_mhz, :left)
        isapprox(calc_monotone_cubic_frequency(spline, time_target_s), target_freq_mhz; atol=1e-7) ||
            throw(ArgumentError("Left spline tail turns before $target_freq_mhz MHz"))
        return time_target_s
    elseif (target_freq_mhz - val_freq_mhz[end]) * sign(val_freq_mhz[end] - val_freq_mhz[1]) >= 0
        time_target_s = find_monotone_cubic_tail_limit(spline, target_freq_mhz, :right)
        isapprox(calc_monotone_cubic_frequency(spline, time_target_s), target_freq_mhz; atol=1e-7) ||
            throw(ArgumentError("Right spline tail turns before $target_freq_mhz MHz"))
        return time_target_s
    end
    idx = only(filter(i -> min(val_freq_mhz[i], val_freq_mhz[i + 1]) <= target_freq_mhz <=
        max(val_freq_mhz[i], val_freq_mhz[i + 1]), 1:(length(val_freq_mhz) - 1)))
    time_left_s, time_right_s = val_time_s[idx], val_time_s[idx + 1]
    for _ in 1:60
        time_mid_s = (time_left_s + time_right_s) / 2
        if (calc_monotone_cubic_frequency(spline, time_mid_s) - target_freq_mhz) *
            sign(val_freq_mhz[end] - val_freq_mhz[1]) < 0
            time_left_s = time_mid_s
        else
            time_right_s = time_mid_s
        end
    end
    return (time_left_s + time_right_s) / 2
end

function plot_monotone_cubic_calibration(
    path_png::AbstractString, tag_calibration::AbstractString,
    points_reference::AbstractVector, ids_file::AbstractVector{<:Integer},
    extrapolation::Symbol,
)
    splines_by_file = Dict{Int,Any}()
    time_zero_by_file = Dict{Int,Float64}()
    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1]; width=1200, height=850,
        title="Spectrum 626 time–frequency calibration: $tag_calibration",
        xlabel="Time from zero-frequency t₀ (s)", ylabel="Reference frequency (MHz)",
        titlesize=32, xlabelsize=28, ylabelsize=28,
        xticklabelsize=22, yticklabelsize=22)
    gl_text = GridLayout()
    fig[1, 2] = gl_text
    name_extrapolation = extrapolation == :linear ?
        "Linear endpoint-tangent extrapolation" :
        "Three-reference quadratic tails"
    Label(gl_text[1, 1], "Monotone cubic interpolation\n$name_extrapolation";
        color=:black, fontsize=22, halign=:left, valign=:top)
    marker_by_file = Dict(14 => :utriangle, 19 => :circle)
    for (idx_file, id_file) in enumerate(ids_file)
        ref_file = filter(point -> point.id_file == id_file, points_reference)
        spline = fit_monotone_cubic([point.time_peak_s for point in ref_file],
            [point.freq_mhz for point in ref_file]; extrapolation)
        splines_by_file[id_file] = spline
        t_zero_s = find_monotone_cubic_time(spline, 0.0)
        time_zero_by_file[id_file] = t_zero_s
        color_file = oklch_to_rgb(mod(24 + (idx_file - 1) * 360 / length(ids_file), 360), 0.60, 0.15)
        marker = get(marker_by_file, id_file, :circle)
        time_min_s, time_max_s = find_monotone_cubic_model_limits(
            spline, bounds_spectrum_frequency_mhz...)
        for (time_start_s, time_stop_s) in ((time_min_s, first(spline.time_s)),
            (last(spline.time_s), time_max_s))
            time_start_s < time_stop_s || continue
            val_time_plot_s = range(time_start_s, time_stop_s; length=400)
            lines!(ax, val_time_plot_s .- t_zero_s,
                [calc_monotone_cubic_frequency(spline, t) for t in val_time_plot_s];
                color=color_file, linewidth=2, linestyle=:dash)
        end
        val_time_plot_s = range(first(spline.time_s), last(spline.time_s); length=700)
        lines!(ax, val_time_plot_s .- t_zero_s,
            [calc_monotone_cubic_frequency(spline, t) for t in val_time_plot_s];
            color=color_file, linewidth=2,
            label="SDS$(lpad(string(id_file), 5, '0')) monotone cubic")
        scatter!(ax, [point.time_peak_s - t_zero_s for point in ref_file],
            [point.freq_mhz for point in ref_file];
            color=color_file, marker, markersize=24, strokecolor=:white, strokewidth=0.4)
        Label(gl_text[idx_file + 1, 1],
            "SDS$(lpad(string(id_file), 5, '0')): t₀=$(round(t_zero_s; digits=6)) s\n" *
            "endpoint slopes $(round(first(spline.tangent_mhz_per_s); digits=2)), " *
            "$(round(last(spline.tangent_mhz_per_s); digits=2)) MHz/s";
            color=color_file, fontsize=20, halign=:left, valign=:top)
        println("Calibration '$tag_calibration' SDS$(lpad(string(id_file), 5, '0')): monotone cubic through $(length(ref_file)) references; endpoint slopes $(first(spline.tangent_mhz_per_s)), $(last(spline.tangent_mhz_per_s)) MHz/s; t₀=$t_zero_s s")
    end
    Label(gl_text[length(ids_file) + 2, 1],
        "Solid: reference span; dashed: $extrapolation extrapolation";
        color=:black, fontsize=18, halign=:left, valign=:top)
    Legend(gl_text[length(ids_file) + 3, 1], ax;
        framevisible=false, labelsize=20, patchsize=(36, 24))
    resize_to_layout!(fig)
    mkpath(dirname(abspath(path_png)))
    save(path_png, fig; px_per_unit=2)
    println("Saved monotone cubic calibration PNG to $path_png")
    return (; tag=String(tag_calibration), model=:cubicspline, fig, points_reference, ids_file,
        extrapolation, splines_by_file, time_zero_by_file)
end

function plot_spectrum_frequency_calibration(
    path_data::AbstractString,
    path_png::AbstractString,
    entry_calibration::AbstractDict,
)
    tag_calibration = String(entry_calibration["tag"])
    model = Symbol(lowercase(String(get(entry_calibration, "model", "p4"))))
    model in (:p4, :cubicspline) ||
        throw(ArgumentError("Calibration '$tag_calibration' has unknown model '$model'; expected p4 or cubicspline"))
    entries_source = entry_calibration["source"]
    isempty(entries_source) && throw(ArgumentError("Calibration '$tag_calibration' has no sources"))

    time_uncertainty_s = 0.001
    points_reference = NamedTuple[]
    for entry in entries_source
        id_file = Int(entry["index"])
        name_file = "SDS$(lpad(string(id_file), 5, '0'))"
        path_jld2 = joinpath(path_data, "$name_file.jld2")
        isfile(path_jld2) || throw(ArgumentError("JLD2 cache does not exist: $path_jld2"))
        payload = JLD2.load(path_jld2)
        ch1_v = smooth_scope_channel(payload["ch1_v"])
        ch2_v = smooth_scope_channel(payload["ch2_v"])
        time_s = calc_time_axis(payload["metadata_lines"], length(ch2_v))
        time_alignment_zero_s = find_ch1_alignment_zero_crossing(ch1_v, time_s)
        time_aligned_s = time_s .- time_alignment_zero_s
        val_baseline_v = calc_mean_excluding_top_fraction(ch2_v, 0.01)

        for peak_config in entry["peaks"]
            freq_mhz = Float64(peak_config["freq"])
            time_config_s = Float64(peak_config["time"])
            idx_peak = find_local_peak_near_time(ch2_v, time_aligned_s, time_config_s, time_uncertainty_s)
            time_peak_s = Float64(time_aligned_s[idx_peak])
            time_peak_raw_s = Float64(time_s[idx_peak])
            push!(points_reference, (; id_file, freq_mhz, time_config_s, time_peak_s, time_peak_raw_s))
            val_prominence_v = ch2_v[idx_peak] - val_baseline_v
            println("Calibration '$tag_calibration' $name_file: config aligned $(time_config_s) s, $(freq_mhz) MHz → detected at aligned $(round(time_peak_s; digits=6)) s (raw $(round(time_peak_raw_s; digits=6)) s; CH1 zero $(round(time_alignment_zero_s; digits=6)) s); CH2=$(round(ch2_v[idx_peak]; sigdigits=6)) V, prominence over trimmed mean=$(round(val_prominence_v; sigdigits=5)) V")
        end
        GC.gc()
    end

    ids_file = unique([point.id_file for point in points_reference])
    if model == :cubicspline
        extrapolation = Symbol(lowercase(String(get(entry_calibration, "extrapolation", "linear"))))
        return plot_monotone_cubic_calibration(
            path_png, tag_calibration, points_reference, ids_file, extrapolation)
    end
    time_zero_by_file = Dict{Int,Float64}()
    ids_fitted_zero = Int[]
    for id_file in ids_file
        points_zero = filter(point -> point.id_file == id_file && iszero(point.freq_mhz), points_reference)
        length(points_zero) <= 1 || throw(ArgumentError("Expected at most one 0 MHz reference for SDS$(lpad(string(id_file), 5, '0')); found $(length(points_zero))"))
        if isempty(points_zero)
            ref_file = filter(point -> point.id_file == id_file, points_reference)
            time_zero_by_file[id_file] = mean(point.time_peak_s - point.freq_mhz / -472.0 for point in ref_file)
            push!(ids_fitted_zero, id_file)
        else
            time_zero_by_file[id_file] = only(points_zero).time_peak_s
        end
    end

    idx_group_by_file = Dict(id_file => idx for (idx, id_file) in enumerate(ids_file))
    idx_zero_by_file = Dict(id_file => idx for (idx, id_file) in enumerate(ids_fitted_zero))
    n_parameters = 3 + length(ids_file) + length(ids_fitted_zero)
    length(points_reference) > n_parameters ||
        throw(ArgumentError("Need more calibration references than the $n_parameters fitted parameters"))
    mat_x = vcat(
        reshape([point.time_peak_s for point in points_reference], 1, :),
        reshape([idx_group_by_file[point.id_file] for point in points_reference], 1, :),
    )
    val_frequency_mhz = [point.freq_mhz for point in points_reference]
    val_p_initial_all = vcat(0.08, -0.0013, 0.0, fill(-472.0, length(ids_file)),
        [time_zero_by_file[id_file] for id_file in ids_fitted_zero])
    function model_calibration(x, p)
        return [begin
            idx_group = Int(x[2, idx])
            id_file = ids_file[idx_group]
            t_zero_s = haskey(idx_zero_by_file, id_file) ?
                p[3 + length(ids_file) + idx_zero_by_file[id_file]] : time_zero_by_file[id_file]
            alpha_mhz_per_s = p[3 + idx_group]
            calc_frequency_from_time(alpha_mhz_per_s, p[1:3], x[1, idx] - t_zero_s)
        end for idx in axes(x, 2)]
    end
    function jacobian_calibration(x, p)
        jac = zeros(Float64, size(x, 2), length(p))
        for idx in axes(x, 2)
            idx_group = Int(x[2, idx])
            id_file = ids_file[idx_group]
            t_zero_s = haskey(idx_zero_by_file, id_file) ?
                p[3 + length(ids_file) + idx_zero_by_file[id_file]] : time_zero_by_file[id_file]
            u_s = x[1, idx] - t_zero_s
            alpha_mhz_per_s = p[3 + idx_group]
            shape = 1 + p[1] * u_s + p[2] * u_s^3 + p[3] * u_s^4
            jac[idx, 1] = alpha_mhz_per_s * u_s^2
            jac[idx, 2] = alpha_mhz_per_s * u_s^4
            jac[idx, 3] = alpha_mhz_per_s * u_s^5
            jac[idx, 3 + idx_group] = u_s * shape
            if haskey(idx_zero_by_file, id_file)
                jac[idx, 3 + length(ids_file) + idx_zero_by_file[id_file]] =
                    -alpha_mhz_per_s * (1 + 2p[1] * u_s + 4p[2] * u_s^3 + 5p[3] * u_s^4)
            end
        end
        return jac
    end
    fit_calibration = LsqFit.curve_fit(
        model_calibration, jacobian_calibration, mat_x, val_frequency_mhz, val_p_initial_all)
    fit_calibration.converged || @warn "Calibration fit did not report convergence"
    val_parameters = fit_calibration.param
    val_p_shape = val_parameters[1:3]
    alpha_by_file = Dict(id_file => val_parameters[3 + idx_group_by_file[id_file]] for id_file in ids_file)
    for id_file in ids_fitted_zero
        time_zero_by_file[id_file] = val_parameters[3 + length(ids_file) + idx_zero_by_file[id_file]]
    end

    function color_for_file(id_file::Int)
        idx_file = idx_group_by_file[id_file]
        return oklch_to_rgb(mod(24 + (idx_file - 1) * 360 / length(ids_file), 360), 0.60, 0.15)
    end
    marker_by_file = Dict(14 => :utriangle, 19 => :circle)

    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1];
        width=1200,
        height=850,
        title="Spectrum 626 time–frequency calibration: $tag_calibration",
        xlabel="Time from zero-frequency t₀ (s)",
        ylabel="Reference frequency (MHz)",
        titlesize=32,
        xlabelsize=28,
        ylabelsize=28,
        xticklabelsize=22,
        yticklabelsize=22,
    )

    for (idx_group, id_file) in enumerate(ids_file)
        ref_file = filter(point -> point.id_file == id_file, points_reference)
        t_zero_s = time_zero_by_file[id_file]
        color_file = color_for_file(id_file)
        marker = get(marker_by_file, id_file, :circle)
        val_u_file = [point.time_peak_s - t_zero_s for point in ref_file]
        alpha_mhz_per_s = alpha_by_file[id_file]
        u_min_ref, u_max_ref = extrema(val_u_file)
        u_min_model, u_max_model = find_monotonic_model_limits(
            alpha_mhz_per_s, val_p_shape, bounds_spectrum_frequency_mhz...)
        for (u_start, u_stop) in ((u_min_model, u_min_ref), (u_max_ref, u_max_model))
            u_start < u_stop || continue
            u_extrap = range(u_start, u_stop; length=400)
            lines!(ax, u_extrap,
                [calc_frequency_from_time(alpha_mhz_per_s, val_p_shape, u) for u in u_extrap];
                color=color_file, linewidth=2, linestyle=:dash)
        end
        u_model = range(u_min_ref, u_max_ref; length=700)
        freq_model = [calc_frequency_from_time(alpha_mhz_per_s, val_p_shape, u) for u in u_model]
        lines!(ax, u_model, freq_model;
            color=color_file, linewidth=2, label="SDS$(lpad(string(id_file), 5, '0')) model")
        offset_label_x_s = (idx_group - (length(ids_file) + 1) / 2) * 0.14
        offset_label_y_mhz = (idx_group - (length(ids_file) + 1) / 2) * 70
        scatter!(ax,
            [point.time_peak_s - t_zero_s for point in ref_file],
            [point.freq_mhz for point in ref_file];
            color=color_file, marker, markersize=24,
            strokecolor=:white, strokewidth=0.4,
        )
        for point in ref_file
            time_relative_s = point.time_peak_s - t_zero_s
            freq_model_mhz = calc_frequency_from_time(alpha_by_file[id_file], val_p_shape, time_relative_s)
            freq_error_mhz = point.freq_mhz - freq_model_mhz
            str_sign = freq_error_mhz >= 0 ? "+" : ""
            str_error = "Δf=$(str_sign)$(round(freq_error_mhz; digits=1)) MHz"
            text!(ax, time_relative_s + offset_label_x_s, point.freq_mhz + offset_label_y_mhz;
                text=str_error, color=color_file, fontsize=15,
                align=(:center, :center), strokecolor=:white, strokewidth=0.4)
            println("SDS$(lpad(string(id_file), 5, '0')) $(point.freq_mhz) MHz: Δf = $(round(freq_error_mhz; digits=3)) MHz")
        end
    end
    gl_text = GridLayout()
    fig[1, 2] = gl_text
    str_coefficients = "Shared shape coefficients\n" * join(
        ["p$(idx + 1) = $(round(val_p_shape[idx]; sigdigits=7))" for idx in eachindex(val_p_shape)], "\n")
    Label(gl_text[1, 1], str_coefficients;
        color=:black, fontsize=22, halign=:left, valign=:top)
    Label(gl_text[2, 1], "Zero-frequency t₀ (detected or fitted) and fitted scale";
        color=:black, fontsize=22, halign=:left, valign=:top)
    for (idx_file, id_file) in enumerate(ids_file)
        Label(gl_text[idx_file + 2, 1], "t₀[$id_file] = $(round(time_zero_by_file[id_file]; digits=6)) s; α = $(round(alpha_by_file[id_file]; digits=3)) MHz/s";
            color=color_for_file(id_file), fontsize=22, halign=:left, valign=:top)
    end
    Label(gl_text[length(ids_file) + 3, 1], "Solid: reference span; dashed: extrapolation";
        color=:black, fontsize=18, halign=:left, valign=:top)
    Legend(gl_text[length(ids_file) + 4, 1], ax;
        framevisible=false, labelsize=20, patchsize=(36, 24))
    resize_to_layout!(fig)
    mkpath(dirname(abspath(path_png)))
    save(path_png, fig; px_per_unit=2)
    println("Calibration '$tag_calibration' shared shape coefficients p2–p4: $(collect(val_p_shape))")
    println("Calibration '$tag_calibration' t0 by index: $time_zero_by_file")
    println("Fitted alpha by index: $alpha_by_file")
    println("Joint-fit RSS: $(LsqFit.rss(fit_calibration)); converged=$(fit_calibration.converged); parameters=$(length(val_parameters))")
    println("Saved resized calibration PNG to $path_png")
    u_calibration_range_s = extrema(
        point.time_peak_s - time_zero_by_file[point.id_file] for point in points_reference)
    return (; tag=tag_calibration, model=:p4, fig, val_p_shape, points_reference, ids_file,
        time_zero_by_file, alpha_by_file, u_calibration_range_s)
end

function draw_fine_scan_panels!(gl_fine_scan::GridLayout, traces_fine_scan::AbstractVector)
    for (idx_trace, trace) in enumerate(traces_fine_scan)
        pos_row = fld(idx_trace - 1, 3) + 1
        pos_col = mod(idx_trace - 1, 3) + 1
        ax_fine = Axis(gl_fine_scan[pos_row, pos_col];
            width=450, height=260,
            title="$(trace.tag): $(round(trace.target_frequency_mhz; digits=1)) MHz",
            xlabel="Frequency (MHz)", ylabel="Signal (V)",
            titlesize=15, xlabelsize=12, ylabelsize=12,
            xticklabelsize=11, yticklabelsize=11)
        if !isempty(trace.frequency_coarse_fine)
            lines!(ax_fine, trace.frequency_coarse_fine, trace.signal_coarse_fine_v;
                color=(trace.color_line, 0.35), linewidth=1.2, linestyle=:dash)
        end
        lines!(ax_fine, trace.frequency_fine_mhz, trace.signal_fine_v;
            color=(trace.color_line, 0.5), linewidth=2)
        xlims!(ax_fine, extrema(trace.frequency_fine_mhz)...)
    end
    return nothing
end

function plot_calibrated_spectrum_results(
    path_data::AbstractString,
    path_png::AbstractString,
    calibrations_by_tag::AbstractDict;
    normalize_amplitude::Bool=false,
)
    config = YAML.load_file(joinpath(path_data, "config.yaml"))
    entries_result = get(config, "result", nothing)
    isnothing(entries_result) && throw(ArgumentError("No 'result' section in config.yaml"))

    freq_min_mhz, freq_max_mhz = bounds_result_frequency_mhz
    time_uncertainty_s = 0.001
    axis_width_px = 1440 # 54 / 25.4 * 96
    axis_height_px = 960 # 34 / 25.4 * 96
    fig = Figure(fontsize=14)
    ax = Axis(fig[1, 1];
        width=axis_width_px,
        height=axis_height_px,
        xlabel="Frequency (MHz)",
        ylabel=normalize_amplitude ? "Normalized spectrum signal" : "Spectrum signal (V)",
        xlabelsize=10.667,
        ylabelsize=10.667,
        xticklabelsize=9.333,
        yticklabelsize=9.333,
        xticks=(vmin, vmax) -> begin
            if vmax - vmin >= 1000
                return (ceil(Int, vmin / 200) * 200):200:(floor(Int, vmax / 200) * 200)
            end
            return GLMakie.Makie.get_tickvalues(WilkinsonTicks(7), vmin, vmax)
        end,
        limits=((freq_min_mhz, freq_max_mhz), nothing),
    )

    traces_result = NamedTuple[]
    traces_fine_scan = NamedTuple[]
    mat_data_by_path = Dict{String,Any}()
    points_annotations = NamedTuple[]
    for (idx_entry, entry) in enumerate(entries_result)
        tag_calibration = String(entry["calibration"])
        haskey(calibrations_by_tag, tag_calibration) ||
            throw(ArgumentError("Result '$(entry["tag"])' selects unknown calibration '$tag_calibration'"))
        calibration = calibrations_by_tag[tag_calibration]
        is_spline = calibration.model == :cubicspline
        val_p_shape = is_spline ? nothing : calibration.val_p_shape
        alpha_by_file = is_spline ? nothing : calibration.alpha_by_file
        id_file = Int(entry["index"])
        name_file = "SDS$(lpad(string(id_file), 5, '0'))"
        tag = String(entry["tag"])
        path_jld2 = joinpath(path_data, "$name_file.jld2")
        isfile(path_jld2) || throw(ArgumentError("JLD2 cache does not exist: $path_jld2"))
        payload = JLD2.load(path_jld2)
        ch1_v = smooth_scope_channel(payload["ch1_v"])
        ch2_v = smooth_scope_channel(payload["ch2_v"])
        time_s = calc_time_axis(payload["metadata_lines"], length(ch2_v))
        time_alignment_zero_s = find_ch1_alignment_zero_crossing(ch1_v, time_s)
        time_aligned_s = time_s .- time_alignment_zero_s
        points_reference = NamedTuple[]
        for peak_config in entry["peaks"]
            get(peak_config, "ref", false) === true || continue
            freq_mhz = Float64(peak_config["freq"])
            haskey(peak_config, "time") && !is_missing_config_value(peak_config["time"]) ||
                throw(ArgumentError("Reference peak for $name_file must have a time: $peak_config"))
            time_config_s = Float64(peak_config["time"])
            idx_peak = find_local_peak_near_time(ch2_v, time_aligned_s, time_config_s, time_uncertainty_s)
            push!(points_reference, (; freq_mhz, time_config_s, idx_peak,
                time_peak_s=Float64(time_aligned_s[idx_peak]),
                time_peak_raw_s=Float64(time_s[idx_peak])))
        end
        isempty(points_reference) && throw(ArgumentError("No ref: true frequency references configured for $name_file"))
        points_zero = filter(point -> iszero(point.freq_mhz), points_reference)
        length(points_zero) <= 1 ||
            throw(ArgumentError("Expected at most one 0 MHz result reference for $name_file; found $(length(points_zero))"))
        if is_spline
            spline = if haskey(calibration.splines_by_file, id_file)
                calibration.splines_by_file[id_file]
            else
                fit_monotone_cubic([point.time_peak_s for point in points_reference],
                    [point.freq_mhz for point in points_reference];
                    extrapolation=calibration.extrapolation)
            end
            time_zero_s = find_monotone_cubic_time(spline, 0.0)
            println("$name_file ($tag): using monotone cubic calibration '$tag_calibration' with t₀=$(round(time_zero_s; digits=6)) s")
        elseif haskey(alpha_by_file, id_file)
            time_zero_s = calibration.time_zero_by_file[id_file]
            alpha_mhz_per_s = Float64(alpha_by_file[id_file])
            println("$name_file ($tag): using t₀=$(round(time_zero_s; digits=6)) s and α=$(round(alpha_mhz_per_s; digits=3)) MHz/s from calibration '$tag_calibration'")
        elseif isempty(points_zero)
            time_center_guess_s = mean(
                point.time_peak_s - point.freq_mhz / -472.0 for point in points_reference)
            fit_result = fit_time_zero_and_scale(
                time_center_guess_s,
                [point.time_peak_s for point in points_reference],
                [point.freq_mhz for point in points_reference],
                val_p_shape)
            time_zero_s = fit_result.time_zero_s
            alpha_mhz_per_s = fit_result.alpha_mhz_per_s
            println("$name_file ($tag): fitted t₀=$(round(time_zero_s; digits=6)) s and α=$(round(alpha_mhz_per_s; digits=3)) MHz/s from $(length(points_reference)) ref: true peaks")
        else
            time_zero_s = only(points_zero).time_peak_s
            fit_result = fit_scale_with_fixed_time_zero(
                time_zero_s,
                [point.time_peak_s for point in points_reference],
                [point.freq_mhz for point in points_reference],
                val_p_shape)
            alpha_mhz_per_s = fit_result.alpha_mhz_per_s
            println("$name_file ($tag): fixed t₀ to detected aligned 0 MHz peak at $(round(time_zero_s; digits=6)) s; fitted α=$(round(alpha_mhz_per_s; digits=3)) MHz/s")
        end
        for point in points_reference
            freq_calibrated_mhz = is_spline ?
                calc_monotone_cubic_frequency(spline, point.time_peak_s) :
                calc_frequency_from_time(alpha_mhz_per_s, val_p_shape, point.time_peak_s - time_zero_s)
            println("$name_file ($tag): reference $(point.freq_mhz) MHz at aligned t=" *
                "$(round(point.time_peak_s; digits=6)) s (raw $(round(point.time_peak_raw_s; digits=6)) s; " *
                "CH2=$(round(ch2_v[point.idx_peak]; sigdigits=6)) V) → " *
                "$(round(freq_calibrated_mhz; digits=3)) MHz " *
                "(residual $(round(point.freq_mhz - freq_calibrated_mhz; digits=3)) MHz)")
        end
        time_relative_s = time_aligned_s .- time_zero_s
        if is_spline
            frequency_mhz = [calc_monotone_cubic_frequency(spline, t) for t in time_aligned_s]
            time_min_plot_s, time_max_plot_s = find_monotone_cubic_model_limits(
                spline, freq_min_mhz, freq_max_mhz)
            idx_first = searchsortedfirst(time_aligned_s, time_min_plot_s)
            idx_last = searchsortedlast(time_aligned_s, time_max_plot_s)
            idx_first <= idx_last ||
                throw(ArgumentError("Spline result for $name_file does not intersect the plotted frequency range"))
            ids_trace = idx_first:idx_last
        else
            frequency_mhz = alpha_mhz_per_s .* time_relative_s .* (
                1 .+ val_p_shape[1] .* time_relative_s .+ val_p_shape[2] .* time_relative_s .^ 3 .+
                val_p_shape[3] .* time_relative_s .^ 4)
            ids_trace = find_monotonic_frequency_indices(
                time_relative_s, frequency_mhz, val_p_shape, freq_min_mhz, freq_max_mhz)
        end
        frequency_trace = frequency_mhz[ids_trace]
        signal_trace_v = ch2_v[ids_trace]
        val_ch2_baseline_v = normalize_amplitude ?
            calc_mean_excluding_top_fraction(signal_trace_v, 0.01) : 0.0
        color_line = parse(Colorant, String(entry["color"]))
        push!(traces_result, (; name_file, tag, frequency_trace, signal_trace_v,
            val_ch2_baseline_v, color_line))
        for peak_config in entry["peaks"]
            replace_config = get(peak_config, "replace", nothing)
            if !normalize_amplitude && !is_missing_config_value(replace_config)
                replace_config isa AbstractVector && length(replace_config) == 2 &&
                    all(value -> value isa AbstractString && !isempty(strip(value)), replace_config) ||
                    throw(ArgumentError("$name_file peak replace must be [MAT filename, group_data field]; got $replace_config"))
                name_mat, name_case = String.(replace_config)
                path_mat = joinpath(path_data, "individual reprocess", name_mat)
                isfile(path_mat) || throw(ArgumentError("Fine-scan MAT file does not exist: $path_mat"))
                mat_data = get!(mat_data_by_path, path_mat) do
                    MAT.matread(path_mat)
                end
                haskey(mat_data, "group_data") ||
                    throw(ArgumentError("$path_mat has no group_data structure"))
                group_data = mat_data["group_data"]
                haskey(group_data, name_case) ||
                    throw(ArgumentError("$path_mat group_data has no case '$name_case'"))
                case_data = group_data[name_case]
                all(haskey(case_data, key) for key in ("frequency", "voltage")) ||
                    throw(ArgumentError("$path_mat group_data.$name_case needs frequency and voltage fields"))
                frequency_fine_mhz = vec(Float64.(case_data["frequency"]))
                signal_fine_v = smooth_scope_channel(vec(Float64.(case_data["voltage"])))
                length(frequency_fine_mhz) == length(signal_fine_v) && length(frequency_fine_mhz) >= 2 ||
                    throw(DimensionMismatch("$path_mat group_data.$name_case has $(length(frequency_fine_mhz)) frequency and $(length(signal_fine_v)) voltage values; expected matching lengths of at least 2"))
                all(isfinite, frequency_fine_mhz) && all(isfinite, signal_fine_v) ||
                    throw(ArgumentError("$path_mat group_data.$name_case contains nonfinite values"))
                target_frequency_mhz = Float64(get(peak_config, "freq", mean(frequency_fine_mhz)))
                freq_min_fine, freq_max_fine = extrema(frequency_fine_mhz)
                ids_coarse_fine = findall(freq -> freq_min_fine <= freq <= freq_max_fine,
                    frequency_trace)
                frequency_coarse_fine = frequency_trace[ids_coarse_fine]
                signal_coarse_fine_v = signal_trace_v[ids_coarse_fine]
                push!(traces_fine_scan, (; frequency_fine_mhz, signal_fine_v,
                    frequency_coarse_fine, signal_coarse_fine_v,
                    color_line, tag, target_frequency_mhz))
                println("$name_file ($tag): overlaying $name_mat group_data.$name_case over $(round(extrema(frequency_fine_mhz)[1]; digits=3))…$(round(extrema(frequency_fine_mhz)[2]; digits=3)) MHz")
            end
            is_reference = get(peak_config, "ref", false) === true
            has_time = haskey(peak_config, "time") && !is_missing_config_value(peak_config["time"])
            if is_reference
                has_time || throw(ArgumentError("Reference peak for $name_file must have a time: $peak_config"))
                time_config_s = Float64(peak_config["time"])
                idx_peak = find_local_peak_near_time(ch2_v, time_aligned_s, time_config_s, time_uncertainty_s)
            elseif haskey(peak_config, "freq")
                target_frequency_mhz = Float64(peak_config["freq"])
                idx_peak = first(ids_trace) - 1 + find_local_peak_near_frequency(
                    @view(ch2_v[ids_trace]),
                    @view(frequency_mhz[ids_trace]),
                    target_frequency_mhz, 5.0)
                println("$name_file ($tag): frequency-only peak target $(target_frequency_mhz) MHz → selected sample at $(round(frequency_mhz[idx_peak]; digits=3)) MHz, $(round(ch2_v[idx_peak]; sigdigits=6)) V")
            elseif has_time
                time_config_s = Float64(peak_config["time"])
                idx_peak = find_local_peak_near_time(ch2_v, time_aligned_s, time_config_s, time_uncertainty_s)
            else
                throw(ArgumentError("Peak needs a frequency or time: $peak_config"))
            end
            time_peak_s = Float64(time_aligned_s[idx_peak])
            freq_peak_mhz = Float64(frequency_mhz[idx_peak])
            freq_marker_mhz = clamp(freq_peak_mhz, freq_min_mhz, freq_max_mhz)
            val_peak_v = ch2_v[idx_peak]
            peak_tag = String(get(peak_config, "tag", ""))
            offset_label_x_mhz = if freq_peak_mhz > 620
                -18.0
            elseif freq_peak_mhz < -620
                18.0
            elseif is_reference && iszero(Float64(peak_config["freq"]))
                val_peak_v > 0.3 ? 12.0 : -12.0
            elseif id_file == 21
                freq_peak_mhz > 0 ? -12.0 : 12.0
            else
                freq_peak_mhz > 0 ? 12.0 : -12.0
            end
            offset_label_y_v = if freq_peak_mhz < -750
                0.28 - val_peak_v
            elseif freq_peak_mhz < -600
                (id_file == 21 ? 0.12 : 0.20) - val_peak_v
            elseif freq_peak_mhz > 750
                0.20 - val_peak_v
            elseif freq_peak_mhz > 600
                0.12 - val_peak_v
            elseif val_peak_v > 0.5
                -0.012 - 0.055 * (idx_entry - 1)
            else
                0.012 + 0.055 * (idx_entry - 1)
            end
            push!(points_annotations, (; freq_marker_mhz, val_peak_v, val_ch2_baseline_v,
                time_peak_s, freq_peak_mhz, peak_tag, color_line,
                marker=is_reference ? :utriangle : :circle,
                text_x=freq_marker_mhz + offset_label_x_mhz,
                offset_label_y_v,
                text_align=(offset_label_x_mhz > 0 ? :left : :right,
                    offset_label_y_v > 0 ? :bottom : :top)))
            println("$name_file ($tag): $(is_reference ? "reference" : "sample") peak at aligned $(round(time_peak_s; digits=6)) s, $(round(freq_peak_mhz; digits=3)) MHz, $(round(val_peak_v; sigdigits=6)) V$(isempty(peak_tag) ? "" : " — $peak_tag")")
        end
        println("$name_file ($tag): using calibration '$tag_calibration' with plotted frequency range $(round(extrema(frequency_trace)[1]; digits=3))…$(round(extrema(frequency_trace)[2]; digits=3)) MHz; " *
            "CH1 alignment zero=$(round(time_alignment_zero_s; digits=6)) s; raw CH2 voltage range=$(round(minimum(signal_trace_v); sigdigits=6))…$(round(maximum(signal_trace_v); sigdigits=6)) V")
        GC.gc()
    end

    val_normalization_v = 1.0
    if normalize_amplitude
        val_peak_by_result_v = [maximum(trace.signal_trace_v) - trace.val_ch2_baseline_v
            for trace in traces_result]
        all(isfinite, val_peak_by_result_v) ||
            throw(ArgumentError("Result traces contain a nonfinite baseline-corrected maximum"))
        pos_normalization = argmax(val_peak_by_result_v)
        val_normalization_v = val_peak_by_result_v[pos_normalization]
        val_normalization_v > 0 ||
            throw(ArgumentError("The largest baseline-corrected result signal must be positive; got $val_normalization_v V"))
        println("Normalizing result spectra by $(traces_result[pos_normalization].name_file) ($(traces_result[pos_normalization].tag)) maximum $(round(val_normalization_v; sigdigits=6)) V above its baseline")
    else
        println("Plotting result spectra in raw CH2 volts")
    end
    for trace in traces_result
        signal_plot = normalize_amplitude ?
            (trace.signal_trace_v .- trace.val_ch2_baseline_v) ./ val_normalization_v : trace.signal_trace_v
        lines!(ax, trace.frequency_trace, signal_plot;
            color=trace.color_line, linewidth=1.2, label=trace.tag)
    end
    for trace in traces_fine_scan
        lines!(ax, trace.frequency_fine_mhz, trace.signal_fine_v;
            color=(trace.color_line, 0.5), linewidth=1.2)
    end
    for point in points_annotations
        val_peak_plot = normalize_amplitude ?
            (point.val_peak_v - point.val_ch2_baseline_v) / val_normalization_v : point.val_peak_v
        offset_label_y_plot = normalize_amplitude ?
            point.offset_label_y_v / val_normalization_v : point.offset_label_y_v
        text_y = val_peak_plot + offset_label_y_plot
        text_value = normalize_amplitude ?
            "A=$(round(val_peak_plot; digits=3))" : "V=$(round(val_peak_plot; digits=3)) V"
        text_label = "t=$(round(point.time_peak_s; digits=3)) s\nf=$(round(point.freq_peak_mhz; digits=1)) MHz\n$text_value" *
            (isempty(point.peak_tag) ? "" : "\n$(point.peak_tag)")
        scatter!(ax, [point.freq_marker_mhz], [val_peak_plot];
            color=point.color_line, marker=point.marker, markersize=8,
            strokecolor=:black, strokewidth=0.8)
        text!(ax, point.text_x, text_y;
            text=text_label, color=:black, fontsize=14,
            align=point.text_align, strokecolor=:white, strokewidth=0.4)
    end

    axislegend(ax; position=:rt, framevisible=false, labelsize=9.333)
    if !isempty(traces_fine_scan)
        Label(fig[2, 1], "Fine-scan detail · dashed: coarse scan · solid: fine scan (50% opacity)";
            halign=:left, fontsize=14)
        gl_fine_scan = GridLayout()
        fig[3, 1] = gl_fine_scan
        draw_fine_scan_panels!(gl_fine_scan, traces_fine_scan)
    end
    resize_to_layout!(fig)
    mkpath(dirname(abspath(path_png)))
    save(path_png, fig; px_per_unit=4)
    println("Saved result spectrum PNG with $(round(axis_width_px; digits=1))×$(round(axis_height_px; digits=1)) px axis to $path_png")
    if !isempty(traces_fine_scan)
        fig_detail = Figure(fontsize=14)
        Label(fig_detail[1, 1], "Fine-scan overlays · dashed: coarse scan · solid: fine scan (50% opacity)";
            halign=:left, fontsize=14)
        gl_detail = GridLayout()
        fig_detail[2, 1] = gl_detail
        draw_fine_scan_panels!(gl_detail, traces_fine_scan)
        resize_to_layout!(fig_detail)
        path_detail_png = splitext(path_png)[1] * "_fine_scan_detail.png"
        save(path_detail_png, fig_detail; px_per_unit=3)
        println("Saved fine-scan detail PNG to $path_detail_png")
    end
    return fig
end

fig = function plot_ch1_spectrum(
    path_jld2::AbstractString,
    path_png::AbstractString,
    path_period_csv::AbstractString,
)
    payload = JLD2.load(path_jld2)
    metadata_lines = payload["metadata_lines"]
    ch1_v = smooth_scope_channel(payload["ch1_v"])
    ch2_v = smooth_scope_channel(payload["ch2_v"])
    time_key = haskey(payload, "time_s") ? "time_s" : "time_us"
    time_data = payload[time_key]
    length(time_data) == length(ch1_v) ||
        throw(DimensionMismatch("JLD2 time and CH1 arrays must have equal lengths; got $(length(time_data)) and $(length(ch1_v))"))
    time_s = calc_time_axis(metadata_lines, length(ch1_v))
    units = parse_scope_metadata_value(metadata_lines, "Horizontal Units")
    scale_per_div = parse(Float64, parse_scope_metadata_value(metadata_lines, "Horizontal Scale"))

    px_width = 3840
    px_height = 2160
    dead_time_s = 2.0
    idx_minima = find_windowed_extrema(ch1_v, time_s, true, dead_time_s)
    idx_maxima = find_windowed_extrema(ch1_v, time_s, false, dead_time_s)
    time_minima = time_s[idx_minima]
    time_maxima = time_s[idx_maxima]
    fig = Figure(size=(px_width, px_height), fontsize=28)
    ax = Axis(fig[1, 1];
        title="SDS00014 — Channel 1",
        xlabel="Time (s)",
        ylabel="CH1 voltage (V)",
        titlesize=36,
        xlabelsize=32,
        ylabelsize=32,
        xticklabelsize=26,
        yticklabelsize=26,
    )
    lines!(ax, time_s, ch1_v; color=:dodgerblue, linewidth=1.5)
    vlines!(ax, time_minima; color=:red, linewidth=2, linestyle=:dash, label="Minima ($(length(time_minima)))")
    vlines!(ax, time_maxima; color=:darkorange, linewidth=2, linestyle=:dash, label="Maxima ($(length(time_maxima)))")
    axislegend(ax; position=:rt, framevisible=false, labelsize=24)
    save(path_png, fig; px_per_unit=1)
    println("Applied Horizontal Scale: $scale_per_div $units/div across $n_horizontal_divisions divisions")
    println("Selected $(length(idx_minima)) minima and $(length(idx_maxima)) maxima with $(dead_time_s) s dead-time")
    println("Minimum times/values: $(collect(zip(round.(time_minima; digits=6), ch1_v[idx_minima])))")
    println("Maximum times/values: $(collect(zip(round.(time_maxima; digits=6), ch1_v[idx_maxima])))")
    println("Minimum separation: $(length(time_minima) > 1 ? minimum(diff(time_minima)) : NaN) s minima; $(length(time_maxima) > 1 ? minimum(diff(time_maxima)) : NaN) s maxima")
    println("Saved $(px_width)×$(px_height) PNG to $path_png")
    write_ch1_periods_csv(path_period_csv, time_s, ch1_v, idx_minima, idx_maxima)
    return nothing
end

args_data = filter(!=("--skip-conversion"), ARGS)
args_data = filter(!=("--calibration-only"), args_data)
args_data = filter(!=("--display-only"), args_data)
args_data = filter(!=("--result-voltage"), args_data)
args_data = filter(!=("--result-normalized"), args_data)
length(args_data) <= 1 || throw(ArgumentError(
    "Usage: julia --project=. multidualmot/spectral626.jl [data-folder] [--skip-conversion] [--calibration-only|--display-only] [--result-voltage|--result-normalized]"))
!("--calibration-only" in ARGS && "--display-only" in ARGS) ||
    throw(ArgumentError("--calibration-only and --display-only cannot be combined"))
!("--result-voltage" in ARGS && "--result-normalized" in ARGS) ||
    throw(ArgumentError("--result-voltage and --result-normalized cannot be combined"))
path_data = isempty(args_data) ? path_data_root : abspath(args_data[1])
path_config = joinpath(path_data, "config.yaml")
isfile(path_config) || throw(ArgumentError("Spectrum config does not exist: $path_config"))
config = YAML.load_file(path_config)
if "--skip-conversion" in ARGS
    println("Skipping CSV conversion; using existing JLD2 caches in $path_data")
else
    convert_missing_scope_csvs(path_data)
end

# The CH1 overview and period table remain specific to SDS00014.
if !("--calibration-only" in ARGS || "--display-only" in ARGS)
    path_jld2 = joinpath(path_data, "SDS00014.jld2")
    isfile(path_jld2) || throw(ArgumentError("SDS00014 JLD2 cache is missing after conversion: $path_jld2"))
    path_png = joinpath(path_data, "SDS00014_CH1_vs_time.png")
    path_period_csv = joinpath(path_data, "SDS00014_CH1_periods.csv")
    plot_ch1_spectrum(path_jld2, path_png, path_period_csv)
end
if !("--calibration-only" in ARGS)
    entries_display = get(config, "display", nothing)
    isnothing(entries_display) && throw(ArgumentError("No 'display' section in $path_config"))
    fig_spectra = []
    for entry in entries_display
        ids_display = Int.(entry["indices"])
        isempty(ids_display) && throw(ArgumentError("A display entry has no indices"))
        name_display = join([lpad(string(id), idx == 1 ? 5 : 2, '0')
            for (idx, id) in enumerate(ids_display)], "_")
        path_ch2_png = joinpath(path_data, "SDS$(name_display)_CH2_aligned_periods.png")
        push!(fig_spectra, plot_ch2_rising_ramps(path_data, ids_display, path_ch2_png; val_min_peak_prominence_v=0.003))
    end
    length(fig_spectra) == 1 && (fig_spectrum = only(fig_spectra))
end

if !("--display-only" in ARGS)
entries_calibration = get(config, "calibration", nothing)
isnothing(entries_calibration) && throw(ArgumentError("No 'calibration' section in $path_config"))
calibration_results = []
calibrations_by_tag = Dict{String,Any}()
names_calibration_png = Set{String}()
for entry in entries_calibration
    tag_calibration = String(entry["tag"])
    isempty(strip(tag_calibration)) && throw(ArgumentError("Calibration tag cannot be empty"))
    haskey(calibrations_by_tag, tag_calibration) &&
        throw(ArgumentError("Duplicate calibration tag '$tag_calibration'"))
    name_tag_png = replace(tag_calibration, r"[^A-Za-z0-9_-]+" => "_")
    name_tag_png in names_calibration_png &&
        throw(ArgumentError("Calibration tags must produce unique output filenames; '$tag_calibration' conflicts with another tag"))
    push!(names_calibration_png, name_tag_png)
    path_calibration_png = joinpath(path_data, "spectrum626_time_frequency_calibration_$(name_tag_png).png")
    calibration_result = plot_spectrum_frequency_calibration(path_data, path_calibration_png, entry)
    push!(calibration_results, calibration_result)
    calibrations_by_tag[tag_calibration] = calibration_result
end
fig_calibrations = [calibration.fig for calibration in calibration_results]
path_result_spectrum_png = joinpath(path_data, "spectrum626_calibrated_results.png")
fig_result_spectrum = plot_calibrated_spectrum_results(
    path_data, path_result_spectrum_png, calibrations_by_tag;
    normalize_amplitude=("--result-normalized" in ARGS))
end
