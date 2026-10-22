using CSV
using GLMakie
using Colors
using DataFrames
using JLD2
using LsqFit
using Peaks
using YAML

GLMakie.activate!()
const path_data_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\spectrum 626"
const n_horizontal_divisions = 14

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

    n_bins = ceil(Int, (time_s[end] - time_s[1]) / dead_time_s)
    val_bin_mean = Float64[]
    idx_bin_first = Int[]
    idx_bin_last = Int[]
    sizehint!(val_bin_mean, n_bins)
    sizehint!(idx_bin_first, n_bins)
    sizehint!(idx_bin_last, n_bins)

    for idx_bin in 1:n_bins
        time_first = time_s[1] + (idx_bin - 1) * dead_time_s
        time_last = min(time_first + dead_time_s, time_s[end])
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
    idx_candidates = findmaxima(y).indices
    idx_selected = select_deadtime_extrema(idx_candidates, y, time_s, false, dead_time_s)
    filter!(idx_peak -> y[idx_peak] - baseline_v > min_prominence_v, idx_selected)
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

    isempty(period_rows) && throw(ArgumentError("At least two minima are required to calculate a completed period"))
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
)::Float64
    idx_start < idx_end || throw(ArgumentError("A period must contain at least two samples"))
    for idx in idx_start:(idx_end - 1)
        val = y[idx]
        val_next = y[idx + 1]
        iszero(val) && return Float64(time_s[idx])
        if val < 0 && val_next >= 0
            frac = -val / (val_next - val)
            return Float64(time_s[idx] + frac * (time_s[idx + 1] - time_s[idx]))
        end
    end
    throw(ArgumentError("No rising CH1 zero-crossing found in period samples $idx_start:$idx_end"))
end

function find_ch1_alignment_zero_crossing(
    ch1_v::AbstractVector{<:Real},
    time_s::AbstractVector{<:Real},
    dead_time_s::Real=2.0,
)::Float64
    length(ch1_v) == length(time_s) || throw(DimensionMismatch("CH1 and time arrays must have equal lengths"))
    idx_minima = find_windowed_extrema(ch1_v, time_s, true, dead_time_s)
    length(idx_minima) >= 2 || throw(ArgumentError("At least two CH1 minima are required to align one period"))
    idx_period = findfirst(1:(length(idx_minima) - 1)) do idx
        time_s[idx_minima[idx]] <= 0 <= time_s[idx_minima[idx + 1]]
    end
    isnothing(idx_period) && throw(ArgumentError("No completed CH1 period spans raw t = 0"))
    return find_first_uphill_zero_crossing(
        ch1_v, time_s, idx_minima[idx_period], idx_minima[idx_period + 1])
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
    ch1_v = payload["ch1_v"]
    ch2_v = payload["ch2_v"]
    metadata_lines = payload["metadata_lines"]
    length(ch1_v) == length(ch2_v) ||
        throw(DimensionMismatch("CH1 and CH2 arrays must have equal lengths in $path_jld2"))
    time_s = calc_time_axis(metadata_lines, length(ch1_v))
    idx_minima = find_windowed_extrema(ch1_v, time_s, true, dead_time_s)
    idx_maxima = find_windowed_extrema(ch1_v, time_s, false, dead_time_s)
    return (; time_s, ch1_v, ch2_v, idx_minima, idx_maxima)
end

function plot_ch2_periods_spanning_zero(
    path_data::AbstractString,
    ids_file::AbstractVector{<:Integer},
    path_png::AbstractString,
)
    isempty(ids_file) && throw(ArgumentError("At least one spectrum file is required"))
    length(ids_file) <= 6 || throw(ArgumentError("The existing OKLCH palette supports at most six files"))

    val_offset_per_file_v = 0.5
    val_min_peak_prominence_v = 0.008
    val_fraction_top_excluded = 0.01
    val_peak_dead_time_s = 0.05
    val_hue_deg_file = [24, 80, 136, 192, 248, 304]
    val_lightness_ch2 = 0.60
    val_lightness_ch1 = 0.88
    val_chroma = 0.15
    colors_ch2 = [oklch_to_rgb(val_hue_deg_file[idx], val_lightness_ch2, val_chroma) for idx in eachindex(ids_file)]
    colors_ch1 = [oklch_to_rgb(val_hue_deg_file[idx], val_lightness_ch1, val_chroma) for idx in eachindex(ids_file)]

    px_width = 3840
    px_height = 2160
    fig = Figure(size=(px_width, px_height), fontsize=28)
    ax = Axis(fig[1, 1];
        title="CH1 and CH2 cycles spanning t = 0, aligned to CH1 zero-crossings",
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

        pos_period = findfirst(1:(length(idx_minima) - 1)) do idx_period
            time_s[idx_minima[idx_period]] <= 0 <= time_s[idx_minima[idx_period + 1]]
        end
        isnothing(pos_period) && throw(ArgumentError("No completed CH1 period in $name_file spans t = 0"))

        idx_start = idx_minima[pos_period]
        idx_end = idx_minima[pos_period + 1]
        any(idx_start < idx_peak < idx_end for idx_peak in idx_maxima) ||
            throw(ArgumentError("No CH1 maximum found inside the selected period in $name_file"))
        time_zero = find_first_uphill_zero_crossing(ch1_v, time_s, idx_start, idx_end)
        val_offset_v = (idx_file - 1) * val_offset_per_file_v
        time_period_offset = @view(time_s[idx_start:idx_end]) .- time_zero
        ch1_period_offset = (@view ch1_v[idx_start:idx_end]) ./ 8 .+ val_offset_v
        time_period = @view time_s[idx_start:idx_end]
        ch2_period = @view ch2_v[idx_start:idx_end]
        ch2_period_offset = ch2_period .+ val_offset_v
        lines!(ax, time_period_offset, ch1_period_offset;
            color=colors_ch1[idx_file], linewidth=1.5)
        lines!(ax, time_period_offset, ch2_period_offset;
            color=colors_ch2[idx_file], linewidth=1.5, label=name_file)
        val_ch2_baseline_v = calc_mean_excluding_top_fraction(ch2_period, val_fraction_top_excluded)
        idx_peaks_ch2 = find_local_peaks_with_deadtime(
            ch2_period, time_period, val_peak_dead_time_s, val_ch2_baseline_v, val_min_peak_prominence_v)
        for idx_peak in idx_peaks_ch2
            time_peak = time_period[idx_peak] - time_zero
            value_peak = ch2_period[idx_peak] + val_offset_v
            scatter!(ax, [time_peak], [value_peak];
                color=colors_ch2[idx_file], markersize=10, strokecolor=:black, strokewidth=1)
            text!(ax, time_peak, value_peak + 0.025;
                text=string(round(time_peak; digits=3)),
                color=colors_ch2[idx_file], fontsize=16,
                rotation=π / 2, align=(:center, :bottom),
                strokecolor=:white, strokewidth=0.4)
        end
        println("$name_file: selected period $(round(time_s[idx_start]; digits=6)) to $(round(time_s[idx_end]; digits=6)) s; rising zero-crossing $(round(time_zero; digits=6)) s; CH2 offset $val_offset_v V; trimmed CH2 baseline $(round(val_ch2_baseline_v; digits=4)) V; annotated $(length(idx_peaks_ch2)) peaks (>$val_min_peak_prominence_v V above baseline, $(val_peak_dead_time_s) s dead-time)")
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

function calc_frequency_from_time(val_p::AbstractVector{<:Real}, u_s::Real)::Float64
    return val_p[1] * u_s + val_p[2] * u_s^2 + val_p[3] * u_s^3 + val_p[4] * u_s^4
end

function fit_time_zero_from_references(
    val_p::AbstractVector{<:Real},
    time_peak_s::AbstractVector{<:Real},
    freq_mhz::AbstractVector{<:Real},
    u_calibration_range_s::Tuple{<:Real,<:Real},
)
    length(time_peak_s) == length(freq_mhz) ||
        throw(DimensionMismatch("Reference times and frequencies must have equal lengths"))
    length(time_peak_s) >= 2 || throw(ArgumentError("At least two references are needed to fit t₀"))
    u_calibration_range_s[1] < u_calibration_range_s[2] ||
        throw(ArgumentError("Calibration time range must be increasing"))
    time_zero_lower_s = maximum(time_peak_s .- u_calibration_range_s[2])
    time_zero_upper_s = minimum(time_peak_s .- u_calibration_range_s[1])
    time_zero_lower_s < time_zero_upper_s ||
        throw(ArgumentError("No t₀ keeps all reference peaks inside the calibrated time range"))

    function calc_sum_squared_frequency_error(time_zero_s::Real)::Float64
        val_error_sum_mhz2 = 0.0
        for idx in eachindex(time_peak_s)
            freq_fit_mhz = calc_frequency_from_time(val_p, time_peak_s[idx] - time_zero_s)
            val_error_sum_mhz2 += (freq_fit_mhz - freq_mhz[idx])^2
        end
        return val_error_sum_mhz2
    end

    # Keep every reference within the time domain supported by the calibration.
    # Within that domain, explicitly minimize the per-entry sum of squared errors.
    val_t_zero_grid_s = range(time_zero_lower_s, time_zero_upper_s; length=1001)
    pos_best = argmin(calc_sum_squared_frequency_error.(val_t_zero_grid_s))
    pos_best == firstindex(val_t_zero_grid_s) && return first(val_t_zero_grid_s)
    pos_best == lastindex(val_t_zero_grid_s) && return last(val_t_zero_grid_s)
    time_left_s = val_t_zero_grid_s[pos_best - 1]
    time_right_s = val_t_zero_grid_s[pos_best + 1]
    val_golden = (sqrt(5) - 1) / 2
    time_x1_s = time_right_s - val_golden * (time_right_s - time_left_s)
    time_x2_s = time_left_s + val_golden * (time_right_s - time_left_s)
    val_error_x1 = calc_sum_squared_frequency_error(time_x1_s)
    val_error_x2 = calc_sum_squared_frequency_error(time_x2_s)
    for _ in 1:100
        if val_error_x1 <= val_error_x2
            time_right_s = time_x2_s
            time_x2_s = time_x1_s
            val_error_x2 = val_error_x1
            time_x1_s = time_right_s - val_golden * (time_right_s - time_left_s)
            val_error_x1 = calc_sum_squared_frequency_error(time_x1_s)
        else
            time_left_s = time_x1_s
            time_x1_s = time_x2_s
            val_error_x1 = val_error_x2
            time_x2_s = time_left_s + val_golden * (time_right_s - time_left_s)
            val_error_x2 = calc_sum_squared_frequency_error(time_x2_s)
        end
    end
    return (time_left_s + time_right_s) / 2
end

function plot_spectrum_frequency_calibration(
    path_data::AbstractString,
    path_png::AbstractString,
)
    path_config = joinpath(path_data, "config.yaml")
    isfile(path_config) || throw(ArgumentError("Calibration config does not exist: $path_config"))
    config = YAML.load_file(path_config)
    entries_calibration = get(config, "calibration", nothing)
    isnothing(entries_calibration) && throw(ArgumentError("No 'calibration' section in $path_config"))

    time_uncertainty_s = 0.001
    points_reference = NamedTuple[]
    for entry in entries_calibration
        id_file = Int(entry["index"])
        name_file = "SDS$(lpad(string(id_file), 5, '0'))"
        path_jld2 = joinpath(path_data, "$name_file.jld2")
        isfile(path_jld2) || throw(ArgumentError("JLD2 cache does not exist: $path_jld2"))
        payload = JLD2.load(path_jld2)
        ch1_v = payload["ch1_v"]
        ch2_v = payload["ch2_v"]
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
            println("$name_file: config aligned $(time_config_s) s, $(freq_mhz) MHz → detected at aligned $(round(time_peak_s; digits=6)) s (raw $(round(time_peak_raw_s; digits=6)) s; CH1 zero $(round(time_alignment_zero_s; digits=6)) s); CH2=$(round(ch2_v[idx_peak]; sigdigits=6)) V, prominence over trimmed mean=$(round(val_prominence_v; sigdigits=5)) V")
        end
        GC.gc()
    end

    ids_file = unique([point.id_file for point in points_reference])
    time_zero_by_file = Dict{Int,Float64}()
    for id_file in ids_file
        points_zero = filter(point -> point.id_file == id_file && iszero(point.freq_mhz), points_reference)
        length(points_zero) == 1 || throw(ArgumentError("Expected exactly one 0 MHz reference for SDS$(lpad(string(id_file), 5, '0')); found $(length(points_zero))"))
        time_zero_by_file[id_file] = only(points_zero).time_peak_s
    end

    length(points_reference) > 4 + length(ids_file) ||
        throw(ArgumentError("Need more calibration references than the $(4 + length(ids_file)) fitted parameters"))
    points_nonzero = filter(point -> !iszero(point.freq_mhz), points_reference)
    val_u_initial = [point.time_peak_s - time_zero_by_file[point.id_file] for point in points_nonzero]
    val_freq_nonzero = [point.freq_mhz for point in points_nonzero]
    mat_design_initial = hcat(val_u_initial, val_u_initial .^ 2, val_u_initial .^ 3, val_u_initial .^ 4)
    val_p_initial = mat_design_initial \ val_freq_nonzero
    idx_group_by_file = Dict(id_file => idx for (idx, id_file) in enumerate(ids_file))
    mat_x = vcat(
        reshape([point.time_peak_s for point in points_reference], 1, :),
        reshape([idx_group_by_file[point.id_file] for point in points_reference], 1, :),
    )
    val_frequency_mhz = [point.freq_mhz for point in points_reference]
    val_p_initial_all = vcat(val_p_initial, [time_zero_by_file[id_file] for id_file in ids_file])
    function model_calibration(x, p)
        return [begin
            u = x[1, idx] - p[4 + Int(x[2, idx])]
            p[1] * u + p[2] * u^2 + p[3] * u^3 + p[4] * u^4
        end for idx in axes(x, 2)]
    end
    function jacobian_calibration(x, p)
        jac = zeros(Float64, size(x, 2), length(p))
        for idx in axes(x, 2)
            u = x[1, idx] - p[4 + Int(x[2, idx])]
            jac[idx, 1:4] .= (u, u^2, u^3, u^4)
            idx_t0 = 4 + Int(x[2, idx])
            jac[idx, idx_t0] = -(p[1] + 2p[2] * u + 3p[3] * u^2 + 4p[4] * u^3)
        end
        return jac
    end
    fit_calibration = LsqFit.curve_fit(
        model_calibration, jacobian_calibration, mat_x, val_frequency_mhz, val_p_initial_all)
    fit_calibration.converged || @warn "Calibration fit did not report convergence"
    val_parameters = fit_calibration.param
    val_p = val_parameters[1:4]
    for id_file in ids_file
        time_zero_by_file[id_file] = val_parameters[4 + idx_group_by_file[id_file]]
    end

    ids_palette = [14, 15, 16, 19, 21]
    val_hue_deg_file = [24, 80, 136, 192, 248]
    function color_for_file(id_file::Int)
        idx_file = findfirst(==(id_file), ids_palette)
        isnothing(idx_file) && throw(ArgumentError("No matching spectrum hue configured for index $id_file"))
        return oklch_to_rgb(val_hue_deg_file[idx_file], 0.60, 0.15)
    end
    marker_by_file = Dict(14 => :utriangle, 19 => :circle)

    fig = Figure(fontsize=24)
    ax = Axis(fig[1, 1];
        width=1200,
        height=850,
        title="Spectrum 626 time–frequency calibration",
        xlabel="Time from fitted t₀ (s)",
        ylabel="Reference frequency (MHz)",
        titlesize=32,
        xlabelsize=28,
        ylabelsize=28,
        xticklabelsize=22,
        yticklabelsize=22,
    )

    val_u_all = [point.time_peak_s - time_zero_by_file[point.id_file] for point in points_reference]
    u_min, u_max = extrema(val_u_all)
    u_model = range(u_min, u_max; length=700)
    freq_model = val_p[1] .* u_model .+ val_p[2] .* u_model .^ 2 .+
        val_p[3] .* u_model .^ 3 .+ val_p[4] .* u_model .^ 4
    lines!(ax, u_model, freq_model;
        color=:black, linewidth=3, label="Shared quartic model")

    for (idx_group, id_file) in enumerate(ids_file)
        ref_file = filter(point -> point.id_file == id_file, points_reference)
        t_zero_s = time_zero_by_file[id_file]
        color_file = color_for_file(id_file)
        marker = get(marker_by_file, id_file, :circle)
        offset_label_x_s = (idx_group - (length(ids_file) + 1) / 2) * 0.14
        offset_label_y_mhz = (idx_group - (length(ids_file) + 1) / 2) * 70
        scatter!(ax,
            [point.time_peak_s - t_zero_s for point in ref_file],
            [point.freq_mhz for point in ref_file];
            color=color_file, marker, markersize=24,
            strokecolor=:white, strokewidth=0.4,
            label="SDS$(lpad(string(id_file), 5, '0'))")
        for point in ref_file
            time_relative_s = point.time_peak_s - t_zero_s
            freq_model_mhz = val_p[1] * time_relative_s + val_p[2] * time_relative_s^2 +
                val_p[3] * time_relative_s^3 + val_p[4] * time_relative_s^4
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
    str_coefficients = "Shared polynomial coefficients\n" * join(
        ["p$idx = $(round(val_p[idx]; sigdigits=7)) MHz/s^$idx" for idx in 1:4], "\n")
    Label(gl_text[1, 1], str_coefficients;
        color=:black, fontsize=22, halign=:left, valign=:top)
    Label(gl_text[2, 1], "Zero-frequency time offsets";
        color=:black, fontsize=22, halign=:left, valign=:top)
    for (idx_file, id_file) in enumerate(ids_file)
        Label(gl_text[idx_file + 2, 1], "t₀[$id_file] = $(round(time_zero_by_file[id_file]; digits=6)) s";
            color=color_for_file(id_file), fontsize=22, halign=:left, valign=:top)
    end
    Legend(gl_text[length(ids_file) + 3, 1], ax;
        framevisible=false, labelsize=20, patchsize=(36, 24))
    resize_to_layout!(fig)
    mkpath(dirname(abspath(path_png)))
    save(path_png, fig; px_per_unit=2)
    println("Shared polynomial coefficients p1–p4: $(collect(val_p))")
    println("Fitted t0 by index: $time_zero_by_file")
    println("Joint-fit RSS: $(LsqFit.rss(fit_calibration)); converged=$(fit_calibration.converged); parameters=$(length(val_parameters))")
    println("Saved resized calibration PNG to $path_png")
    return (; fig, val_p, points_reference, ids_file, time_zero_by_file)
end

function plot_calibrated_spectrum_results(
    path_data::AbstractString,
    path_png::AbstractString,
    val_p::AbstractVector{<:Real},
    u_calibration_range_s::Tuple{<:Real,<:Real},
)
    config = YAML.load_file(joinpath(path_data, "config.yaml"))
    entries_result = get(config, "result", nothing)
    isnothing(entries_result) && throw(ArgumentError("No 'result' section in config.yaml"))

    freq_min_mhz, freq_max_mhz = -750.0, 750.0
    time_uncertainty_s = 0.001
    axis_width_px = 1440 # 54 / 25.4 * 96
    axis_height_px = 960 # 34 / 25.4 * 96
    fig = Figure(fontsize=14)
    ax = Axis(fig[1, 1];
        width=axis_width_px,
        height=axis_height_px,
        xlabel="Frequency (MHz)",
        ylabel="Spectrum signal (V)",
        xlabelsize=10.667,
        ylabelsize=10.667,
        xticklabelsize=9.333,
        yticklabelsize=9.333,
        limits=((freq_min_mhz, freq_max_mhz), nothing),
    )

    for entry in entries_result
        id_file = Int(entry["index"])
        name_file = "SDS$(lpad(string(id_file), 5, '0'))"
        tag = String(entry["tag"])
        path_jld2 = joinpath(path_data, "$name_file.jld2")
        isfile(path_jld2) || throw(ArgumentError("JLD2 cache does not exist: $path_jld2"))
        payload = JLD2.load(path_jld2)
        ch1_v = payload["ch1_v"]
        ch2_v = payload["ch2_v"]
        time_s = calc_time_axis(payload["metadata_lines"], length(ch2_v))
        time_alignment_zero_s = find_ch1_alignment_zero_crossing(ch1_v, time_s)
        time_aligned_s = time_s .- time_alignment_zero_s
        points_result = NamedTuple[]
        for peak_config in entry["peaks"]
            freq_mhz = Float64(peak_config["freq"])
            time_config_s = Float64(peak_config["time"])
            idx_peak = find_local_peak_near_time(ch2_v, time_aligned_s, time_config_s, time_uncertainty_s)
            push!(points_result, (; freq_mhz, time_config_s, idx_peak,
                time_peak_s=Float64(time_aligned_s[idx_peak]),
                time_peak_raw_s=Float64(time_s[idx_peak])))
        end
        points_zero = filter(point -> iszero(point.freq_mhz), points_result)
        if length(points_zero) > 1
            throw(ArgumentError("Expected at most one 0 MHz result reference for $name_file; found $(length(points_zero))"))
        elseif isempty(points_zero)
            time_zero_s = fit_time_zero_from_references(val_p,
                [point.time_peak_s for point in points_result],
                [point.freq_mhz for point in points_result],
                u_calibration_range_s)
            println("$name_file ($tag): fitted t₀=$(round(time_zero_s; digits=6)) s from $(length(points_result)) nonzero-frequency references")
        else
            time_zero_s = only(points_zero).time_peak_s
            println("$name_file ($tag): detected aligned 0 MHz peak at $(round(time_zero_s; digits=6)) s " *
                "(raw $(round(only(points_zero).time_peak_raw_s; digits=6)) s; CH1 zero $(round(time_alignment_zero_s; digits=6)) s)")
        end
        for point in points_result
            freq_calibrated_mhz = calc_frequency_from_time(val_p, point.time_peak_s - time_zero_s)
            println("$name_file ($tag): reference $(point.freq_mhz) MHz at aligned t=" *
                "$(round(point.time_peak_s; digits=6)) s (raw $(round(point.time_peak_raw_s; digits=6)) s; " *
                "CH2=$(round(ch2_v[point.idx_peak]; sigdigits=6)) V) → " *
                "$(round(freq_calibrated_mhz; digits=3)) MHz " *
                "(residual $(round(point.freq_mhz - freq_calibrated_mhz; digits=3)) MHz)")
        end
        time_relative_s = time_aligned_s .- time_zero_s
        frequency_mhz = val_p[1] .* time_relative_s .+ val_p[2] .* time_relative_s .^ 2 .+
            val_p[3] .* time_relative_s .^ 3 .+ val_p[4] .* time_relative_s .^ 4
        mask_frequency = (u_calibration_range_s[1] .<= time_relative_s) .&
            (time_relative_s .<= u_calibration_range_s[2]) .&
            (freq_min_mhz .<= frequency_mhz) .& (frequency_mhz .<= freq_max_mhz)
        any(mask_frequency) || throw(ArgumentError("No samples for $name_file fall in the requested frequency range"))
        frequency_trace = frequency_mhz[mask_frequency]
        signal_trace_v = ch2_v[mask_frequency]
        color_line = parse(Colorant, String(entry["color"]))
        lines!(ax, frequency_trace, signal_trace_v;
            color=color_line, linewidth=1.2, label=tag)
        for point in points_result
            freq_peak_mhz = calc_frequency_from_time(val_p, point.time_peak_s - time_zero_s)
            val_peak_v = ch2_v[point.idx_peak]
            offset_label_x_mhz = if iszero(point.freq_mhz)
                val_peak_v > 0.3 ? 12.0 : -12.0
            elseif id_file == 21
                point.freq_mhz > 0 ? -12.0 : 12.0
            else
                point.freq_mhz > 0 ? 12.0 : -12.0
            end
            offset_label_y_v = val_peak_v > 0.5 ? -0.012 : 0.012
            scatter!(ax, [freq_peak_mhz], [val_peak_v];
                color=color_line, marker=:utriangle, markersize=8,
                strokecolor=:black, strokewidth=0.8)
            text!(ax, freq_peak_mhz + offset_label_x_mhz, val_peak_v + offset_label_y_v;
                text="t=$(round(point.time_peak_s; digits=3)) s\nf=$(round(freq_peak_mhz; digits=1)) MHz\nV=$(round(val_peak_v; digits=3)) V",
                color=:black, fontsize=14,
                align=(offset_label_x_mhz > 0 ? :left : :right,
                    offset_label_y_v > 0 ? :bottom : :top),
                strokecolor=:white, strokewidth=0.4)
        end
        println("$name_file ($tag): using calibration time range $(u_calibration_range_s[1])…$(u_calibration_range_s[2]) s; " *
            "CH1 alignment zero=$(round(time_alignment_zero_s; digits=6)) s; raw CH2 voltage range=$(round(minimum(signal_trace_v); sigdigits=6))…$(round(maximum(signal_trace_v); sigdigits=6)) V")
        GC.gc()
    end

    axislegend(ax; position=:rt, framevisible=false, labelsize=9.333)
    resize_to_layout!(fig)
    mkpath(dirname(abspath(path_png)))
    save(path_png, fig; px_per_unit=4)
    println("Saved result spectrum PNG with $(round(axis_width_px; digits=1))×$(round(axis_height_px; digits=1)) px axis to $path_png")
    return fig
end

fig = function plot_ch1_spectrum(
    path_jld2::AbstractString,
    path_png::AbstractString,
    path_period_csv::AbstractString,
    path_ch2_png::AbstractString,
)
    payload = JLD2.load(path_jld2)
    metadata_lines = payload["metadata_lines"]
    ch1_v = payload["ch1_v"]
    ch2_v = payload["ch2_v"]
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
length(args_data) <= 1 || throw(ArgumentError(
    "Usage: julia --project=. multidualmot/spectral626.jl [data-folder] [--skip-conversion] [--calibration-only]"))
path_data = isempty(args_data) ? path_data_root : abspath(args_data[1])
if "--skip-conversion" in ARGS
    println("Skipping CSV conversion; using existing JLD2 caches in $path_data")
else
    convert_missing_scope_csvs(path_data)
end

# The follow-up plots and period table remain specific to SDS00014.
if !("--calibration-only" in ARGS)
    path_jld2 = joinpath(path_data, "SDS00014.jld2")
    isfile(path_jld2) || throw(ArgumentError("SDS00014 JLD2 cache is missing after conversion: $path_jld2"))
    path_png = joinpath(path_data, "SDS00014_CH1_vs_time.png")
    path_period_csv = joinpath(path_data, "SDS00014_CH1_periods.csv")
    path_ch2_png = joinpath(path_data, "SDS00014_15_16_19_21_CH2_aligned_periods.png")
    plot_ch1_spectrum(path_jld2, path_png, path_period_csv, path_ch2_png)
    fig_spectrum = plot_ch2_periods_spanning_zero(path_data, [14, 15, 16, 19, 21], path_ch2_png)
end

path_calibration_png = joinpath(path_data, "spectrum626_time_frequency_calibration.png")
calibration_result = plot_spectrum_frequency_calibration(path_data, path_calibration_png)
fig_calibration = calibration_result.fig
path_result_spectrum_png = joinpath(path_data, "spectrum626_calibrated_results.png")
u_calibration_range_s = extrema([
    point.time_peak_s - calibration_result.time_zero_by_file[point.id_file]
    for point in calibration_result.points_reference
])
fig_result_spectrum = plot_calibrated_spectrum_results(
    path_data, path_result_spectrum_png, calibration_result.val_p, u_calibration_range_s)
