using CSV
using GLMakie
using Colors
using DataFrames
using JLD2
using Peaks

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
                strokecolor=:white, strokewidth=2)
        end
        println("$name_file: selected period $(round(time_s[idx_start]; digits=6)) to $(round(time_s[idx_end]; digits=6)) s; rising zero-crossing $(round(time_zero; digits=6)) s; CH2 offset $val_offset_v V; trimmed CH2 baseline $(round(val_ch2_baseline_v; digits=4)) V; annotated $(length(idx_peaks_ch2)) peaks (>$val_min_peak_prominence_v V above baseline, $(val_peak_dead_time_s) s dead-time)")
        GC.gc()
    end

    axislegend(ax; position=:rt, framevisible=false, labelsize=24)
    save(path_png, fig; px_per_unit=1)
    println("Saved $(px_width)×$(px_height) multi-file CH1/CH2 PNG to $path_png")
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

length(ARGS) <= 1 || throw(ArgumentError("Usage: julia --project=. multidualmot/spectral626.jl [data-folder]"))
path_data = isempty(ARGS) ? path_data_root : abspath(ARGS[1])
convert_missing_scope_csvs(path_data)

# The follow-up plots and period table remain specific to SDS00014.
path_jld2 = joinpath(path_data, "SDS00014.jld2")
isfile(path_jld2) || throw(ArgumentError("SDS00014 JLD2 cache is missing after conversion: $path_jld2"))
path_png = joinpath(path_data, "SDS00014_CH1_vs_time.png")
path_period_csv = joinpath(path_data, "SDS00014_CH1_periods.csv")
path_ch2_png = joinpath(path_data, "SDS00014_15_16_19_21_CH2_aligned_periods.png")
plot_ch1_spectrum(path_jld2, path_png, path_period_csv, path_ch2_png)
fig_spectrum = plot_ch2_periods_spanning_zero(path_data, [14, 15, 16, 19, 21], path_ch2_png)
