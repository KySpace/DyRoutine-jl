using CSV
using CairoMakie
using DataFrames
using JLD2

const path_data_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\spectrum 626"
const path_csv_default = joinpath(path_data_root, "SDS00014.csv")
const path_jld2_default = joinpath(path_data_root, "SDS00014.jld2")
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

function calc_minmax_envelope(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, n_bins::Int)::Tuple{Vector{Float64},Vector{Float64}}
    length(x) == length(y) || throw(DimensionMismatch("x and y must have equal lengths; got $(length(x)) and $(length(y))"))
    isempty(x) && throw(ArgumentError("Cannot plot empty data"))
    n_bins > 0 || throw(ArgumentError("n_bins must be positive; got $n_bins"))

    n_bins = min(n_bins, length(x))
    x_envelope = Float64[]
    y_envelope = Float64[]
    sizehint!(x_envelope, 2 * n_bins)
    sizehint!(y_envelope, 2 * n_bins)

    for idx_bin in 1:n_bins
        idx_first = fld((idx_bin - 1) * length(x), n_bins) + 1
        idx_last = fld(idx_bin * length(x), n_bins)
        y_bin = @view y[idx_first:idx_last]
        y_min, pos_min = findmin(y_bin)
        y_max, pos_max = findmax(y_bin)
        idx_min = idx_first + pos_min - 1
        idx_max = idx_first + pos_max - 1

        if idx_min <= idx_max
            append!(x_envelope, (Float64(x[idx_min]), Float64(x[idx_max])))
            append!(y_envelope, (Float64(y_min), Float64(y_max)))
        else
            append!(x_envelope, (Float64(x[idx_max]), Float64(x[idx_min])))
            append!(y_envelope, (Float64(y_max), Float64(y_min)))
        end
    end

    return x_envelope, y_envelope
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

function plot_ch1_spectrum(path_jld2::AbstractString, path_png::AbstractString, path_period_csv::AbstractString)::Nothing
    payload = JLD2.load(path_jld2)
    metadata_lines = payload["metadata_lines"]
    ch1_v = payload["ch1_v"]
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
    x_plot, y_plot = calc_minmax_envelope(time_s, ch1_v, 2 * px_width)

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
    lines!(ax, x_plot, y_plot; color=:dodgerblue, linewidth=1.5)
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

length(ARGS) <= 4 || throw(ArgumentError("Usage: julia --project=. multidualmot/spectral626.jl [input.csv] [output.jld2] [output.png] [periods.csv]"))
path_csv = isempty(ARGS) ? path_csv_default : abspath(ARGS[1])
path_jld2 = length(ARGS) < 2 ? splitext(path_csv)[1] * ".jld2" : abspath(ARGS[2])
path_png = length(ARGS) < 3 ? splitext(path_jld2)[1] * "_CH1_vs_time.png" : abspath(ARGS[3])
path_period_csv = length(ARGS) < 4 ? splitext(path_jld2)[1] * "_CH1_periods.csv" : abspath(ARGS[4])

if isfile(path_jld2)
    println("Using existing JLD2 file: $path_jld2")
else
    convert_scope_csv_to_jld2(path_csv, path_jld2)
end
plot_ch1_spectrum(path_jld2, path_png, path_period_csv)
