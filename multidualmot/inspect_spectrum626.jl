using GLMakie
using JLD2

GLMakie.activate!()

const path_spectrum626_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\spectrum 626"
const n_horizontal_divisions = 14
const n_envelope_bins = 3000

function resolve_spectrum_jld2(input::AbstractString)::String
    if isfile(input)
        path_input = abspath(input)
        if endswith(lowercase(path_input), ".csv")
            path_cache = splitext(path_input)[1] * ".jld2"
            isfile(path_cache) || throw(ArgumentError("No JLD2 cache found beside CSV: $path_cache"))
            return path_cache
        end
        endswith(lowercase(path_input), ".jld2") ||
            throw(ArgumentError("Expected an SDS spectrum name, CSV, or JLD2 path; got $input"))
        return path_input
    end

    name = uppercase(splitext(basename(input))[1])
    occursin(r"^SDS\d{5}$", name) ||
        throw(ArgumentError("Expected an SDS name such as SDS00014 or a path; got $input"))
    path_jld2 = joinpath(path_spectrum626_root, "$name.jld2")
    isfile(path_jld2) || throw(ArgumentError("JLD2 cache does not exist: $path_jld2"))
    return path_jld2
end

function parse_metadata_value(metadata_lines::AbstractVector{<:AbstractString}, key::AbstractString)::String
    idx = findfirst(line -> startswith(line, "$key,"), metadata_lines)
    isnothing(idx) && throw(ArgumentError("Spectrum metadata is missing '$key'"))
    fields = split(metadata_lines[idx], ',')
    length(fields) >= 2 || throw(ArgumentError("Could not parse metadata line: $(metadata_lines[idx])"))
    return strip(fields[2])
end

function spectrum_time_axis(metadata_lines::AbstractVector{<:AbstractString}, n_samples::Int)::AbstractRange{Float64}
    units = parse_metadata_value(metadata_lines, "Horizontal Units")
    unit_scale = units == "us" ? 1e-6 : units == "ms" ? 1e-3 : units == "s" ? 1.0 :
        throw(ArgumentError("Unsupported Horizontal Units '$units'; expected us, ms, or s"))
    scale_per_div = parse(Float64, parse_metadata_value(metadata_lines, "Horizontal Scale")) * unit_scale
    half_span_s = (n_horizontal_divisions / 2) * scale_per_div
    return range(-half_span_s, stop=half_span_s, length=n_samples)
end

function minmax_envelope(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, n_bins::Int)
    length(x) == length(y) || throw(DimensionMismatch("Time and signal arrays must have equal lengths"))
    n_bins = min(n_bins, length(x))
    x_plot = Float64[]
    y_plot = Float64[]
    sizehint!(x_plot, 2n_bins)
    sizehint!(y_plot, 2n_bins)

    for idx_bin in 1:n_bins
        idx_first = fld((idx_bin - 1) * length(x), n_bins) + 1
        idx_last = fld(idx_bin * length(x), n_bins)
        y_bin = @view y[idx_first:idx_last]
        y_min, pos_min = findmin(y_bin)
        y_max, pos_max = findmax(y_bin)
        idx_min = idx_first + pos_min - 1
        idx_max = idx_first + pos_max - 1
        if idx_min <= idx_max
            append!(x_plot, (x[idx_min], x[idx_max]))
            append!(y_plot, (y_min, y_max))
        else
            append!(x_plot, (x[idx_max], x[idx_min]))
            append!(y_plot, (y_max, y_min))
        end
    end
    return x_plot, y_plot
end

length(ARGS) <= 1 || throw(ArgumentError("Usage: julia --project=. multidualmot/inspect_spectrum626.jl [SDS00024 | path/to/file.jld2]"))
input_spectrum = isempty(ARGS) ? "SDS00014" : only(ARGS)
path_jld2 = resolve_spectrum_jld2(input_spectrum)
payload = JLD2.load(path_jld2)
metadata_lines = payload["metadata_lines"]
ch1_v = payload["ch1_v"]
ch2_v = payload["ch2_v"]
length(ch1_v) == length(ch2_v) || throw(DimensionMismatch("CH1 and CH2 arrays must have equal lengths"))
time_s = spectrum_time_axis(metadata_lines, length(ch1_v))
x_ch1, y_ch1 = time_s, ch1_v # minmax_envelope(time_s, ch1_v, n_envelope_bins)
x_ch2, y_ch2 = time_s, ch2_v # minmax_envelope(time_s, ch2_v, n_envelope_bins)

fig = Figure(size=(1500, 900), fontsize=22)
ax_ch1 = Axis(fig[1, 1]; title=basename(path_jld2), xlabel="", ylabel="CH1 voltage (V)")
ax_ch2 = Axis(fig[2, 1]; xlabel="Time (s)", ylabel="CH2 voltage (V)")
lines!(ax_ch1, x_ch1, y_ch1; color=:dodgerblue, linewidth=1)
lines!(ax_ch2, x_ch2, y_ch2; color=:darkorange, linewidth=1)
linkxaxes!(ax_ch1, ax_ch2)
display(fig)
println("Showing $path_jld2; close the figure window to exit. No files were saved.")
