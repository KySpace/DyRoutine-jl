using Printf
using JLD2
using CairoMakie

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))

const ROOT_ANLZ = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine"
const OUTPUT_FOLDER = "[06.21].98.FullTime.PrflModlEvol.Tailess.sansMask"
const OUTPUT_DIR = joinpath(ROOT_ANLZ, OUTPUT_FOLDER, "comparison 2")

const SOURCES = [
    (id=96, folder="[06.21].96.FullTime.PrflModlEvol.KeptTail.sansMask"),
    (id=95, folder="[06.21].95.FullTime.PrflModlEvol.KeptTail.avecMask"),
    (id=98, folder=OUTPUT_FOLDER),
    (id=97, folder="[06.21].97.FullTime.PrflModlEvol.Tailess.avecMask"),
]

const SOURCE_PAIR_IDS = ((96, 95), (98, 97))
const VAL_ISTP_PLOT = ["162", "164"]
const VAL_T_HOLD = collect(6:2:200)
const SMWH_CORE = (30, 60)
const PX_IN_UM = 6.5 / 22.06
const STEP_MODL = 1 ./ (2 .* SMWH_CORE .* PX_IN_UM)
const Y_MODL = collect((-SMWH_CORE[2]):SMWH_CORE[2]) .* STEP_MODL[2]

const TAG_CONFIGS = Dict(
    "CFNM" => (
        IB=[5.311, 5.313, 5.316, 5.318, 5.322, 5.325, 5.326, 5.328, 5.332, 5.333, 5.336, 5.338],
        runid=[95, 82, 52, 80, 67, 96, 68, 50, 81, 51, 79, 53],
    ),
    "NTRC" => (
        IB=[5.314, 5.316, 5.318, 5.322, 5.326, 5.332, 5.336, 5.340, 5.343],
        runid=[29, 28, 27, 26, 25, 61, 62, 63, 64],
    ),
)

# Per-IB/per-istp colorrange controls.
# Each matrix must be n_IB x 2, with columns matching VAL_ISTP_PLOT = ["162", "164"].
# Cell values may be:
# - nothing: use the current auto pairwise per-IB limit, shared for source pairs (96,95) and (98,97)
# - upper::Real: use (0, upper)
# - (lower, upper): use that exact range
const COLORRANGE_BY_TAG = Dict(
    "CFNM" => Matrix{Any}(fill(nothing, 12, length(VAL_ISTP_PLOT))),
    "NTRC" => Matrix{Any}(fill(nothing, 9, length(VAL_ISTP_PLOT))),
)

struct CorrCache
    source_id::Int
    path_corr::String
    prfl_evol_stacked
end

function load_corr_cache(source, tag::AbstractString)
    path_corr = joinpath(ROOT_ANLZ, source.folder, @sprintf("%s_corr.jld2", tag))
    isfile(path_corr) || throw(ArgumentError("missing correlation cache: $path_corr"))
    cache = JLD2.load(path_corr)
    haskey(cache, "prfl_evol_stacked") || throw(KeyError("prfl_evol_stacked"))
    return CorrCache(source.id, path_corr, cache["prfl_evol_stacked"])
end

function tag_config(tag::AbstractString)
    haskey(TAG_CONFIGS, tag) || throw(ArgumentError("unknown tag $tag; expected one of $(keys(TAG_CONFIGS))"))
    cfg = TAG_CONFIGS[tag]
    length(cfg.IB) == length(cfg.runid) || throw(DimensionMismatch("$tag IB/runid count mismatch"))
    return cfg
end

function tag_IBs_for(tag::AbstractString, cfg)
    return [@sprintf("%s_%.3f_r%02d", tag, IB, runid) for (IB, runid) in zip(cfg.IB, cfg.runid)]
end

function validate_cache_shapes!(caches::AbstractVector{CorrCache}, tag::AbstractString, cfg)
    n_IB = length(cfg.IB)
    n_istp = length(VAL_ISTP_PLOT)
    for cache in caches
        size(cache.prfl_evol_stacked) == (n_IB, n_istp) ||
            throw(DimensionMismatch("$(cache.path_corr) prfl_evol_stacked size $(size(cache.prfl_evol_stacked)); expected $((n_IB, n_istp))"))
        for idx in eachindex(cache.prfl_evol_stacked)
            size(cache.prfl_evol_stacked[idx]) == (length(Y_MODL), length(VAL_T_HOLD)) ||
                throw(DimensionMismatch("$(cache.path_corr) stacked profile size $(size(cache.prfl_evol_stacked[idx])); expected $((length(Y_MODL), length(VAL_T_HOLD)))"))
        end
    end
    return nothing
end

function finite_max(xs)
    vals = Float64[]
    for x in xs
        append!(vals, filter(isfinite, Float64.(vec(x))))
    end
    isempty(vals) && return 1.0
    val_max = maximum(vals)
    return val_max > 0 ? val_max : 1.0
end

function calc_colorrange_by_pair(caches_by_id, idx_ib::Integer, ids_istp::AbstractVector{<:Integer})
    ranges = Dict{Int,Vector{Tuple{Float64,Float64}}}()
    for pair in SOURCE_PAIR_IDS
        upper_by_istp = [
            finite_max(
                caches_by_id[id].prfl_evol_stacked[idx_ib, idx_istp]
                for id in pair
            )
            for idx_istp in ids_istp
        ]
        for id in pair
            ranges[id] = [(0.0, upper) for upper in upper_by_istp]
        end
    end
    return ranges
end

function parse_manual_colorrange(value)
    isnothing(value) && return nothing
    value isa Real && return (0.0, Float64(value))
    value isa Tuple && length(value) == 2 && return (Float64(value[1]), Float64(value[2]))
    throw(ArgumentError("colorrange entries must be nothing, an upper limit, or a (lower, upper) tuple; got $value"))
end

function manual_colorranges_for(tag::AbstractString, idx_ib::Integer, ids_istp::AbstractVector{<:Integer})
    haskey(COLORRANGE_BY_TAG, tag) || return fill(nothing, length(ids_istp))
    ranges = COLORRANGE_BY_TAG[tag]
    size(ranges, 2) == length(VAL_ISTP_PLOT) ||
        throw(DimensionMismatch("$tag manual colorrange columns $(size(ranges, 2)); expected $(length(VAL_ISTP_PLOT))"))
    idx_ib <= size(ranges, 1) ||
        throw(DimensionMismatch("$tag manual colorrange rows $(size(ranges, 1)); expected at least $idx_ib"))
    return [parse_manual_colorrange(ranges[idx_ib, idx_istp]) for idx_istp in ids_istp]
end

function merge_manual_colorranges(auto_ranges, manual_ranges)
    return [
        isnothing(manual_ranges[i]) ? auto_ranges[i] : manual_ranges[i]
        for i in eachindex(auto_ranges)
    ]
end

function set_compact_stacked_axes!(
    fig::Figure,
    row::Integer,
    idx_ib::Integer,
    idx_source::Integer,
    source_id::Integer;
    width::Real=720,
    height::Real=100,
)
    row_label = 5 * (idx_ib - 1) + idx_source + 1
    Label(fig[row_label, 0], string(source_id); tellwidth=true, tellheight=false, fontsize=10)
    ax_162 = Axis(fig[row_label, 1]; width, height, yticklabelspace=26.0)
    ax_164 = Axis(fig[row_label, 3]; width, height, yticklabelspace=26.0)
    return Dict(
        "repeats" => Array{Axis}(undef, 0, length(VAL_ISTP_PLOT)),
        "stacked" => [ax_162, ax_164],
        "colorbars" => [fig[row_label, 2], fig[row_label, 4]],
    )
end

function tune_compact_axes!(axs::AbstractVector{Axis}, idx_source::Integer)
    for ax in axs
        ax.xlabel = ""
        ax.ylabel = ""
        ax.xticklabelsize = 8
        ax.yticklabelsize = 8
        ax.xticksize = 3
        ax.yticksize = 3
        ax.xgridvisible = false
        ax.ygridvisible = false
        idx_source < length(SOURCES) && hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false)
    end
    return nothing
end

function plot_source_row!(
    fig::Figure,
    cache::CorrCache,
    idx_ib::Integer,
    idx_source::Integer,
    ids_istp::AbstractVector{<:Integer},
    colorranges::AbstractVector{<:Tuple{<:Real,<:Real}},
)
    axs = set_compact_stacked_axes!(fig, 0, idx_ib, idx_source, cache.source_id)
    for (col, idx_istp) in enumerate(ids_istp)
        axs_istp = Dict(
            "repeats" => Array{Axis}(undef, 0, 1),
            "stacked" => [axs["stacked"][col]],
            "colorbars" => [axs["colorbars"][col]],
        )
        plot_prfl_modl_evol!(
            axs_istp,
            Array{Any}(undef, 0, 1),
            [cache.prfl_evol_stacked[idx_ib, idx_istp]],
            VAL_T_HOLD,
            Y_MODL,
            [VAL_ISTP_PLOT[idx_istp]];
            colorrange=colorranges[col],
        )
    end
    tune_compact_axes!(axs["stacked"], idx_source)
    return nothing
end

function make_tag_comparison(tag::AbstractString)
    cfg = tag_config(tag)
    caches = [load_corr_cache(source, tag) for source in SOURCES]
    validate_cache_shapes!(caches, tag, cfg)
    caches_by_id = Dict(cache.source_id => cache for cache in caches)
    ids_istp = collect(eachindex(VAL_ISTP_PLOT))
    tag_IBs = tag_IBs_for(tag, cfg)

    mkpath(OUTPUT_DIR)
    fig = Figure()
    Label(fig[0, 1:4], "$tag stacked modulation sidepeak profile"; tellwidth=false, tellheight=true, halign=:left)
    Label(fig[1, 1], "162"; tellwidth=false, tellheight=true, halign=:center)
    Label(fig[1, 3], "164"; tellwidth=false, tellheight=true, halign=:center)

    for idx_ib in eachindex(tag_IBs)
        row_label = 5 * (idx_ib - 1) + 2
        Label(
            fig[row_label, 1:4],
            string(tag_IBs[idx_ib]);
            tellwidth=false,
            tellheight=true,
            halign=:center,
            fontsize=12,
        )
        auto_ranges = calc_colorrange_by_pair(caches_by_id, idx_ib, ids_istp)
        manual_ranges = manual_colorranges_for(tag, idx_ib, ids_istp)
        for (idx_source, source) in enumerate(SOURCES)
            ranges_source = merge_manual_colorranges(auto_ranges[source.id], manual_ranges)
            plot_source_row!(fig, caches_by_id[source.id], idx_ib, idx_source, ids_istp, ranges_source)
        end
    end

    colgap!(fig.layout, 4)
    rowgap!(fig.layout, 1)
    colsize!(fig.layout, 0, Fixed(24))
    colsize!(fig.layout, 2, Fixed(18))
    colsize!(fig.layout, 4, Fixed(18))
    resize_to_layout!(fig)

    path_pdf = joinpath(OUTPUT_DIR, @sprintf("comparison.[%s].pdf", tag))
    save(path_pdf, fig; backend=CairoMakie)
    println("Wrote $path_pdf")
    return path_pdf
end

function main()
    CairoMakie.activate!()
    for tag in ("CFNM", "NTRC")
        make_tag_comparison(tag)
    end
end

main()
