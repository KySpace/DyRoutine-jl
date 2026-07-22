using Printf
using JLD2
using CairoMakie

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))

const path_root_anlz = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine"
const folder_output = "[06.21].98.FullTime.PrflModlEvol.Tailess.sansMask"
const path_output = joinpath(path_root_anlz, folder_output, "comparison 2")

const sources = [
    (id=96, folder="[06.21].96.FullTime.PrflModlEvol.KeptTail.sansMask"),
    (id=95, folder="[06.21].95.FullTime.PrflModlEvol.KeptTail.avecMask"),
    (id=98, folder=folder_output),
    (id=97, folder="[06.21].97.FullTime.PrflModlEvol.Tailess.avecMask"),
]
const id_source_monosrc = 97

const ids_source_pair = ((96, 95), (98, 97))
const val_istp_plot = ["162", "164"]
const val_t_hold = collect(6:2:200)
const smwh_core = (30, 60)
const px_in_um = 6.5 / 22.06
const step_modl = 1 ./ (2 .* smwh_core .* px_in_um)
const y_modl = collect((-smwh_core[2]):smwh_core[2]) .* step_modl[2]

const configs_tag = Dict(
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
# Each matrix must be n_IB x 2, with columns matching val_istp_plot = ["162", "164"].
# Cell values may be:
# - nothing: use the current auto pairwise per-IB limit, shared for source pairs (96,95) and (98,97)
# - upper::Real: use (0, upper)
# - (lower, upper): use that exact range
const colorrange_by_tag = Dict(
    "CFNM" => Any[
        (0.0, 0.8) (0.0, 1.0)  # IB 5.311, 162 / 164
        (0.0, 0.8) (0.0, 1.0)  # IB 5.313, 162 / 164
        (0.0, 0.8) (0.0, 1.0)  # IB 5.316, 162 / 164
        (0.0, 0.8) (0.0, 1.0)  # IB 5.318, 162 / 164
        (0.0, 0.8) (0.0, 1.0)  # IB 5.322, 162 / 164
        (0.0, 0.8) (0.0, 1.0)  # IB 5.325, 162 / 164
        (0.0, 0.8) (0.0, 1.0)  # IB 5.326, 162 / 164
        (0.0, 0.8) (0.0, 0.8)  # IB 5.328, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.332, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.333, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.336, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.338, 162 / 164
    ],
    "NTRC" => Any[
        (0.0, 1.0) (0.0, 1.8)  # IB 5.314, 162 / 164
        (0.0, 1.0) (0.0, 1.8)  # IB 5.316, 162 / 164
        (0.0, 1.0) (0.0, 1.8)  # IB 5.318, 162 / 164
        (0.0, 1.0) (0.0, 1.5)  # IB 5.322, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.326, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.332, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.336, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.340, 162 / 164
        (0.0, 0.6) (0.0, 0.8)  # IB 5.343, 162 / 164
    ],
)
# Examples:
# colorrange_by_tag["CFNM"][1, :] .= [(0.0, 0.8), (0.0, 1.2)]  # first CFNM IB, 162/164
# colorrange_by_tag["CFNM"][:, 1] .= [(0.0, 0.9)]              # all CFNM 162 panels
# colorrange_by_tag["NTRC"][3, 2] = 1.4                        # NTRC third IB, 164 => (0, 1.4)
# colorrange_by_tag["NTRC"][5, :] .= [nothing]                 # use auto ranges for one IB

# One multiplier per source row, applied after the IB/istp colorrange is chosen.
const factor_colorrange_source = Dict(
    96 => 2.0,
    95 => 2.0,
    98 => 1.0,
    97 => 1.0,
)

struct CorrCache
    source_id::Int
    path_corr::String
    prfl_evol
    prfl_evol_stacked
end

function load_corr_cache(source, tag::AbstractString)
    path_corr = joinpath(path_root_anlz, source.folder, @sprintf("%s_corr.jld2", tag))
    isfile(path_corr) || throw(ArgumentError("missing correlation cache: $path_corr"))
    cache = JLD2.load(path_corr)
    haskey(cache, "prfl_evol") || throw(KeyError("prfl_evol"))
    haskey(cache, "prfl_evol_stacked") || throw(KeyError("prfl_evol_stacked"))
    return CorrCache(source.id, path_corr, cache["prfl_evol"], cache["prfl_evol_stacked"])
end

function tag_config(tag::AbstractString)
    haskey(configs_tag, tag) || throw(ArgumentError("unknown tag $tag; expected one of $(keys(configs_tag))"))
    cfg = configs_tag[tag]
    length(cfg.IB) == length(cfg.runid) || throw(DimensionMismatch("$tag IB/runid count mismatch"))
    return cfg
end

function tag_IBs_for(tag::AbstractString, cfg)
    return [@sprintf("%s_%.3f_r%02d", tag, IB, runid) for (IB, runid) in zip(cfg.IB, cfg.runid)]
end

function validate_cache_shapes!(caches::AbstractVector{CorrCache}, tag::AbstractString, cfg)
    n_IB = length(cfg.IB)
    n_istp = length(val_istp_plot)
    for cache in caches
        size(cache.prfl_evol_stacked) == (n_IB, n_istp) ||
            throw(DimensionMismatch("$(cache.path_corr) prfl_evol_stacked size $(size(cache.prfl_evol_stacked)); expected $((n_IB, n_istp))"))
        size(cache.prfl_evol, 1) == n_IB ||
            throw(DimensionMismatch("$(cache.path_corr) prfl_evol IB count $(size(cache.prfl_evol, 1)); expected $n_IB"))
        size(cache.prfl_evol, 3) == n_istp ||
            throw(DimensionMismatch("$(cache.path_corr) prfl_evol istp count $(size(cache.prfl_evol, 3)); expected $n_istp"))
        for idx in eachindex(cache.prfl_evol_stacked)
            size(cache.prfl_evol_stacked[idx]) == (length(y_modl), length(val_t_hold)) ||
                throw(DimensionMismatch("$(cache.path_corr) stacked profile size $(size(cache.prfl_evol_stacked[idx])); expected $((length(y_modl), length(val_t_hold)))"))
        end
        for idx in eachindex(cache.prfl_evol)
            size(cache.prfl_evol[idx]) == (length(y_modl), length(val_t_hold)) ||
                throw(DimensionMismatch("$(cache.path_corr) repeat profile size $(size(cache.prfl_evol[idx])); expected $((length(y_modl), length(val_t_hold)))"))
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
    for pair in ids_source_pair
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
    haskey(colorrange_by_tag, tag) || return fill(nothing, length(ids_istp))
    ranges = colorrange_by_tag[tag]
    size(ranges, 2) == length(val_istp_plot) ||
        throw(DimensionMismatch("$tag manual colorrange columns $(size(ranges, 2)); expected $(length(val_istp_plot))"))
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

function scale_colorranges(colorranges, factor::Real)
    factor > 0 || throw(ArgumentError("colorrange factor must be positive, got $factor"))
    return [(lower * factor, upper * factor) for (lower, upper) in colorranges]
end

function set_axis_prfl_compare!(ax::Axis, idx_row::Integer, n_row::Integer)
    ax.xlabel = ""
    ax.ylabel = ""
    ax.xticklabelsize = 7
    ax.yticklabelsize = 7
    ax.xticksize = 3
    ax.yticksize = 3
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xticks = 0:50:200
    ax.yticks = 0:0.2:0.6
    ylims!(ax, (0, 0.6))
    idx_row < n_row && hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false)
    return nothing
end

function draw_prfl_inner_stack!(
    gl::GridLayout,
    caption::AbstractString,
    labels::AbstractVector,
    prfls::AbstractVector,
    istp::AbstractString,
    colorranges::AbstractVector{<:Tuple{<:Real,<:Real}};
    width::Real=2.5 * (maximum(val_t_hold) - minimum(val_t_hold)),
    height::Real=100,
)
    length(labels) == length(prfls) ||
        throw(DimensionMismatch("label count $(length(labels)) does not match profile count $(length(prfls))"))
    length(colorranges) == length(prfls) ||
        throw(DimensionMismatch("colorrange count $(length(colorranges)) does not match profile count $(length(prfls))"))
    Label(gl[1, 1:2], caption; tellwidth=false, tellheight=true, halign=:center, fontsize=10)
    clrmap = gen_clrmap_solo(hue_theme_istp[istp])
    n_row = length(prfls)
    for idx_row in eachindex(prfls)
        row = idx_row + 1
        Label(gl[row, 0], string(labels[idx_row]); tellwidth=true, tellheight=false, fontsize=8)
        ax = Axis(gl[row, 1]; width, height, yticklabelspace=20.0)
        hm = heatmap!(
            ax,
            val_t_hold,
            y_modl,
            prfls[idx_row]';
            colorrange=colorranges[idx_row],
            colormap=clrmap,
        )
        Colorbar(gl[row, 2], hm; width=8, ticklabelsize=7)
        set_axis_prfl_compare!(ax, idx_row, n_row)
    end
    rowgap!(gl, 1)
    colgap!(gl, 2)
    colsize!(gl, 0, Fixed(22))
    colsize!(gl, 2, Fixed(16))
    return nothing
end

function set_outer_prfl_layout!(
    fig::Figure,
    tag::AbstractString,
    tag_IBs::AbstractVector,
    title::AbstractString,
)
    Label(fig[0, 1:length(val_istp_plot)], title; tellwidth=false, tellheight=true, halign=:left)
    for (idx_istp, istp) in enumerate(val_istp_plot)
        Label(fig[1, idx_istp], istp; tellwidth=false, tellheight=true, halign=:center)
    end
    grids = Array{GridLayout}(undef, length(tag_IBs), length(val_istp_plot))
    for idx_ib in eachindex(tag_IBs), idx_istp in eachindex(val_istp_plot)
        idx_istp == 1 && Label(fig[idx_ib+1, 0], string(tag_IBs[idx_ib]); tellwidth=true, tellheight=false, fontsize=10)
        gl = GridLayout()
        fig[idx_ib+1, idx_istp] = gl
        grids[idx_ib, idx_istp] = gl
    end
    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 10)
    return grids
end

function load_tag_context(tag::AbstractString)
    cfg = tag_config(tag)
    caches = [load_corr_cache(source, tag) for source in sources]
    validate_cache_shapes!(caches, tag, cfg)
    caches_by_id = Dict(cache.source_id => cache for cache in caches)
    haskey(caches_by_id, id_source_monosrc) ||
        throw(ArgumentError("id_source_monosrc=$id_source_monosrc is not listed in sources"))
    tag_IBs = tag_IBs_for(tag, cfg)
    return (; cfg, caches_by_id, tag_IBs)
end

function calc_scaled_colorranges(tag::AbstractString, caches_by_id, idx_ib::Integer)
    ids_istp = collect(eachindex(val_istp_plot))
    auto_ranges = calc_colorrange_by_pair(caches_by_id, idx_ib, ids_istp)
    manual_ranges = manual_colorranges_for(tag, idx_ib, ids_istp)
    return Dict(
        source.id => scale_colorranges(
            merge_manual_colorranges(auto_ranges[source.id], manual_ranges),
            get(factor_colorrange_source, source.id, 1.0),
        )
        for source in sources
    )
end

function save_comparison_figure(fig::Figure, name::AbstractString, tag::AbstractString)
    resize_to_layout!(fig)
    path_pdf = joinpath(path_output, @sprintf("%s.[%s].pdf", name, tag))
    save(path_pdf, fig; backend=CairoMakie)
    println("Wrote $path_pdf")
    return path_pdf
end

function make_tag_comparison(tag::AbstractString)
    ctx = load_tag_context(tag)
    mkpath(path_output)
    fig = Figure()
    grids = set_outer_prfl_layout!(fig, tag, ctx.tag_IBs, "$tag | stacked modulation sidepeak profile")
    for idx_ib in eachindex(ctx.tag_IBs), idx_istp in eachindex(val_istp_plot)
        ranges_by_source = calc_scaled_colorranges(tag, ctx.caches_by_id, idx_ib)
        prfls = [ctx.caches_by_id[source.id].prfl_evol_stacked[idx_ib, idx_istp] for source in sources]
        labels = [source.id for source in sources]
        draw_prfl_inner_stack!(
            grids[idx_ib, idx_istp],
            string(ctx.tag_IBs[idx_ib]),
            labels,
            prfls,
            val_istp_plot[idx_istp],
            [ranges_by_source[source.id][idx_istp] for source in sources];
        )
    end
    return save_comparison_figure(fig, "comparison.processing", tag)
end

function make_tag_monosrc_stack(tag::AbstractString)
    ctx = load_tag_context(tag)
    mkpath(path_output)
    fig = Figure()
    grids = set_outer_prfl_layout!(fig, tag, ctx.tag_IBs, "$tag | source $id_source_monosrc stacked profile")
    cache = ctx.caches_by_id[id_source_monosrc]
    for idx_ib in eachindex(ctx.tag_IBs), idx_istp in eachindex(val_istp_plot)
        ranges_by_source = calc_scaled_colorranges(tag, ctx.caches_by_id, idx_ib)
        draw_prfl_inner_stack!(
            grids[idx_ib, idx_istp],
            string(ctx.tag_IBs[idx_ib]),
            ["stack"],
            [cache.prfl_evol_stacked[idx_ib, idx_istp]],
            val_istp_plot[idx_istp],
            [ranges_by_source[id_source_monosrc][idx_istp]],
        )
    end
    return save_comparison_figure(fig, "comparison.monosrc_stack", tag)
end

function make_tag_monosrc_reps_stack(tag::AbstractString)
    ctx = load_tag_context(tag)
    mkpath(path_output)
    fig = Figure()
    grids = set_outer_prfl_layout!(fig, tag, ctx.tag_IBs, "$tag | source $id_source_monosrc repeats and stack")
    cache = ctx.caches_by_id[id_source_monosrc]
    for idx_ib in eachindex(ctx.tag_IBs), idx_istp in eachindex(val_istp_plot)
        ranges_by_source = calc_scaled_colorranges(tag, ctx.caches_by_id, idx_ib)
        labels = [@sprintf("rep %d", r) for r in axes(cache.prfl_evol, 2)]
        push!(labels, "stack")
        prfls = [cache.prfl_evol[idx_ib, r, idx_istp] for r in axes(cache.prfl_evol, 2)]
        push!(prfls, cache.prfl_evol_stacked[idx_ib, idx_istp])
        draw_prfl_inner_stack!(
            grids[idx_ib, idx_istp],
            string(ctx.tag_IBs[idx_ib]),
            labels,
            prfls,
            val_istp_plot[idx_istp],
            fill(ranges_by_source[id_source_monosrc][idx_istp], length(prfls)),
        )
    end
    return save_comparison_figure(fig, "comparison.monosrc_reps_stack", tag)
end

function main()
    CairoMakie.activate!()
    for tag in ("CFNM", "NTRC")
        make_tag_comparison(tag)
        make_tag_monosrc_stack(tag)
        make_tag_monosrc_reps_stack(tag)
    end
end

main()
