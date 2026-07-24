using Printf
using JLD2
using CairoMakie

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))

const path_root_simu = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
const folder_simu_selected = "Anlz.17.Simu-03.[2025.07.22].[←16]"
const path_simu_selected = joinpath(path_root_simu, folder_simu_selected)
const path_output = joinpath(path_simu_selected, "Prfl Evol")

const tag_simu = "SIMU-NTRC"
const name_param = "as_s"
const val_as_s = [80, 85, 90, 95, 100, 105]
const val_istp_plot = ["162", "164"]

const height_prfl_default = 100.0

struct SimuCorrCache
    path_corr::String
    meta
    prfl_modl_evol
    prfl_modl_evol_stacked
    prfl_axial_evol
    prfl_axial_evol_stacked
    prfl_radial_evol
    prfl_radial_evol_stacked
end

function load_simu_corr_cache(tag::AbstractString)
    path_corr = joinpath(path_simu_selected, @sprintf("%s_corr.jld2", tag))
    isfile(path_corr) || throw(ArgumentError("missing correlation cache: $path_corr"))
    cache = JLD2.load(path_corr)
    for key in (
        "meta_corr",
        "prfl_evol",
        "prfl_evol_stacked",
        "prfl_axial_evol",
        "prfl_axial_evol_stacked",
        "prfl_radial_evol",
        "prfl_radial_evol_stacked",
    )
        haskey(cache, key) || throw(KeyError(key))
    end
    return SimuCorrCache(
        path_corr,
        cache["meta_corr"],
        cache["prfl_evol"],
        cache["prfl_evol_stacked"],
        cache["prfl_axial_evol"],
        cache["prfl_axial_evol_stacked"],
        cache["prfl_radial_evol"],
        cache["prfl_radial_evol_stacked"],
    )
end

function validate_simu_cache!(cache::SimuCorrCache)
    collect(cache.meta.val_vars.IB) == val_as_s ||
        throw(ArgumentError("cache $(cache.path_corr) has $(cache.meta.val_vars.IB), expected $val_as_s"))
    collect(cache.meta.val_vars.istp) == val_istp_plot ||
        throw(ArgumentError("cache $(cache.path_corr) has istp $(cache.meta.val_vars.istp), expected $val_istp_plot"))
    n_param = length(val_as_s)
    n_rep = length(cache.meta.val_vars.rep)
    n_istp = length(val_istp_plot)
    for (name, evol) in (
        (:modl, cache.prfl_modl_evol),
        (:axial, cache.prfl_axial_evol),
        (:radial, cache.prfl_radial_evol),
    )
        size(evol) == (n_param, n_rep, n_istp) ||
            throw(DimensionMismatch("$name evol size $(size(evol)); expected $((n_param, n_rep, n_istp))"))
    end
    for (name, stacked) in (
        (:modl, cache.prfl_modl_evol_stacked),
        (:axial, cache.prfl_axial_evol_stacked),
        (:radial, cache.prfl_radial_evol_stacked),
    )
        size(stacked) == (n_param, n_istp) ||
            throw(DimensionMismatch("$name stacked size $(size(stacked)); expected $((n_param, n_istp))"))
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

function calc_colorrange_auto(prfls::AbstractVector)
    return (0.0, finite_max(prfls))
end

function calc_height_same_unit(pos::AbstractVector, pos_ref::AbstractVector; height_ref::Real=height_prfl_default)
    span = maximum(pos) - minimum(pos)
    span_ref = maximum(pos_ref) - minimum(pos_ref)
    span_ref > 0 || throw(ArgumentError("reference position span must be positive"))
    return height_ref * span / span_ref
end

function set_axis_prfl_evol!(
    ax::Axis,
    idx_row::Integer,
    n_row::Integer,
    pos::AbstractVector;
    y_lims=extrema(pos),
)
    ax.xlabel = ""
    ax.ylabel = ""
    ax.xticklabelsize = 7
    ax.yticklabelsize = 7
    ax.xticksize = 3
    ax.yticksize = 3
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xticks = 0:50:250
    ax.yticks = LinearTicks(3)
    ylims!(ax, y_lims)
    idx_row < n_row && hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false)
    return nothing
end

function set_outer_layout!(fig::Figure, cache::SimuCorrCache, title::AbstractString)
    Label(fig[0, 1:length(val_istp_plot)], title; tellwidth=false, tellheight=true, halign=:left)
    for (idx_istp, istp) in enumerate(val_istp_plot)
        Label(fig[1, idx_istp], istp; tellwidth=false, tellheight=true, halign=:center)
    end
    grids = Array{GridLayout}(undef, length(val_as_s), length(val_istp_plot))
    for idx_as_s in eachindex(val_as_s), idx_istp in eachindex(val_istp_plot)
        idx_istp == 1 && Label(
            fig[idx_as_s+1, 0],
            @sprintf("%s=%s", name_param, val_as_s[idx_as_s]);
            tellwidth=true,
            tellheight=false,
            fontsize=10,
        )
        gl = GridLayout()
        fig[idx_as_s+1, idx_istp] = gl
        grids[idx_as_s, idx_istp] = gl
    end
    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 10)
    return grids
end

function draw_prfl_rows!(
    gl::GridLayout,
    caption::AbstractString,
    rows::AbstractVector,
    istp::AbstractString,
)
    Label(gl[1, 1:2], caption; tellwidth=false, tellheight=true, halign=:center, fontsize=10)
    clrmap = gen_clrmap_solo(hue_theme_istp[istp])
    n_row = length(rows)
    for (idx_row, row_spec) in enumerate(rows)
        row = idx_row + 1
        Label(gl[row, 0], row_spec.label; tellwidth=true, tellheight=false, fontsize=8)
        ax = Axis(gl[row, 1]; width=row_spec.width, height=row_spec.height, yticklabelspace=28.0)
        hm = heatmap!(
            ax,
            row_spec.val_t,
            row_spec.pos,
            row_spec.prfl';
            colorrange=row_spec.colorrange,
            colormap=clrmap,
        )
        Colorbar(gl[row, 2], hm; width=8, ticklabelsize=7)
        set_axis_prfl_evol!(ax, idx_row, n_row, row_spec.pos; y_lims=row_spec.y_lims)
    end
    rowgap!(gl, 1)
    colgap!(gl, 2)
    colsize!(gl, 0, Fixed(36))
    colsize!(gl, 2, Fixed(16))
    return nothing
end

function make_modl_rows(prfls::AbstractVector, labels::AbstractVector, val_t, pos_modl)
    colorrange = calc_colorrange_auto(prfls)
    width_prfl = 2.5 * (maximum(val_t) - minimum(val_t))
    return [
        (
            label=string(labels[i]),
            prfl=prfls[i],
            val_t,
            pos=pos_modl,
            colorrange,
            width=width_prfl,
            height=height_prfl_default,
            y_lims=(0.0, 0.6),
        )
        for i in eachindex(prfls)
    ]
end

function make_core_rows(prfls_axial::AbstractVector, prfls_radial::AbstractVector, labels::AbstractVector, val_t, pos_axial, pos_radial)
    length(prfls_axial) == length(prfls_radial) ||
        throw(DimensionMismatch("axial count $(length(prfls_axial)) does not match radial count $(length(prfls_radial))"))
    colorrange_axial = calc_colorrange_auto(prfls_axial)
    colorrange_radial = calc_colorrange_auto(prfls_radial)
    width_prfl = 2.5 * (maximum(val_t) - minimum(val_t))
    height_radial = calc_height_same_unit(pos_radial, pos_axial)
    rows = Any[]
    for i in eachindex(prfls_axial)
        push!(
            rows,
            (
                label="$(labels[i]) axial",
                prfl=prfls_axial[i],
                val_t,
                pos=pos_axial,
                colorrange=colorrange_axial,
                width=width_prfl,
                height=height_prfl_default,
                y_lims=extrema(pos_axial),
            ),
        )
        push!(
            rows,
            (
                label="$(labels[i]) radial",
                prfl=prfls_radial[i],
                val_t,
                pos=pos_radial,
                colorrange=colorrange_radial,
                width=width_prfl,
                height=height_radial,
                y_lims=extrema(pos_radial),
            ),
        )
    end
    return rows
end

function save_profile_figure(fig::Figure, name::AbstractString, tag::AbstractString)
    mkpath(path_output)
    resize_to_layout!(fig)
    path_pdf = joinpath(path_output, @sprintf("%s.[%s].pdf", name, tag))
    save(path_pdf, fig; backend=CairoMakie)
    println("Wrote $path_pdf")
    return path_pdf
end

function make_profile_pdf!(cache::SimuCorrCache, name::AbstractString, title::AbstractString, row_builder::Function)
    fig = Figure()
    grids = set_outer_layout!(fig, cache, title)
    for idx_as_s in eachindex(val_as_s), idx_istp in eachindex(val_istp_plot)
        rows = row_builder(idx_as_s, idx_istp)
        draw_prfl_rows!(
            grids[idx_as_s, idx_istp],
            @sprintf("%s=%s", name_param, val_as_s[idx_as_s]),
            rows,
            val_istp_plot[idx_istp],
        )
    end
    return save_profile_figure(fig, name, cache.meta.tag)
end

function make_modl_stack_pdf!(cache::SimuCorrCache)
    val_t = collect(cache.meta.val_vars.t_hold)
    pos_modl = collect(cache.meta.y_modl)
    row_builder = (idx_as_s, idx_istp) -> make_modl_rows(
        [cache.prfl_modl_evol_stacked[idx_as_s, idx_istp]],
        ["stack"],
        val_t,
        pos_modl,
    )
    return make_profile_pdf!(cache, "prfl_modl_stack", "$(cache.meta.tag) | modulation profile stack", row_builder)
end

function make_modl_reps_stack_pdf!(cache::SimuCorrCache)
    val_t = collect(cache.meta.val_vars.t_hold)
    pos_modl = collect(cache.meta.y_modl)
    row_builder = (idx_as_s, idx_istp) -> begin
        labels = [@sprintf("rep %d", r) for r in axes(cache.prfl_modl_evol, 2)]
        push!(labels, "stack")
        prfls = [cache.prfl_modl_evol[idx_as_s, r, idx_istp] for r in axes(cache.prfl_modl_evol, 2)]
        push!(prfls, cache.prfl_modl_evol_stacked[idx_as_s, idx_istp])
        make_modl_rows(prfls, labels, val_t, pos_modl)
    end
    return make_profile_pdf!(cache, "prfl_modl_reps_stack", "$(cache.meta.tag) | modulation profile reps and stack", row_builder)
end

function make_core_stack_pdf!(cache::SimuCorrCache)
    val_t = collect(cache.meta.val_vars.t_hold)
    pos_axial = collect(((-cache.meta.smwh_core[2]):cache.meta.smwh_core[2]) .* cache.meta.px_in_um[2])
    pos_radial = collect(((-cache.meta.smwh_core[1]):cache.meta.smwh_core[1]) .* cache.meta.px_in_um[1])
    row_builder = (idx_as_s, idx_istp) -> make_core_rows(
        [cache.prfl_axial_evol_stacked[idx_as_s, idx_istp]],
        [cache.prfl_radial_evol_stacked[idx_as_s, idx_istp]],
        ["stack"],
        val_t,
        pos_axial,
        pos_radial,
    )
    return make_profile_pdf!(cache, "prfl_core_stack", "$(cache.meta.tag) | axial/radial profile stack", row_builder)
end

function make_core_reps_stack_pdf!(cache::SimuCorrCache)
    val_t = collect(cache.meta.val_vars.t_hold)
    pos_axial = collect(((-cache.meta.smwh_core[2]):cache.meta.smwh_core[2]) .* cache.meta.px_in_um[2])
    pos_radial = collect(((-cache.meta.smwh_core[1]):cache.meta.smwh_core[1]) .* cache.meta.px_in_um[1])
    row_builder = (idx_as_s, idx_istp) -> begin
        labels = [@sprintf("rep %d", r) for r in axes(cache.prfl_axial_evol, 2)]
        push!(labels, "stack")
        prfls_axial = [cache.prfl_axial_evol[idx_as_s, r, idx_istp] for r in axes(cache.prfl_axial_evol, 2)]
        prfls_radial = [cache.prfl_radial_evol[idx_as_s, r, idx_istp] for r in axes(cache.prfl_radial_evol, 2)]
        push!(prfls_axial, cache.prfl_axial_evol_stacked[idx_as_s, idx_istp])
        push!(prfls_radial, cache.prfl_radial_evol_stacked[idx_as_s, idx_istp])
        make_core_rows(prfls_axial, prfls_radial, labels, val_t, pos_axial, pos_radial)
    end
    return make_profile_pdf!(cache, "prfl_core_reps_stack", "$(cache.meta.tag) | axial/radial profile reps and stack", row_builder)
end

function main()
    CairoMakie.activate!()
    cache = load_simu_corr_cache(tag_simu)
    validate_simu_cache!(cache)
    make_modl_stack_pdf!(cache)
    make_modl_reps_stack_pdf!(cache)
    make_core_stack_pdf!(cache)
    make_core_reps_stack_pdf!(cache)
end

main()
