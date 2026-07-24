using GLMakie: Grid
using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using CairoMakie: Axis, linewidth
using Colors: Oklch
using LaTeXStrings

function gen_clrmap_posneg(hue_pos, hue_neg)
    return [Oklch(1 - abs(t), 0.24 * abs(t), t > 0 ? hue_pos : hue_neg) |> c -> RGBAf(c) for t in range(-1, 1; length=256)]
end

function plot_num_stat_evo!(
    ax::Axis,
    val_t::AbstractVector,
    stat_number::AbstractVector,
    val_istp;
    hue_shift=0.0,
    label=nothing
)
    hue = hue_theme_istp[val_istp] + hue_shift
    val_mean = map(s -> s[1], stat_number)
    val_err = map(s -> s[2], stat_number)
    clr_line = Oklch(0.62, 0.18, hue) |> c -> RGBAf(c, 0.95)
    clr_bar = Oklch(0.62, 0.18, hue) |> c -> RGBAf(c, 0.40)

    errorbars!(ax, val_t, val_mean, val_err; color=clr_bar)
    lines!(ax, val_t, val_mean; color=clr_line, linewidth=2.6, label=val_istp)
    scatter!(ax, val_t, val_mean; color=clr_bar, markersize=9)
    # return path_plot
end

function set_axis_sidepeak_nvlp!(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer}, panel_setter::Function, runinfo)
    length(n_dim_vars) == 3 || throw(ArgumentError("n_dim_vars must be a 3-tuple"))
    fig = Figure()
    axs_repeats = Array{Dict}(undef, n_dim_vars[1])
    fig[0, 1] = Label(fig, text="$(runinfo.date) $(@sprintf("run%02d", runinfo.runid)) IB=$(@sprintf("%.3f", runinfo.IB))A $(runinfo.tag_head)"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    for r in 1:n_dim_vars[1]
        fig[1, r] = Label(fig, text="repeat $r"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
        print("\r\033[2K\rbuilding axes for side peak trend for repeat $r")
        gl = GridLayout()
        fig[2, r] = gl
        axs_repeats[r] = panel_setter(gl, r)
    end
    println("\r\033[2K\raxes built for trends.")
    Box(fig[2, n_dim_vars[1]+1]; color=:black, strokewidth=0)
    colsize!(fig.layout, n_dim_vars[1] + 1, Fixed(2))
    gl = GridLayout()
    fig[1, n_dim_vars[1]+2] = Label(fig, text="Processed after stacked"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+2] = gl
    axs_stacked = panel_setter(gl, 1)
    gl = GridLayout()
    fig[1, n_dim_vars[1]+3] = Label(fig, text="Reps overlayed"; tellwidth=false, tellheight=false, halign=:center, valign=:bottom)
    fig[2, n_dim_vars[1]+3] = gl
    axs_all = panel_setter(gl, 1; extra=true)
    return fig, Dict("repeats" => axs_repeats, "stacked" => axs_stacked, "all" => axs_all)
end

function set_axis_prfl_modl_evol!(
    val_rep::AbstractVector,
    val_istp::AbstractVector,
    title::AbstractString;
    width::Real=240,
    height::Real=400,
)
    isempty(val_rep) && throw(ArgumentError("val_rep must not be empty"))
    isempty(val_istp) && throw(ArgumentError("val_istp must not be empty"))

    fig = Figure()
    n_rep = length(val_rep)
    n_istp = length(val_istp)
    n_col = n_rep + 2
    axs_repeats = Array{Axis}(undef, n_rep, n_istp)
    axs_stacked = Vector{Axis}(undef, n_istp)
    slots_colorbar = Vector{Any}(undef, n_istp)

    Label(fig[0, 1:n_col]; text=title, tellwidth=false, tellheight=true, halign=:left, valign=:bottom)
    rowsize!(fig.layout, 0, 12)
    for (r, rep) in enumerate(val_rep)
        Label(fig[1, r]; text="repeat $rep", tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
    end
    Box(fig[2:n_istp+1, n_rep+1]; color=:black, strokewidth=0)
    colsize!(fig.layout, n_rep + 1, Fixed(2))
    Label(fig[1, n_rep+2]; text="Processed after stacked", tellwidth=false, tellheight=true, halign=:center, valign=:bottom)

    for (i, istp) in enumerate(val_istp)
        row = i + 1
        Label(fig[row, 0]; text="istp=$istp", tellwidth=true, tellheight=false)
        for r in eachindex(val_rep)
            axs_repeats[r, i] = Axis(fig[row, r]; width, height)
        end
        axs_stacked[i] = Axis(fig[row, n_rep+2]; width, height)
        slots_colorbar[i] = fig[row, n_rep+3]
    end
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 8)
    return fig, Dict("repeats" => axs_repeats, "stacked" => axs_stacked, "colorbars" => slots_colorbar)
end

function plot_prfl_modl_evol!(
    axs::Dict,
    prfl_evol::AbstractArray,
    prfl_evol_stacked::AbstractVector,
    val_t::AbstractVector,
    y_modl::AbstractVector,
    val_istp::AbstractVector;
    colorrange=nothing,
    x_ticks=0:10:210,
    y_ticks=0:0.1:0.6,
    y_lims=(0, 0.6),
)
    axs_repeats = axs["repeats"]
    axs_stacked = axs["stacked"]
    slots_colorbar = axs["colorbars"]
    size(prfl_evol, 1) == size(axs_repeats, 1) ||
        throw(DimensionMismatch("repeat profile count $(size(prfl_evol, 1)) does not match repeat axes $(size(axs_repeats, 1))"))
    size(prfl_evol, 2) == length(val_istp) ||
        throw(DimensionMismatch("profile istp count $(size(prfl_evol, 2)) does not match val_istp length $(length(val_istp))"))
    length(prfl_evol_stacked) == length(val_istp) ||
        throw(DimensionMismatch("stacked profile count $(length(prfl_evol_stacked)) does not match val_istp length $(length(val_istp))"))

    for (i, val_istp) in enumerate(val_istp)
        clrmap = gen_clrmap_solo(hue_theme_istp[val_istp])
        hm = nothing
        for r in axes(prfl_evol, 1)
            ax = axs_repeats[r, i]
            hm = isnothing(colorrange) ?
                heatmap!(ax, val_t, y_modl, prfl_evol[r, i]'; colormap=clrmap) :
                heatmap!(ax, val_t, y_modl, prfl_evol[r, i]'; colorrange, colormap=clrmap)
            ylims!(ax, y_lims)
            ax.xticks = x_ticks
            ax.yticks = y_ticks
        end
        ax = axs_stacked[i]
        hm = isnothing(colorrange) ?
            heatmap!(ax, val_t, y_modl, prfl_evol_stacked[i]'; colormap=clrmap) :
            heatmap!(ax, val_t, y_modl, prfl_evol_stacked[i]'; colorrange, colormap=clrmap)
        Colorbar(slots_colorbar[i], hm)
        ylims!(ax, y_lims)
        ax.xticks = x_ticks
        ax.yticks = y_ticks
    end
    return nothing
end

function plot_prfl_core_evol!(
    axs::Dict,
    prfl_evol::AbstractArray,
    prfl_evol_stacked::AbstractVector,
    val_t::AbstractVector,
    pos::AbstractVector,
    val_istp::AbstractVector;
    colorrange=nothing,
    x_ticks=0:10:210,
    pos_ticks=nothing,
    pos_lims=(-14, 14),
    colormap=nothing,
)
    axs_repeats = axs["repeats"]
    axs_stacked = axs["stacked"]
    slots_colorbar = axs["colorbars"]
    size(prfl_evol) == size(axs_repeats) ||
        throw(DimensionMismatch("core profile size $(size(prfl_evol)) does not match repeat axes $(size(axs_repeats))"))
    length(prfl_evol_stacked) == length(val_istp) ||
        throw(DimensionMismatch("stacked core profile count $(length(prfl_evol_stacked)) does not match val_istp length $(length(val_istp))"))
    for i in eachindex(val_istp)
        clrmap = isnothing(colormap) ? gen_clrmap_solo(hue_theme_istp[val_istp[i]]) : colormap
        draw_heatmap(ax, prfl) = isnothing(colorrange) ?
            heatmap!(ax, val_t, pos, prfl'; colormap=clrmap) :
            heatmap!(ax, val_t, pos, prfl'; colorrange, colormap=clrmap)
        hm = nothing
        for r in axes(prfl_evol, 1)
            ax = axs_repeats[r, i]
            hm = draw_heatmap(ax, prfl_evol[r, i])
            ylims!(ax, pos_lims)
            ax.xticks = x_ticks
            isnothing(pos_ticks) || (ax.yticks = pos_ticks)
        end
        ax = axs_stacked[i]
        hm = draw_heatmap(ax, prfl_evol_stacked[i])
        Colorbar(slots_colorbar[i], hm)
        ylims!(ax, pos_lims)
        ax.xticks = x_ticks
        isnothing(pos_ticks) || (ax.yticks = pos_ticks)
    end
    return nothing
end

"""Build the compact stack/repeat comparison used for excitation profiles.

Each row specification must contain `label`, `evol`, and `pos`, where `evol`
is indexed `(row, istp)` and each entry is a profile matrix.
"""
function set_axis_prfl_comparison!(
    rows::AbstractVector,
    val_istp::AbstractVector,
    title::AbstractString;
    width::Real,
    height::Real,
)
    isempty(rows) && throw(ArgumentError("rows must not be empty"))
    isempty(val_istp) && throw(ArgumentError("val_istp must not be empty"))
    n_col = 2 * length(val_istp)
    fig = Figure()
    Label(fig[0, 1:n_col], title; tellwidth=false, tellheight=true, halign=:left)
    axs = Array{Axis}(undef, length(rows), length(val_istp))
    colorbars = Array{Any}(undef, length(rows), length(val_istp))
    for (idx_row, row_spec) in enumerate(rows)
        Label(fig[idx_row, 0], row_spec.label; tellwidth=true, tellheight=false, fontsize=8)
        for idx_istp in eachindex(val_istp)
            row_height = Float64(hasproperty(row_spec, :height) ? row_spec.height : height)
            idx_col = 2 * idx_istp - 1
            axs[idx_row, idx_istp] = Axis(fig[idx_row, idx_col]; width=Float64(width), height=row_height, yticklabelspace=28.0)
            colorbars[idx_row, idx_istp] = fig[idx_row, idx_col + 1]
        end
    end
    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 2)
    return fig, axs, colorbars
end

function plot_prfl_comparison!(
    fig::Figure,
    axs::AbstractMatrix,
    colorbars::AbstractMatrix,
    rows::AbstractVector,
    val_t::AbstractVector,
    val_istp::AbstractVector;
    ylims=nothing,
    colorrange=nothing,
    x_ticks=0:10:210,
    pos_ticks=nothing,
)
    size(axs) == (length(rows), length(val_istp)) ||
        throw(DimensionMismatch("profile axes size $(size(axs)) does not match rows/isotopes $((length(rows), length(val_istp)))"))
    for (idx_row, row_spec) in enumerate(rows), idx_istp in eachindex(val_istp)
        prfl = row_spec.evol[1, idx_istp]
        clrmap = gen_clrmap_solo(hue_theme_istp[val_istp[idx_istp]])
        hm = isnothing(colorrange) ?
            heatmap!(axs[idx_row, idx_istp], val_t, row_spec.pos, prfl'; colormap=clrmap) :
            heatmap!(axs[idx_row, idx_istp], val_t, row_spec.pos, prfl'; colorrange, colormap=clrmap)
        Colorbar(colorbars[idx_row, idx_istp], hm; width=8.0, ticklabelsize=7.0)
        isnothing(ylims) ? ylims!(axs[idx_row, idx_istp], extrema(row_spec.pos)) : ylims!(axs[idx_row, idx_istp], ylims)
        axs[idx_row, idx_istp].xticks = x_ticks
        isnothing(pos_ticks) || (axs[idx_row, idx_istp].yticks = pos_ticks)
    end
    resize_to_layout!(fig)
    return nothing
end

function set_axis_prfl_comparison_table!(
    val_group::AbstractVector,
    val_istp::AbstractVector,
    title::AbstractString,
)
    isempty(val_group) && throw(ArgumentError("val_group must not be empty"))
    isempty(val_istp) && throw(ArgumentError("val_istp must not be empty"))

    fig = Figure()
    Label(fig[0, 1:length(val_istp)], title; tellwidth=false, tellheight=true, halign=:left)
    grids = Array{GridLayout}(undef, length(val_group), length(val_istp))
    for (idx_group, group) in enumerate(val_group), (idx_istp, istp) in enumerate(val_istp)
        idx_istp == 1 && Label(
            fig[idx_group + 1, 0],
            string(group);
            tellwidth=true,
            tellheight=false,
            fontsize=9,
        )
        idx_group == 1 && Label(
            fig[1, idx_istp],
            string(istp);
            tellwidth=false,
            tellheight=true,
            halign=:center,
            fontsize=9,
        )
        gl = GridLayout()
        fig[idx_group + 1, idx_istp] = gl
        grids[idx_group, idx_istp] = gl
    end
    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 10)
    return fig, grids
end

function calc_prfl_colorrange_auto(
    prfl::AbstractMatrix,
    val_t::AbstractVector,
    pos::AbstractVector;
    selector_t_hold::Function=t -> true,
    selector_pos::Function=x -> true,
)
    size(prfl) == (length(pos), length(val_t)) || throw(DimensionMismatch(
        "profile size $(size(prfl)) must match position/time lengths $((length(pos), length(val_t)))"
    ))
    mask_t = Bool[selector_t_hold(t) for t in val_t]
    mask_pos = Bool[selector_pos(x) for x in pos]
    any(mask_t) || throw(ArgumentError("selector_t_hold selected no time points"))
    any(mask_pos) || throw(ArgumentError("selector_pos selected no profile positions"))
    vals = Float64.(prfl[mask_pos, mask_t])
    vals = filter(isfinite, vals)
    isempty(vals) && return (0.0, 1.0)
    upper = maximum(vals)
    return (0.0, upper > 0 ? upper : 1.0)
end

function plot_prfl_comparison_table!(
    fig::Figure,
    grids::AbstractMatrix,
    rows_by_group::AbstractVector,
    val_t::AbstractVector,
    val_istp::AbstractVector;
    width::Real,
    height::Real,
    ylims=nothing,
    colorrange=nothing,
    selector_t_hold::Function=t -> true,
    selector_pos::Function=x -> true,
    x_ticks=0:10:210,
    pos_ticks=nothing,
)
    size(grids) == (length(rows_by_group), length(val_istp)) || throw(DimensionMismatch(
        "profile table grid size $(size(grids)) does not match groups/isotopes $((length(rows_by_group), length(val_istp)))"
    ))
    for (idx_group, rows) in enumerate(rows_by_group), (idx_istp, istp) in enumerate(val_istp)
        isempty(rows) && throw(ArgumentError("profile table group $idx_group has no rows"))
        gl = grids[idx_group, idx_istp]
        Label(gl[1, 1:2], string(rows[1].group_caption); tellwidth=false, tellheight=true, halign=:center, fontsize=8)
        for (idx_row, row_spec) in enumerate(rows)
            row = idx_row + 1
            Label(gl[row, 0], row_spec.label; tellwidth=true, tellheight=false, fontsize=8)
            profile_config = hasproperty(row_spec, :profile_config) ? row_spec.profile_config : nothing
            width_row = isnothing(profile_config) ? width : profile_config.width_to_time * (maximum(val_t) - minimum(val_t))
            row_height = Float64(hasproperty(row_spec, :height) ? row_spec.height : isnothing(profile_config) ? height : profile_config.height)
            ax = Axis(gl[row, 1]; width=Float64(width_row), height=row_height, yticklabelspace=28.0)
            ax.xtickalign = 1
            ax.ytickalign = 1
            clrmap = gen_clrmap_solo(hue_theme_istp[istp]; thres_alpha=0.1, alpha_base=0.1)
            prfl = row_spec.evol[1, idx_istp]
            colorrange_row = isnothing(profile_config) ? colorrange : profile_config.colorrange
            selector_t_hold_row = hasproperty(row_spec, :selector_t_hold) ? row_spec.selector_t_hold : selector_t_hold
            selector_pos_row = hasproperty(row_spec, :selector_pos) ? row_spec.selector_pos : selector_pos
            colorrange_use = isnothing(colorrange_row) ?
                calc_prfl_colorrange_auto(
                    prfl,
                    val_t,
                    row_spec.pos;
                    selector_t_hold=selector_t_hold_row,
                    selector_pos=selector_pos_row,
                ) : colorrange_row
            hm = heatmap!(ax, val_t, row_spec.pos, prfl'; colorrange=colorrange_use, colormap=clrmap, rasterize=true)
            Colorbar(gl[row, 2], hm; width=8.0, ticklabelsize=7.0, tickalign=1)
            ylims_row = isnothing(profile_config) ? ylims : profile_config.ylims
            isnothing(ylims_row) ? ylims!(ax, extrema(row_spec.pos)) : ylims!(ax, ylims_row)
            ax.xticks = x_ticks
            if hasproperty(row_spec, :hide_x_ticklabels) && row_spec.hide_x_ticklabels
                ax.xticklabelsvisible = false
                ax.xlabelvisible = false
            end
            isnothing(pos_ticks) || (ax.yticks = pos_ticks)
        end
        rowgap!(gl, 1)
        colgap!(gl, 2)
        colsize!(gl, 0, Fixed(48))
        colsize!(gl, 2, Fixed(16))
    end
    resize_to_layout!(fig)
    return nothing
end

function set_axes_2axes!(vals::NamedTuple, panel_setter::Function, runinfo)
    fig = Figure()
    n_dim_vars = vals |> vs -> map(length, vs) |> Tuple
    name_vars = vals |> propertynames |> ns -> string.(ns)
    length(n_dim_vars) == 2 || throw(ArgumentError("Select only 2 vars for to make this table."))
    axs = Array{Dict}(undef, n_dim_vars)
    for (r, v_r) in enumerate(vals[1]), (c, v_c) in enumerate(vals[2])
        print("\r\033[2K\rbuilding axes for $r-$c")
        gl_rc = GridLayout()
        fig[r, c] = gl_rc
        Label(gl_rc[0, 1], text="$(name_vars[1])=$(v_r) | $(name_vars[2])=$(v_c)"; tellwidth=false, tellheight=false, halign=:center, valign=:top)
        gl = GridLayout()
        gl_rc[1, 1] = gl
        axs[r, c] = panel_setter(gl; row_cmpl=n_dim_vars[1] - r, col=c)
        gl_rc |> l -> colgap!(l, 0)
        gl_rc |> l -> rowgap!(l, 0)
        gl_rc |> l -> rowsize!(l, 0, 20)
    end
    fig.layout |> l -> colgap!(l, 8)
    fig.layout |> l -> rowgap!(l, 4)
    println("\r\033[2K\raxes built for the trends")
    return fig, axs
end

function set_panel_single_axis(gl::GridLayout; row_cmpl=0, col=1)
    gl |> clean_gridlayout!
    ax = Axis(gl[1, 1]; width=400, height=200)
    hidedecorations!(ax; label=false, ticklabels=false, ticks=true, grid=true, minorticks=true, minorgrid=false)
    if col == 1
        ax.yticklabelsvisible = true
        ax.ylabelvisible = true
    end
    if row_cmpl == 0
        ax.xticklabelsvisible = true
        ax.xlabelvisible = true
    end
    return Dict("ax" => ax)
end

function set_axis_stack_all!(_, panel_setter::Function, runinfo)
    fig = Figure()
    fig[0, 1] = Label(fig, text="$(runinfo.date) $(@sprintf("run%02d", runinfo.runid)) IB=$(@sprintf("%.3f", runinfo.IB))A $(runinfo.tag_head)"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    gl = GridLayout()
    fig[1, 1] = Label(fig, text="Processed after stacked"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
    fig[2, 1] = gl
    axs_stacked = panel_setter(gl; col=1)
    gl = GridLayout()
    fig[1, 2] = Label(fig, text="Reps overlayed"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
    fig[2, 2] = gl
    axs_all = panel_setter(gl; col=2)
    return fig, Dict("stacked" => axs_stacked, "all" => axs_all)
end

function default_trend_property_specs()
    return [
        (
            name="number",
            ylabel="density sum",
            ylim=nothing,
            selection_key="t_vec_sel_number",
            overlay_evol_col=1,
            variants=[(name="dens-sum", evol_spct=("all", "sel"), color=:theme, label="sum", extra=false)],
        ),
        (
            name="weight",
            ylabel="side peak \nweight",
            ylim=(-0.02, 0.17),
            selection_key="t_vec_sel_sp_weight",
            overlay_evol_col=1,
            variants=[
                (name="fit-weight", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-weight", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="height",
            ylabel="side peak \nheight",
            ylim=(-0.1, 1.1),
            selection_key="t_vec_sel_sp_height",
            overlay_evol_col=1,
            variants=[
                (name="fit-height", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-height", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="width",
            ylabel="side peak \nwidth (μm⁻¹)",
            ylim=(0.02, 0.205),
            selection_key="t_vec_sel_sp_width",
            overlay_evol_col=1,
            variants=[
                (name="fit-width", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-width", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="wavenum",
            ylabel="side peak \nwavenum (μm⁻¹)",
            ylim=(0.22, 0.38),
            selection_key="t_vec_sel_sp_wavenum",
            overlay_evol_col=1,
            variants=[
                (name="fit-wavenum", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-wavenum", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="nvlp-size",
            ylabel="envelope size (μm)",
            ylim=(1, 8),
            selection_key="t_vec_sel_nvlp_size",
            overlay_evol_col=2,
            variants=[
                (name="fit-size-x", evol_spct=("all", "sel"), color=:variant_low, label="fit size radial", extra=false),
                (name="fit-size-y", evol_spct=("all", "sel"), color=:variant_high, label="fit size axial", extra=false),
            ],
        ),
        (
            name="nvlp-cent",
            ylabel="envelope cent (μm)",
            ylim=(1, 8),
            selection_key="t_vec_sel_nvlp_cent",
            overlay_evol_col=2,
            variants=[
                (name="fit-cent-x", evol_spct=("all", "sel"), color=:variant_low, label="fit cent radial", extra=false),
                (name="fit-cent-y", evol_spct=("all", "sel"), color=:variant_high, label="fit cent axial", extra=false),
            ],
        ),
    ]
end

function set_panel_trend_properties!(
    gl::GridLayout,
    property_specs::AbstractVector;
    col::Int=1,
    row_cmpl::Int=0,
    extra::Bool=false,
    align_overlay::Bool=true,
    width_evol::Real=400,
    width_spct::Real=400,
    width_freq::Union{Nothing,Real}=nothing,
    height::Real=200,
)
    isnothing(width_freq) || (width_spct = width_freq)
    gl |> clean_gridlayout!
    isempty(property_specs) && throw(ArgumentError("property_specs must not be empty"))
    kwargs_evol = (width=width_evol, height=height, yticklabelspace=40.0)
    kwargs_spct = (width=width_spct, height=height, yticklabelspace=40.0)
    has_extra_variant = any(spec -> any(v.extra for v in spec.variants), property_specs)
    use_extra_col = extra && (align_overlay || has_extra_variant)
    col_spct = use_extra_col ? 3 : 2
    dict_axs = Dict{String,Axis}()
    for (row, spec) in enumerate(property_specs)
        name = spec.name
        col_evol = use_extra_col ? spec.overlay_evol_col : 1
        dict_axs["evol-$name"] = Axis(gl[row, col_evol]; kwargs_evol..., ylabel=spec.ylabel)
        dict_axs["spct-$name"] = Axis(gl[row, col_spct]; kwargs_spct...)
        if extra && any(v.extra for v in spec.variants)
            dict_axs["evol-extra-$name"] = Axis(gl[row, 2]; kwargs_evol..., ylabel=spec.ylabel)
        end
        rowsize!(gl, row, Fixed(height))
    end
    colsize!(gl, 1, Fixed(width_evol))
    if use_extra_col
        colsize!(gl, 2, Fixed(width_evol))
        colsize!(gl, 3, Fixed(width_spct))
    else
        colsize!(gl, 2, Fixed(width_spct))
    end
    rowgap!(gl, 4)
    colgap!(gl, 4)
    for ax in values(dict_axs)
        hideydecorations!(ax; label=true, ticklabels=false, ticks=false, grid=false, minorticks=false, minorgrid=false)
        hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
        if col == 1
            ax.ylabelvisible = true
        end
    end
    last_name = last(property_specs).name
    if row_cmpl == 0
        for key in ("evol-$last_name", "spct-$last_name", "evol-extra-$last_name")
            haskey(dict_axs, key) || continue
            dict_axs[key].xticklabelsvisible = true
            dict_axs[key].xlabelvisible = true
        end
    end
    dict_axs["evol-$last_name"].xlabel = "t hold (ms)"
    dict_axs["spct-$last_name"].xlabel = "freq (Hz)"
    haskey(dict_axs, "evol-extra-$last_name") && (dict_axs["evol-extra-$last_name"].xlabel = "t hold (ms)")
    return dict_axs
end

function set_panel_trend_sidepeak_nvlp!(gl::GridLayout, col::Int; extra=false, kwargs...)
    return set_panel_trend_properties!(gl, default_trend_property_specs(); col, extra, kwargs...)
end

function set_panel_trend_nvlp!(gl::GridLayout; col::Int=1, row_cmpl::Int=0, extra=false)
    spec = filter(s -> s.name in ("nvlp-size", "sizes"), default_trend_property_specs())
    return set_panel_trend_properties!(gl, spec; col, row_cmpl, extra)
end

function plot_shade_range!(axs, sample, clr)
    step = sample |> diff |> s -> filter(!iszero, s) |> minimum
    head = minimum(sample) - step / 2
    tail = maximum(sample) + step / 2
    for ax in axs
        vspan!(ax, head, tail; color=clr)
    end
end

function trend_variant_color(color_id::Symbol, istp)
    hue_theme = hue_theme_istp[istp]
    color_id == :fit && return :springgreen3
    color_id == :moment && return Oklch(0.52, 0.14, hue_theme)
    color_id == :variant_low && return Oklch(0.52, 0.14, hue_theme - 20)
    color_id == :variant_high && return Oklch(0.52, 0.14, hue_theme + 20)
    color_id == :theme && return Oklch(0.52, 0.14, hue_theme)
    throw(ArgumentError("unknown trend color id $color_id"))
end

function split_property_variant(property_variant::AbstractString)
    parts = split(property_variant, "-"; limit=2)
    length(parts) == 2 || return ("", property_variant)
    return Tuple(parts)
end

function trend_variant_evol_spct(variant)
    hasproperty(variant, :evol_spct) && return variant.evol_spct
    hasproperty(variant, :evol_freq) && return variant.evol_freq
    throw(ArgumentError("trend variant $(variant.name) must define evol_spct"))
end

function plot_trend_property_variant!(
    axs::AbstractDict,
    trend::AbstractDict,
    property_variant::AbstractString,
    istp;
    axis_property::Union{Nothing,AbstractString}=nothing,
    evol_spct::Tuple{<:AbstractString,<:AbstractString}=("all", "sel"),
    color=:theme,
    alpha::Real=1.0,
    label::AbstractString=property_variant,
    evol_extra::Bool=false,
)
    _, property = split_property_variant(property_variant)
    property_axis = isnothing(axis_property) ? property : axis_property
    key_ax_evol = evol_extra ? "evol-extra-$property_axis" : "evol-$property_axis"
    key_ax_spct = "spct-$property_axis"
    key_evol = "evol-$(evol_spct[1])-$property_variant"
    key_spct = "spct-$(evol_spct[2])-$property_variant"
    for key in (key_ax_evol, key_ax_spct, key_evol, key_spct, "t_vec", "freq_query")
        source = startswith(key, "evol-") || startswith(key, "spct-") ? (haskey(axs, key) ? axs : trend) : trend
        haskey(source, key) || throw(KeyError(key))
    end
    clr = color isa Symbol ? trend_variant_color(color, istp) : color
    lines!(axs[key_ax_evol], trend["t_vec"], trend[key_evol]; color=(clr, alpha), label)
    lines!(axs[key_ax_spct], trend["freq_query"], trend[key_spct]; color=(clr, alpha), label)
    return nothing
end

function set_trend_tick_grid!(axs::AbstractDict)
    for (key, ax) in axs
        if startswith(key, "evol-")
            ax.xticks = vcat(10:10:100, 120:20:220)
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(5)
        elseif startswith(key, "spct-")
            ax.xticks = 0:10:200
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(2)
        end
    end
    return nothing
end

function plot_trend_property!(
    axs::AbstractDict,
    trend::AbstractDict,
    spec,
    istp;
    to_clean::Bool=false,
    alpha::Real=1.0,
    to_legend::Bool=false,
    to_overlay::Bool=false,
)
    key_evol = "evol-$(spec.name)"
    key_spct = "spct-$(spec.name)"
    keys_axes = [key_evol, key_spct]
    haskey(axs, "evol-extra-$(spec.name)") && push!(keys_axes, "evol-extra-$(spec.name)")
    if to_clean
        clear_axes!([axs[key] for key in keys_axes])
        hue_theme = hue_theme_istp[istp]
        clr_shade = RGBAf(Oklch(0.95, 0.1, hue_theme), 0.2)
        plot_shade_range!([axs[key_evol]], trend[spec.selection_key], clr_shade)
        haskey(axs, "evol-extra-$(spec.name)") &&
            plot_shade_range!([axs["evol-extra-$(spec.name)"]], trend[spec.selection_key], clr_shade)
    end
    for variant in spec.variants
        plot_trend_property_variant!(
            axs,
            trend,
            variant.name,
            istp;
            axis_property=spec.name,
            evol_spct=trend_variant_evol_spct(variant),
            color=variant.color,
            alpha,
            label=variant.label,
            evol_extra=to_overlay && variant.extra,
        )
    end
    if !isnothing(spec.ylim)
        ylims!(axs[key_evol], spec.ylim...)
        haskey(axs, "evol-extra-$(spec.name)") && ylims!(axs["evol-extra-$(spec.name)"], spec.ylim...)
    end
    to_legend && axislegend(axs[key_spct]; position=:lt, framevisible=false, labelsize=14)
    return nothing
end

function set_axis_trend_property_IB_istp!(
    val_IB::AbstractVector,
    val_istp::AbstractVector,
    n_rep::Integer,
    spec,
    title::AbstractString;
    groups::Tuple=(:repeats, :stacked, :all),
    align_overlay::Bool=false,
    width_evol::Real=400,
    width_spct::Real=400,
    width_freq::Union{Nothing,Real}=nothing,
    height::Real=200,
)
    isnothing(width_freq) || (width_spct = width_freq)
    n_rep > 0 || throw(ArgumentError("n_rep must be positive, got $n_rep"))
    isempty(val_istp) && throw(ArgumentError("val_istp must not be empty"))
    isempty(groups) && throw(ArgumentError("groups must not be empty"))
    valid_groups = (:repeats, :stacked, :all)
    all(group -> group in valid_groups, groups) ||
        throw(ArgumentError("groups must contain only $valid_groups, got $groups"))
    length(unique(groups)) == length(groups) ||
        throw(ArgumentError("groups must not contain duplicates, got $groups"))

    fig = Figure()
    fig[0, 1:length(val_istp)] = Label(fig, title; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    axs_IB_istp = Array{Dict{String,Any}}(undef, length(val_IB), length(val_istp))

    for (idx_istp, istp) in enumerate(val_istp)
        gl_istp = GridLayout()
        fig[1, idx_istp] = gl_istp

        column_plan = NamedTuple[]
        has_repeat_separator = :repeats in groups && any(group -> group in (:stacked, :all), groups)
        for (idx_group, group) in enumerate(groups)
            if group == :repeats && has_repeat_separator && idx_group > 1
                push!(column_plan, (; group=:separator, rep=0, label=""))
            end
            if group == :repeats
                append!(column_plan, [(; group, rep=r, label="repeat $r") for r in 1:n_rep])
            elseif group == :stacked
                push!(column_plan, (; group, rep=0, label="Processed after stacked"))
            else
                push!(column_plan, (; group, rep=0, label="Reps overlayed"))
            end
            if group == :repeats && has_repeat_separator && idx_group < length(groups)
                push!(column_plan, (; group=:separator, rep=0, label=""))
            end
        end
        n_col_block = length(column_plan)
        Label(
            gl_istp[0, 1:n_col_block],
            "istp=$istp";
            tellwidth=false,
            tellheight=true,
            halign=:center,
            valign=:bottom,
        )
        for (idx_col, column) in enumerate(column_plan)
            if column.group == :separator
                Box(
                    gl_istp[2:length(val_IB)+1, idx_col];
                    color=:black,
                    strokewidth=0,
                )
                colsize!(gl_istp, idx_col, Fixed(2))
            else
                Label(
                    gl_istp[1, idx_col],
                    column.label;
                    tellwidth=false,
                    tellheight=true,
                    halign=:center,
                    valign=:bottom,
                )
            end
        end

        for (row_IB, IB) in enumerate(val_IB)
            row = row_IB + 1
            Label(
                gl_istp[row, 0],
                "IB=$(@sprintf("%.3f", IB)) A";
                tellwidth=true,
                tellheight=false,
            )
            axs_groups = Dict{String,Any}()
            axs_repeats = :repeats in groups ? Vector{Dict}(undef, n_rep) : nothing
            for (idx_col, column) in enumerate(column_plan)
                column.group == :separator && continue
                gl_panel = GridLayout()
                gl_istp[row, idx_col] = gl_panel
                if column.group == :repeats
                    axs_repeats[column.rep] = set_panel_trend_properties!(
                        gl_panel,
                        [spec];
                        col=column.rep,
                        row_cmpl=length(val_IB) - row_IB,
                        width_evol,
                        width_spct,
                        height,
                    )
                else
                    axs_groups[string(column.group)] = set_panel_trend_properties!(
                        gl_panel,
                        [spec];
                        col=1,
                        row_cmpl=length(val_IB) - row_IB,
                        extra=column.group == :all,
                        align_overlay,
                        width_evol,
                        width_spct,
                        height,
                    )
                end
            end
            isnothing(axs_repeats) || (axs_groups["repeats"] = axs_repeats)
            axs_IB_istp[row_IB, idx_istp] = axs_groups
        end
        rowgap!(gl_istp, 4)
        colgap!(gl_istp, 8)
    end
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 32)
    return fig, axs_IB_istp
end

function set_axis_trend_variant_IB_istp!(
    val_IB::AbstractVector,
    val_istp::AbstractVector,
    spec,
    variant,
    title::AbstractString;
    width_evol::Real=400,
    width_spct::Real=400,
    width_freq::Union{Nothing,Real}=nothing,
    height::Real=120,
)
    isnothing(width_freq) || (width_spct = width_freq)
    isempty(val_IB) && throw(ArgumentError("val_IB must not be empty"))
    isempty(val_istp) && throw(ArgumentError("val_istp must not be empty"))
    fig = Figure()
    fig[0, 1:(2*length(val_istp))] = Label(
        fig,
        title;
        tellwidth=false,
        tellheight=true,
        halign=:left,
        valign=:top,
    )
    axs_IB_istp = Array{Dict{String,Axis}}(undef, length(val_IB), length(val_istp))
    for (idx_istp, istp) in enumerate(val_istp)
        col_evol = 2 * idx_istp - 1
        col_spct = 2 * idx_istp
        Label(fig[1, col_evol], "istp=$istp evol"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
        Label(fig[1, col_spct], "istp=$istp spct"; tellwidth=false, tellheight=true, halign=:center, valign=:bottom)
        for (idx_IB, IB) in enumerate(val_IB)
            row = idx_IB + 1
            idx_istp == 1 && Label(
                fig[row, 0],
                "IB=$(@sprintf("%.3f", IB)) A";
                tellwidth=true,
                tellheight=false,
            )
            ax_evol = Axis(
                fig[row, col_evol];
                width=width_evol,
                height,
                ylabel=idx_istp == 1 ? spec.ylabel : "",
                yticklabelspace=40.0,
            )
            ax_spct = Axis(
                fig[row, col_spct];
                width=width_spct,
                height,
                yticklabelspace=40.0,
            )
            if idx_IB < length(val_IB)
                hidexdecorations!(ax_evol; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
                hidexdecorations!(ax_spct; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
            end
            idx_istp > 1 && hideydecorations!(ax_evol; label=true, ticklabels=false, ticks=false, grid=false, minorticks=false, minorgrid=false)
            hideydecorations!(ax_spct; label=true, ticklabels=false, ticks=false, grid=false, minorticks=false, minorgrid=false)
            ax_evol.xlabel = "t hold (ms)"
            ax_spct.xlabel = "freq (Hz)"
            axs_IB_istp[idx_IB, idx_istp] = Dict("evol-$(spec.name)" => ax_evol, "spct-$(spec.name)" => ax_spct)
            rowsize!(fig.layout, row, Fixed(height))
        end
        colsize!(fig.layout, col_evol, Fixed(width_evol))
        colsize!(fig.layout, col_spct, Fixed(width_spct))
    end
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 8)
    return fig, axs_IB_istp
end

function trend_variant_fidl_key(property_name::AbstractString, property_variant::AbstractString)
    startswith(property_name, "nvlp-") && return nothing
    startswith(property_variant, "fit-") && return "evol-all-fit-sp-fidl"
    startswith(property_variant, "moment-") && property_name in ("width", "wavenum") && return "evol-all-moment-sp-fidl"
    return nothing
end

function plot_trend_variant_overlay!(
    axs::AbstractDict,
    trend_reps::AbstractVector,
    spec,
    variant,
    istp;
    fit_evol=nothing,
    to_legend::Bool=false,
)
    key_ax_evol = "evol-$(spec.name)"
    key_ax_spct = "spct-$(spec.name)"
    evol_kind, spct_kind = trend_variant_evol_spct(variant)
    key_evol = "evol-$evol_kind-$(variant.name)"
    key_spct = "spct-$spct_kind-$(variant.name)"
    fidl_key = trend_variant_fidl_key(spec.name, variant.name)
    for key in (key_ax_evol, key_ax_spct)
        haskey(axs, key) || throw(KeyError(key))
    end
    clear_axes!([axs[key_ax_evol], axs[key_ax_spct]])
    if !isempty(trend_reps)
        plot_shade_range!([axs[key_ax_evol]], first(trend_reps)[spec.selection_key], RGBAf(Oklch(0.95, 0.1, hue_theme_istp[istp]), 0.2))
    end
    for r in eachindex(trend_reps)
        trend = trend_reps[r]
        for key in (key_evol, key_spct, "t_vec", "freq_query")
            haskey(trend, key) || throw(KeyError(key))
        end
        clr = Oklch(0.86, 0.053, mod(r / 6 - 0.1, 1) * 360)
        if isnothing(fidl_key)
            scatter!(axs[key_ax_evol], trend["t_vec"], trend[key_evol]; color=RGBAf(clr, 1.0), label="rep $r", markersize=5)
        else
            haskey(trend, fidl_key) || throw(KeyError(fidl_key))
            clr_points = [RGBAf(clr, 0.2 + 0.8 * clamp(fidl, 0, 1)) for fidl in trend[fidl_key]]
            scatter!(axs[key_ax_evol], trend["t_vec"], trend[key_evol]; color=clr_points, label="rep $r", markersize=5)
        end
        lines!(axs[key_ax_spct], trend["freq_query"], trend[key_spct]; color=RGBAf(clr, 1.0), label="rep $r")
    end
    if !isnothing(fit_evol)
        !isnothing(fit_evol) && fit_evol.model == :oscillation_decay ||
            throw(ArgumentError("unsupported evol fit model $(fit_evol.model)"))
        t_fit = range(fit_evol.t_fit[1], fit_evol.t_fit[2]; length=256)
        evol_fit = fit_evol_oscillation_decay_model(t_fit, fit_evol.params)
        lines!(axs[key_ax_evol], t_fit, evol_fit; color=(:black, 1), linewidth=1)
        vlines!(axs[key_ax_spct], fit_evol.ν; color=(:black, 1 - fit_evol.rel_residue))
        text!(
            axs[key_ax_evol],
            0.98,
            0.96;
            text=@sprintf("ν %.1f Hz\nλ %.1f ms\nrel_rss %.2f", fit_evol.ν, fit_evol.λ, fit_evol.rel_residue),
            space=:relative,
            align=(:right, :top),
            color=:black,
            fontsize=11,
        )
    end
    if !isnothing(spec.ylim)
        ylims!(axs[key_ax_evol], spec.ylim...)
    end
    set_trend_tick_grid!(axs)
    to_legend && axislegend(axs[key_ax_spct]; position=:lt, framevisible=false, labelsize=14)
    return nothing
end

function set_axis_spectrum_property_IB_istp!(
    val_IB::AbstractVector,
    val_istp::AbstractVector,
    spec,
    title::AbstractString;
    groups::Tuple=(:stacked, :all),
    width::Real=360,
    height::Real=180,
)
    isempty(val_IB) && throw(ArgumentError("val_IB must not be empty"))
    isempty(val_istp) && throw(ArgumentError("val_istp must not be empty"))
    isempty(spec.variants) && throw(ArgumentError("spec $(spec.name) must contain at least one variant"))
    valid_groups = (:stacked, :all)
    all(group -> group in valid_groups, groups) ||
        throw(ArgumentError("spectrum groups must contain only $valid_groups, got $groups"))
    length(unique(groups)) == length(groups) ||
        throw(ArgumentError("spectrum groups must not contain duplicates, got $groups"))

    fig = Figure()
    fig[0, 1:length(val_istp)] = Label(fig, title; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    axs_istp = Dict{Any,Dict{Symbol,Vector{Axis}}}()
    for (idx_istp, istp) in enumerate(val_istp)
        gl_istp = GridLayout()
        fig[1, idx_istp] = gl_istp
        Label(
            gl_istp[0, 1:length(groups)],
            "istp=$istp";
            tellwidth=false,
            tellheight=true,
            halign=:center,
            valign=:bottom,
        )
        for (idx_group, group) in enumerate(groups)
            label = group == :stacked ? "Processed after stacked" : "Mean of rep spectra"
            Label(
                gl_istp[1, idx_group],
                label;
                tellwidth=false,
                tellheight=true,
                halign=:center,
                valign=:bottom,
            )
        end

        axs_group = Dict{Symbol,Vector{Axis}}()
        for (idx_group, group) in enumerate(groups)
            axs_variant = Axis[]
            for (idx_variant, variant) in enumerate(spec.variants)
                ax = Axis(
                    gl_istp[idx_variant+1, idx_group];
                    width,
                    height,
                    xlabel="IB (A)",
                    ylabel=idx_group == 1 ? variant.label : "",
                )
                if idx_variant < length(spec.variants)
                    hidexdecorations!(ax; label=true, ticklabels=true, ticks=false, grid=false, minorticks=false, minorgrid=false)
                end
                if idx_group > 1
                    hideydecorations!(ax; label=true, ticklabels=false, ticks=false, grid=false, minorticks=false, minorgrid=false)
                end
                push!(axs_variant, ax)
                rowsize!(gl_istp, idx_variant + 1, Fixed(height))
            end
            axs_group[group] = axs_variant
            colsize!(gl_istp, idx_group, Fixed(width))
        end
        rowgap!(gl_istp, 4)
        colgap!(gl_istp, 8)
        axs_istp[istp] = axs_group
    end
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 32)
    return fig, axs_istp
end

function calc_freq_IB_matrix(
    trends::AbstractMatrix{<:AbstractDict},
    variant;
    spct_kind::AbstractString="sel",
    freq_kind::Union{Nothing,AbstractString}=nothing,
)
    isnothing(freq_kind) || (spct_kind = freq_kind)
    key_spct = "spct-$spct_kind-$(variant.name)"
    return [
        begin
            haskey(trends[c, i], key_spct) || throw(KeyError(key_spct))
            trends[c, i][key_spct]
        end
        for c in axes(trends, 1), i in axes(trends, 2)
    ] |> stack
end

function plot_spectrum_property_IB_istp!(
    axs_istp::AbstractDict,
    trends_by_group::NamedTuple,
    val_IB::AbstractVector,
    val_istp::AbstractVector,
    spec;
    spct_kind::AbstractString="sel",
    freq_kind::Union{Nothing,AbstractString}=nothing,
    colorrange=(0.3, 1.0),
)
    isnothing(freq_kind) || (spct_kind = freq_kind)
    for (idx_istp, istp) in enumerate(val_istp)
        hue_theme = hue_theme_istp[istp]
        clrmap = gen_clrmap_solo(hue_theme)
        for (group, axs_variant) in axs_istp[istp]
            trends = trends_by_group[group]
            size(trends, 1) == length(val_IB) ||
                throw(DimensionMismatch("group $group has $(size(trends, 1)) IB rows, expected $(length(val_IB))"))
            size(trends, 2) == length(val_istp) ||
                throw(DimensionMismatch("group $group has $(size(trends, 2)) istp columns, expected $(length(val_istp))"))
            for (idx_variant, variant) in enumerate(spec.variants)
                ax = axs_variant[idx_variant]
                mat_spct_IB_istp = calc_freq_IB_matrix(trends, variant; spct_kind)
                freq_query = trends[1, idx_istp]["freq_query"]
                heatmap!(
                    ax,
                    val_IB,
                    freq_query,
                    mat_spct_IB_istp[:, :, idx_istp]';
                    colormap=clrmap,
                    colorrange,
                )
            end
        end
    end
    return nothing
end

function plot_trends_sidepeak!(axs::Dict, trend::Dict, istp; to_clean=false, alpha=1.0, is_stacked=false, to_legend=false, to_overlay=false)
    hue_theme = hue_theme_istp[istp]
    clr_mmt = Oklch(0.52, 0.14, hue_theme)
    clr_fit = (:springgreen3, 1.0)
    clr_theme = Oklch(0.52, 0.14, hue_theme)
    clr_shade_selected = RGBAf(Oklch(0.95, 0.1, hue_theme), 0.2)
    axs_sidepeaks_evol = axs |> a -> matching_axes(a, r"(evol(-extra)?)-(weight|width|height|wavenum)")
    axs_sidepeaks_spct = axs |> a -> matching_axes(a, r"spct-(weight|width|height|wavenum)")
    if to_clean
        axs |> a -> matching_axes(a, r"(freq|(evol(-extra)?))-(number|dens-sum|weight|width|height|wavenum)") |> clear_axes!
        plot_shade_range!(axs_sidepeaks_evol, trend["t_vec_sel_sp"], clr_shade_selected)
    end
    key_number_evol = haskey(axs, "evol-number") ? "evol-number" : "evol-dens-sum"
    lines!(axs[key_number_evol], trend["t_vec"], trend["evol-all-dens-sum"]; color=(clr_theme, alpha), label="sum")
    lines!(axs["evol-weight"], trend["t_vec"], trend["evol-all-fit-weight"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-height"], trend["t_vec"], trend["evol-all-fit-height"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-width"], trend["t_vec"], trend["evol-all-fit-width"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["evol-wavenum"], trend["t_vec"], trend["evol-all-fit-wavenum"]; color=(clr_fit, alpha), label="fit")
    if to_overlay # overlay plots are plotted separately for fit and moment
        lines!(axs["evol-extra-weight"], trend["t_vec"], trend["evol-all-moment-weight"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-extra-height"], trend["t_vec"], trend["evol-all-moment-height"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-extra-width"], trend["t_vec"], trend["evol-all-moment-width"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-extra-wavenum"], trend["t_vec"], trend["evol-all-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    else
        lines!(axs["evol-weight"], trend["t_vec"], trend["evol-all-moment-weight"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-height"], trend["t_vec"], trend["evol-all-moment-height"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-width"], trend["t_vec"], trend["evol-all-moment-width"]; color=(clr_mmt, alpha), label="moment")
        lines!(axs["evol-wavenum"], trend["t_vec"], trend["evol-all-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    end
    lines!(axs["spct-weight"], trend["freq_query"], trend["spct-sel-fit-weight"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["spct-height"], trend["freq_query"], trend["spct-sel-fit-height"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["spct-width"], trend["freq_query"], trend["spct-sel-fit-width"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["spct-wavenum"], trend["freq_query"], trend["spct-sel-fit-wavenum"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["spct-weight"], trend["freq_query"], trend["spct-sel-moment-weight"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["spct-height"], trend["freq_query"], trend["spct-sel-moment-height"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["spct-width"], trend["freq_query"], trend["spct-sel-moment-width"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["spct-wavenum"], trend["freq_query"], trend["spct-sel-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    ylims!(axs["evol-weight"], -0.02, 0.52)
    ylims!(axs["evol-height"], -0.1, 3.1)
    ylims!(axs["evol-width"], 0.02, 0.205)
    ylims!(axs["evol-wavenum"], 0.22, 0.38)
    if to_legend
        axislegend(axs[key_number_evol]; position=:rt, framevisible=false, labelsize=14)
        axislegend(axs["spct-weight"]; position=:lt, framevisible=false, labelsize=14)
    end
    if to_overlay
        ylims!(axs["evol-extra-weight"], -0.02, 0.52)
        ylims!(axs["evol-extra-height"], -0.1, 3.1)
        ylims!(axs["evol-extra-width"], 0.02, 0.205)
        ylims!(axs["evol-extra-wavenum"], 0.22, 0.38)
    end
end

function plot_trends_nvlp!(axs::Dict, trend::Dict, istp; to_clean=false, alpha=1.0, is_stacked=false, to_legend=false, to_overlay=false)
    hue_theme = hue_theme_istp[istp]
    clr_theme = Oklch(0.52, 0.14, hue_theme)
    clr_shade_selected = RGBAf(Oklch(0.95, 0.1, hue_theme), 0.2)
    clr_theme1 = Oklch(0.52, 0.14, hue_theme - 20)
    clr_theme2 = Oklch(0.52, 0.14, hue_theme + 20)
    key_nvlp = haskey(axs, "evol-nvlp-size") ? "nvlp-size" : "sizes"
    if to_clean
        axs |> a -> matching_axes(a, r"(spct|evol)-(nvlp-size|sizes)") |> clear_axes!
        plot_shade_range!([axs["evol-$key_nvlp"]], trend["t_vec_sel_nvlp"], clr_shade_selected)
    end
    lines!(axs["evol-$key_nvlp"], trend["t_vec"], trend["evol-all-fit-size-x"]; color=(clr_theme1, alpha), label="fit")
    lines!(axs["evol-$key_nvlp"], trend["t_vec"], trend["evol-all-fit-size-y"]; color=(clr_theme2, alpha), label="fit")
    lines!(axs["spct-$key_nvlp"], trend["freq_query"], trend["spct-sel-fit-size-x"]; color=(clr_theme1, alpha), label="fit size x")
    lines!(axs["spct-$key_nvlp"], trend["freq_query"], trend["spct-sel-fit-size-y"]; color=(clr_theme2, alpha), label="fit size y")
    ylims!(axs["evol-$key_nvlp"], 1, 11)
    if to_legend
        axislegend(axs["spct-$key_nvlp"]; position=:lt, framevisible=false, labelsize=14)
    end
end

function plot_trend_all!(
    axs_trend::Dict,
    trend_reps::AbstractVector,
    trend_stacked_over_rep::Dict,
    istp;
    property_specs::AbstractVector=default_trend_property_specs(),
)
    if haskey(axs_trend, "repeats")
        length(axs_trend["repeats"]) == length(trend_reps) ||
            throw(DimensionMismatch("repeat axes count $(length(axs_trend["repeats"])) does not match trend count $(length(trend_reps))"))
        for r in eachindex(trend_reps)
            axs = axs_trend["repeats"][r]
            set_trend_tick_grid!(axs)
            for spec in property_specs
                plot_trend_property!(
                    axs,
                    trend_reps[r],
                    spec,
                    istp;
                    to_clean=true,
                    alpha=1.0,
                    to_legend=true,
                    to_overlay=false,
                )
            end
        end
    end

    if haskey(axs_trend, "all")
        axs = axs_trend["all"]
        set_trend_tick_grid!(axs)
        for r in eachindex(trend_reps), spec in property_specs
            plot_trend_property!(
                axs,
                trend_reps[r],
                spec,
                istp;
                to_clean=r == firstindex(trend_reps),
                alpha=0.5,
                to_legend=r == firstindex(trend_reps),
                to_overlay=true,
            )
        end
    end

    if haskey(axs_trend, "stacked")
        axs = axs_trend["stacked"]
        for spec in property_specs
            plot_trend_property!(axs, trend_stacked_over_rep, spec, istp; to_clean=true, alpha=1.0, to_legend=true)
        end
        set_trend_tick_grid!(axs)
    end
    return nothing
end

function plot_trend_nvlp!(axs_trend::Dict, trend_reps::AbstractVector, trend_stacked_over_rep::Dict, istp)
    function set_tick_grid!(axs)
        axs_sidepeaks_evol = axs |> a -> matching_axes(a, r"(evol(-extra)?)-(weight|width|height|wavenum|number|dens-sum|nvlp-size|nvlp-cent|sizes)")
        axs_sidepeaks_spct = axs |> a -> matching_axes(a, r"spct-(weight|width|height|wavenum|number|dens-sum|nvlp-size|nvlp-cent|sizes)")
        for ax in axs_sidepeaks_evol
            ax.xticks = vcat(10:10:100, 120:20:220)
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(5)
        end
        for ax in axs_sidepeaks_spct
            ax.xticks = 0:10:200
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(2)
        end
    end
    axs = axs_trend["all"]
    for r = axes(trend_reps, 1)
        trend = trend_reps[r]
        # plot on both the individual reps and all reps combined
        alpha = 0.5
        # shade the selected time points on all reps combined only once
        to_clean = r == 1
        to_legend = r == 1
        axs |> set_tick_grid!
        plot_trends_nvlp!(axs, trend, istp; to_clean, alpha, to_legend)
    end
    # plot on the stacked axes
    axs = axs_trend["stacked"]
    plot_trends_nvlp!(axs, trend_stacked_over_rep, istp; to_clean=true, alpha=1.0, to_legend=true)
    axs |> set_tick_grid!
end
