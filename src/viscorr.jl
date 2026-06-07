using GLMakie: Grid
using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using CairoMakie: Axis
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
            name="dens-sum",
            ylabel="density sum",
            ylim=nothing,
            selection_key="t_vec_sel_sp",
            overlay_evol_col=1,
            variants=[(name="dens-sum", evol_freq=("all", "sel"), color=:theme, label="sum", extra=false)],
        ),
        (
            name="weight",
            ylabel="side peak \nweight",
            ylim=(-0.02, 0.22),
            selection_key="t_vec_sel_sp",
            overlay_evol_col=1,
            variants=[
                (name="fit-weight", evol_freq=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-weight", evol_freq=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="height",
            ylabel="side peak \nheight",
            ylim=(-0.1, 1.1),
            selection_key="t_vec_sel_sp",
            overlay_evol_col=1,
            variants=[
                (name="fit-height", evol_freq=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-height", evol_freq=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="width",
            ylabel="side peak \nwidth (μm⁻¹)",
            ylim=(0.02, 0.205),
            selection_key="t_vec_sel_sp",
            overlay_evol_col=1,
            variants=[
                (name="fit-width", evol_freq=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-width", evol_freq=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="wavenum",
            ylabel="side peak \nwavenum (μm⁻¹)",
            ylim=(0.22, 0.38),
            selection_key="t_vec_sel_sp",
            overlay_evol_col=1,
            variants=[
                (name="fit-wavenum", evol_freq=("all", "sel"), color=:fit, label="fit", extra=false),
                (name="moment-wavenum", evol_freq=("all", "sel"), color=:moment, label="moment", extra=true),
            ],
        ),
        (
            name="sizes",
            ylabel="envelope size (μm)",
            ylim=(1, 11),
            selection_key="t_vec_sel_nvlp",
            overlay_evol_col=2,
            variants=[
                (name="fit-size-x", evol_freq=("all", "sel"), color=:variant_low, label="fit size x", extra=false),
                (name="fit-size-y", evol_freq=("all", "sel"), color=:variant_high, label="fit size y", extra=false),
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
    width_freq::Real=400,
    height::Real=200,
)
    gl |> clean_gridlayout!
    isempty(property_specs) && throw(ArgumentError("property_specs must not be empty"))
    kwargs_evol = (width=width_evol, height=height, yticklabelspace=40.0)
    kwargs_freq = (width=width_freq, height=height, yticklabelspace=40.0)
    has_extra_variant = any(spec -> any(v.extra for v in spec.variants), property_specs)
    use_extra_col = extra && (align_overlay || has_extra_variant)
    col_freq = use_extra_col ? 3 : 2
    dict_axs = Dict{String,Axis}()
    for (row, spec) in enumerate(property_specs)
        name = spec.name
        col_evol = use_extra_col ? spec.overlay_evol_col : 1
        dict_axs["evol-$name"] = Axis(gl[row, col_evol]; kwargs_evol..., ylabel=spec.ylabel)
        dict_axs["freq-$name"] = Axis(gl[row, col_freq]; kwargs_freq...)
        if extra && any(v.extra for v in spec.variants)
            dict_axs["evol-extra-$name"] = Axis(gl[row, 2]; kwargs_evol..., ylabel=spec.ylabel)
        end
        rowsize!(gl, row, Fixed(height))
    end
    colsize!(gl, 1, Fixed(width_evol))
    if use_extra_col
        colsize!(gl, 2, Fixed(width_evol))
        colsize!(gl, 3, Fixed(width_freq))
    else
        colsize!(gl, 2, Fixed(width_freq))
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
        for key in ("evol-$last_name", "freq-$last_name", "evol-extra-$last_name")
            haskey(dict_axs, key) || continue
            dict_axs[key].xticklabelsvisible = true
            dict_axs[key].xlabelvisible = true
        end
    end
    dict_axs["evol-$last_name"].xlabel = "t hold (ms)"
    dict_axs["freq-$last_name"].xlabel = "freq (Hz)"
    haskey(dict_axs, "evol-extra-$last_name") && (dict_axs["evol-extra-$last_name"].xlabel = "t hold (ms)")
    return dict_axs
end

function set_panel_trend_sidepeak_nvlp!(gl::GridLayout, col::Int; extra=false, kwargs...)
    return set_panel_trend_properties!(gl, default_trend_property_specs(); col, extra, kwargs...)
end

function set_panel_trend_nvlp!(gl::GridLayout; col::Int=1, row_cmpl::Int=0, extra=false)
    spec = filter(s -> s.name == "sizes", default_trend_property_specs())
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

function plot_trend_property_variant!(
    axs::AbstractDict,
    trend::AbstractDict,
    property_variant::AbstractString,
    istp;
    axis_property::Union{Nothing,AbstractString}=nothing,
    evol_freq::Tuple{<:AbstractString,<:AbstractString}=("all", "sel"),
    color=:theme,
    alpha::Real=1.0,
    label::AbstractString=property_variant,
    evol_extra::Bool=false,
)
    _, property = split_property_variant(property_variant)
    property_axis = isnothing(axis_property) ? property : axis_property
    key_ax_evol = evol_extra ? "evol-extra-$property_axis" : "evol-$property_axis"
    key_ax_freq = "freq-$property_axis"
    key_evol = "evol-$(evol_freq[1])-$property_variant"
    key_freq = "freq-$(evol_freq[2])-$property_variant"
    for key in (key_ax_evol, key_ax_freq, key_evol, key_freq, "t_vec", "freq_query")
        source = startswith(key, "evol-") || startswith(key, "freq-") ? (haskey(axs, key) ? axs : trend) : trend
        haskey(source, key) || throw(KeyError(key))
    end
    clr = color isa Symbol ? trend_variant_color(color, istp) : color
    lines!(axs[key_ax_evol], trend["t_vec"], trend[key_evol]; color=(clr, alpha), label)
    lines!(axs[key_ax_freq], trend["freq_query"], trend[key_freq]; color=(clr, alpha), label)
    return nothing
end

function set_trend_tick_grid!(axs::AbstractDict)
    for (key, ax) in axs
        if startswith(key, "evol-")
            ax.xticks = 0:50:200
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(5)
        elseif startswith(key, "freq-")
            ax.xticks = 0:10:100
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
    key_freq = "freq-$(spec.name)"
    keys_axes = [key_evol, key_freq]
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
            evol_freq=variant.evol_freq,
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
    to_legend && axislegend(axs[key_freq]; position=:lt, framevisible=false, labelsize=14)
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
    width_freq::Real=400,
    height::Real=200,
)
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
                        width_freq,
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
                        width_freq,
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

function plot_trends_sidepeak!(axs::Dict, trend::Dict, istp; to_clean=false, alpha=1.0, is_stacked=false, to_legend=false, to_overlay=false)
    hue_theme = hue_theme_istp[istp]
    clr_mmt = Oklch(0.52, 0.14, hue_theme)
    clr_fit = (:springgreen3, 1.0)
    clr_theme = Oklch(0.52, 0.14, hue_theme)
    clr_shade_selected = RGBAf(Oklch(0.95, 0.1, hue_theme), 0.2)
    axs_sidepeaks_evol = axs |> a -> matching_axes(a, r"(evol(-extra)?)-(weight|width|height|wavenum)")
    axs_sidepeaks_freq = axs |> a -> matching_axes(a, r"freq-(weight|width|height|wavenum)")
    if to_clean
        axs |> a -> matching_axes(a, r"(freq|(evol(-extra)?))-(dens-sum|weight|width|height|wavenum)") |> clear_axes!
        plot_shade_range!(axs_sidepeaks_evol, trend["t_vec_sel_sp"], clr_shade_selected)
    end
    lines!(axs["evol-dens-sum"], trend["t_vec"], trend["evol-all-dens-sum"]; color=(clr_theme, alpha), label="sum")
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
    lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-fit-weight"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-fit-height"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-fit-width"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-fit-wavenum"]; color=(clr_fit, alpha), label="fit")
    lines!(axs["freq-weight"], trend["freq_query"], trend["freq-sel-moment-weight"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-height"], trend["freq_query"], trend["freq-sel-moment-height"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-width"], trend["freq_query"], trend["freq-sel-moment-width"]; color=(clr_mmt, alpha), label="moment")
    lines!(axs["freq-wavenum"], trend["freq_query"], trend["freq-sel-moment-wavenum"]; color=(clr_mmt, alpha), label="moment")
    ylims!(axs["evol-weight"], -0.02, 0.52)
    ylims!(axs["evol-height"], -0.1, 3.1)
    ylims!(axs["evol-width"], 0.02, 0.205)
    ylims!(axs["evol-wavenum"], 0.22, 0.38)
    if to_legend
        axislegend(axs["evol-dens-sum"]; position=:rt, framevisible=false, labelsize=14)
        axislegend(axs["freq-weight"]; position=:lt, framevisible=false, labelsize=14)
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
    if to_clean
        axs |> a -> matching_axes(a, r"(freq|evol)-(sizes)") |> clear_axes!
        plot_shade_range!([axs["evol-sizes"]], trend["t_vec_sel_nvlp"], clr_shade_selected)
    end
    lines!(axs["evol-sizes"], trend["t_vec"], trend["evol-all-fit-size-x"]; color=(clr_theme1, alpha), label="fit")
    lines!(axs["evol-sizes"], trend["t_vec"], trend["evol-all-fit-size-y"]; color=(clr_theme2, alpha), label="fit")
    lines!(axs["freq-sizes"], trend["freq_query"], trend["freq-sel-fit-size-x"]; color=(clr_theme1, alpha), label="fit size x")
    lines!(axs["freq-sizes"], trend["freq_query"], trend["freq-sel-fit-size-y"]; color=(clr_theme2, alpha), label="fit size y")
    ylims!(axs["evol-sizes"], 1, 11)
    if to_legend
        axislegend(axs["freq-sizes"]; position=:lt, framevisible=false, labelsize=14)
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
        axs_sidepeaks_evol = axs |> a -> matching_axes(a, r"(evol(-extra)?)-(weight|width|height|wavenum|dens-sum|sizes)")
        axs_sidepeaks_freq = axs |> a -> matching_axes(a, r"freq-(weight|width|height|wavenum|dens-sum|sizes)")
        for ax in axs_sidepeaks_evol
            ax.xticks = 0:50:200
            ax.xminorticksvisible = true
            ax.xminorgridvisible = true
            ax.xminorticks = IntervalsBetween(5)
        end
        for ax in axs_sidepeaks_freq
            ax.xticks = 0:10:100
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
