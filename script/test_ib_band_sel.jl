using GLMakie
using CairoMakie
using JLD2
using Printf
using Dates
using Statistics: mean
using Colors: Oklch

GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "fitmodels.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))

# Select the extraction and correlation folders for one tag. SVGs use path_output.
# For simulation caches, path_essn and path_corr may point to the same folder.
# path_essn = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine\[07.24].100.Extr"
# path_corr = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine\[07.24].103.PrflProc.NumberMask.[←100]"
path_essn = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations\Anlz.20.Simu-05.[2026.09.04].Cache"                      
path_corr = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations\Anlz.20.Simu-05.[2026.09.04].Cache"
path_output = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine\[10.13].104.DynRangeSel.[←100,103]"
tag = "SIMU-NTRC"
initial_property = "nvlp-size-axial"
initial_IB = nothing # Use the middle IB when nothing.
ylim_property_overrides = (
    nvlp_size_axial=(3.0, 12.0),  # Example: (4.0, 9.0) for a fixed range in μm.
    nvlp_size_radial=(1.0, 4.0),
)

function load_ib_band_data(path_essn::AbstractString, path_corr::AbstractString, tag::AbstractString)
    path_extr = joinpath(path_essn, "$(tag)_essn_extr.jld2")
    isfile(path_extr) || throw(ArgumentError("extraction cache does not exist: $path_extr"))
    meta_extr, extr_fmt, extr_stacked_over_rep = JLD2.jldopen(path_extr, "r") do file
        file["meta_extr"], file["extr_fmt"], file["extr_stacked_over_rep"]
    end
    val_vars = meta_extr.val_vars
    n_ib, n_rep, n_t, n_istp = size(extr_fmt)
    size(extr_stacked_over_rep) == (n_ib, n_t, n_istp) || throw(DimensionMismatch(
        "stacked extraction size $(size(extr_stacked_over_rep)) must be $((n_ib, n_t, n_istp))"
    ))
    length(val_vars.IB) == n_ib || throw(DimensionMismatch("IB count does not match extraction"))
    length(val_vars.t_hold) == n_t || throw(DimensionMismatch("time count does not match extraction"))
    length(val_vars.istp) == n_istp || throw(DimensionMismatch("isotope count does not match extraction"))
    n_istp == 2 || throw(ArgumentError("this window requires two istp values, found $n_istp"))
    n_t >= 2 || throw(ArgumentError("at least two t_hold points are required"))
    issorted(val_vars.t_hold) || throw(ArgumentError("t_hold must be sorted for the range slider"))

    path_corr_file = joinpath(path_corr, "$(tag)_corr.jld2")
    isfile(path_corr_file) || throw(ArgumentError("correlation cache does not exist: $path_corr_file"))
    meta_corr = JLD2.jldopen(file -> file["meta_corr"], path_corr_file, "r")
    for key in (:IB, :rep, :t_hold, :istp)
        getproperty(meta_corr.val_vars, key) == getproperty(val_vars, key) ||
            throw(DimensionMismatch("$key values differ between $path_extr and $path_corr_file"))
    end
    config_corr = hasproperty(meta_corr, :config_corr) ? meta_corr.config_corr : nothing
    return (; meta_extr, extr_fmt, extr_stacked_over_rep, config_corr, path_corr_file)
end

function profile_kind(spec, variant)
    spec.name == "number" && return :axial
    if spec.name == "nvlp-center"
        return variant.name == "fit-cent-x" ? :radial : :axial
    end
    startswith(spec.name, "nvlp-size-") &&
        return endswith(spec.name, "radial") ? :radial_norm : :axial_norm
    return :modl
end

function load_profile_evolution(path_corr_file::AbstractString, kind::Symbol,
    n_ib::Integer, n_rep::Integer, n_istp::Integer)
    prefix = kind === :modl ? "prfl_evol" : "prfl_$(replace(string(kind), "_norm" => ""))_evol"
    key_rep = kind in (:axial_norm, :radial_norm) ? "$(prefix)_norm" : prefix
    key_stack = "$(key_rep)_stacked"
    reps, stack = JLD2.jldopen(path_corr_file, "r") do file
        haskey(file, key_rep) || throw(ArgumentError("$path_corr_file is missing $key_rep"))
        haskey(file, key_stack) || throw(ArgumentError("$path_corr_file is missing $key_stack"))
        file[key_rep], file[key_stack]
    end
    size(reps) == (n_ib, n_rep, n_istp) || throw(DimensionMismatch(
        "$key_rep size $(size(reps)) must be $((n_ib, n_rep, n_istp))"))
    size(stack) == (n_ib, n_istp) || throw(DimensionMismatch(
        "$key_stack size $(size(stack)) must be $((n_ib, n_istp))"))
    return (; reps, stack)
end

function spectrum_selector_key(spec)
    name = spec.name
    startswith(name, "nvlp-") && return :nvlp
    name == "number" && return :number
    name == "wavenum" && return :sp_wavenum
    name == "width" && return :sp_width
    name == "height" && return :sp_height
    name == "weight" && return :sp_weight
    throw(ArgumentError("no spectrum selector for property $name"))
end

function default_time_range(config_corr, spec, val_t::AbstractVector)
    isnothing(config_corr) && return (first(val_t), last(val_t))
    hasproperty(config_corr, :selector_t_spectrum_val) || return (first(val_t), last(val_t))
    selectors = config_corr.selector_t_spectrum_val
    key = spectrum_selector_key(spec)
    hasproperty(selectors, key) || return (first(val_t), last(val_t))
    selected = getproperty(selectors, key)
    isempty(selected) && return (first(val_t), last(val_t))
    return (
        val_t[argmin(abs.(val_t .- minimum(selected)))],
        val_t[argmin(abs.(val_t .- maximum(selected)))],
    )
end

function property_value(extr::SoloExtract, name::AbstractString)
    name == "dens-sum" && return extr.essentials.sum_dens_full
    if startswith(name, "fit-")
        field = Symbol(replace(name[5:end], "-" => "_"))
        if field in (:size_x, :size_y, :cent_x, :cent_y)
            isnothing(extr.envelope) && throw(ArgumentError("missing envelope fit"))
            envelope = extr.envelope
            pair = startswith(string(field), "size") ? envelope.params_asymm.size : envelope.params_asymm.cent
            return pair[endswith(string(field), "x") ? 1 : 2]
        end
        isnothing(extr.sidepeak) && throw(ArgumentError("missing sidepeak fit"))
        return getproperty(extr.sidepeak.params_tailess, field)
    end
    if startswith(name, "moment-")
        isnothing(extr.sidepeak) && throw(ArgumentError("missing sidepeak moments"))
        return getproperty(extr.sidepeak.moments, Symbol(name[8:end]))
    end
    throw(ArgumentError("unsupported property variant $name"))
end

function property_fidelity(extr::SoloExtract, spec, variant)
    key = spec.name in ("width", "wavenum") ? trend_variant_fidl_key(spec.name, variant.name) : nothing
    isnothing(key) && return nothing
    isnothing(extr.sidepeak) && throw(ArgumentError("missing sidepeak fidelity"))
    sidepeak = extr.sidepeak
    if key == "evol-all-fit-sp-fidl"
        return sidepeak.params_tailess.rel_residue > 0.6 ? 0.0 : 1.0
    end
    return sidepeak.moments.weight < 0.05 ? 0.0 : 1.0
end

function build_property_evolutions(extr_fmt::AbstractArray{<:SoloExtract,4},
    extr_stacked::AbstractArray{<:SoloExtract,3}, spec, variant)
    evol_rep = map(extr -> Float64(property_value(extr, variant.name)), extr_fmt)
    evol_stack = map(extr -> Float64(property_value(extr, variant.name)), extr_stacked)
    fidl_rep = !(spec.name in ("width", "wavenum")) || isnothing(trend_variant_fidl_key(spec.name, variant.name)) ? nothing :
        map(extr -> property_fidelity(extr, spec, variant), extr_fmt)
    fidl_stack = isnothing(fidl_rep) ? nothing :
        map(extr -> property_fidelity(extr, spec, variant), extr_stacked)
    return (; evol_rep, evol_stack, fidl_rep, fidl_stack)
end

function calc_ib_band(data, val_t::AbstractVector, freq_query::AbstractVector,
    time_range::Tuple{<:Real,<:Real}, mode::Symbol; query_weight_kwargs::NamedTuple=NamedTuple())
    lo, hi = extrema(time_range)
    mask_t = (lo .<= val_t) .& (val_t .<= hi)
    any(mask_t) || throw(ArgumentError("selected time range $time_range contains no time points"))
    n_ib, n_rep, _, n_istp = size(data.evol_rep)
    band_weight = Array{Float64}(undef, n_ib, length(freq_query), n_istp)
    band_phase = similar(band_weight)
    for i in 1:n_istp, c in 1:n_ib
        if mode === :stacked
            weight = isnothing(data.fidl_stack) ? nothing : @view(data.fidl_stack[c, :, i])
            coeff = query_spectrum_complex(@view(data.evol_stack[c, :, i]), mask_t, val_t, freq_query;
                weight, query_weight_kwargs...)
            band_weight[c, :, i] .= abs2.(coeff)
            band_phase[c, :, i] .= angle.(coeff)
        elseif mode === :mean_reps
            coeff_reps = [
                query_spectrum_complex(@view(data.evol_rep[c, r, :, i]), mask_t, val_t, freq_query;
                    weight=isnothing(data.fidl_rep) ? nothing : @view(data.fidl_rep[c, r, :, i]),
                    query_weight_kwargs...)
                for r in 1:n_rep
            ]
            band_weight[c, :, i] .= mean([abs2.(coeff) for coeff in coeff_reps])
            # Average phases on the circle, weighted by each rep's normalized spectral power.
            phase_vector = reduce(+, (abs.(coeff) .* coeff for coeff in coeff_reps))
            band_phase[c, :, i] .= angle.(phase_vector)
        else
            throw(ArgumentError("unknown spectrum mode $mode"))
        end
    end
    return (; weight=band_weight, phase=band_phase)
end

function color_ib_band_phase(weight::AbstractMatrix{<:Real}, phase::AbstractMatrix{<:Real})
    size(weight) == size(phase) || throw(DimensionMismatch(
        "weight size $(size(weight)) must match phase size $(size(phase))"
    ))
    return [
        clr_parabola(clamp(weight[idx], 0, 1), mod(rad2deg(phase[idx]), 360), 0.75, 0.12, 0.75)
        for idx in CartesianIndices(weight)
    ] |> colors -> reshape(colors, size(weight))
end

function save_svg_preserve_axis_limits(path::AbstractString, fig::Figure, axs::AbstractVector{<:Axis})
    limits_saved = [ax.limits[] for ax in axs]
    targets_saved = [ax.targetlimits[] for ax in axs]
    try
        for (ax, target) in zip(axs, targets_saved)
            lo, hi = minimum(target), maximum(target)
            ax.limits[] = ((lo[1], hi[1]), (lo[2], hi[2]))
        end
        CairoMakie.save(path, fig; backend=CairoMakie)
    finally
        for (ax, limits) in zip(axs, limits_saved)
            ax.limits[] = limits
        end
        for (ax, target) in zip(axs, targets_saved)
            ax.targetlimits[] = target
        end
    end
end

function run_ib_band_sel(; path_essn::AbstractString=path_essn,
    path_corr::AbstractString=path_corr, path_output::AbstractString=path_output,
    tag::AbstractString=tag, ylim_property_overrides::NamedTuple=ylim_property_overrides)
    cache = load_ib_band_data(path_essn, path_corr, tag)
    meta = cache.meta_extr
    val_ib = Float64.(meta.val_vars.IB)
    val_t = Float64.(meta.val_vars.t_hold)
    val_istp = string.(meta.val_vars.istp)
    freq_query = collect(meta.freq_query)
    specs = [(spec, variant) for spec in meta.trend_property_specs for variant in spec.variants]
    isempty(specs) && throw(ArgumentError("extraction metadata contains no property variants"))
    n_rep = size(cache.extr_fmt, 2)
    n_ib = length(val_ib)
    n_istp = length(val_istp)
    pos_axial = ((-meta.smwh_core[2]):meta.smwh_core[2]) .* meta.px_in_um
    pos_radial = ((-meta.smwh_core[1]):meta.smwh_core[1]) .* meta.px_in_um
    pos_by_kind = Dict(
        :axial => pos_axial, :radial => pos_radial,
        :axial_norm => pos_axial, :radial_norm => pos_radial,
        :modl => meta.y_modl,
    )
    idx_property = findfirst(pair -> pair[1].name == initial_property, specs)
    idx_property = isnothing(idx_property) ? 1 : idx_property
    idx_ib = isnothing(initial_IB) ? cld(length(val_ib), 2) : findfirst(==(initial_IB), val_ib)
    isnothing(idx_ib) && throw(ArgumentError("initial_IB=$initial_IB is absent"))
    idx_rep = n_rep + 1 # Last position is the fitted stack.
    idx_band_mode = 1
    band_modes = (:stacked, :mean_reps)
    ranges_by_property = Dict(
        spec.name => Float64.(default_time_range(cache.config_corr, spec, val_t)) |> Tuple
        for (spec, _) in specs
    )
    evol_cache = Dict{String,Any}()
    evol_axis_cache = Dict{String,Any}()
    profile_cache = Dict{Symbol,Any}()
    function get_profile(idx::Int)
        kind = profile_kind(specs[idx]...)
        data = get!(profile_cache, kind) do
            load_profile_evolution(cache.path_corr_file, kind, n_ib, n_rep, n_istp)
        end
        return kind, data, pos_by_kind[kind]
    end
    function get_data(idx::Int)
        spec, variant = specs[idx]
        get!(evol_cache, "$(spec.name)/$(variant.name)") do
            build_property_evolutions(cache.extr_fmt, cache.extr_stacked_over_rep, spec, variant)
        end
    end
    function get_evolution_axis(idx::Int)
        spec, variant = specs[idx]
        get!(evol_axis_cache, "$(spec.name)/$(variant.name)") do
            data = get_data(idx)
            values = filter(isfinite, vcat(vec(data.evol_rep), vec(data.evol_stack)))
            isempty(values) && throw(ArgumentError("property $(spec.name) has no finite fitted values"))
            lo_data, hi_data = extrema(values)
            pad = max(0.05 * (hi_data - lo_data), 1e-6)
            lo_ref, hi_ref = isnothing(spec.ylim) ? (lo_data, hi_data) : spec.ylim
            key = Symbol(replace(spec.name, "-" => "_"))
            override = hasproperty(ylim_property_overrides, key) ? getproperty(ylim_property_overrides, key) : nothing
            limits = if isnothing(override)
                (min(lo_ref, lo_data - pad), max(hi_ref, hi_data + pad))
            else
                override isa Tuple{<:Real,<:Real} || throw(ArgumentError(
                    "ylim_property_overrides.$key must be nothing or a pair of real limits, got $override"))
                lo, hi = Float64.(override)
                isfinite(lo) && isfinite(hi) && lo < hi || throw(ArgumentError(
                    "ylim_property_overrides.$key must be finite and increasing, got $override"))
                (lo, hi)
            end
            ticks = GLMakie.Makie.get_tickvalues(WilkinsonTicks(5; k_min=4), limits...)
            (; limits, ticks)
        end
    end
    function get_evolution(data, c::Int, r::Int, i::Int)
        r > n_rep ? copy(@view(data.evol_stack[c, :, i])) : copy(@view(data.evol_rep[c, r, :, i]))
    end

    time_range = ranges_by_property[specs[idx_property][1].name]
    obs_property = Observable(idx_property)
    obs_ib = Observable(idx_ib)
    obs_rep = Observable(idx_rep)
    obs_band_mode = Observable(idx_band_mode)
    obs_color_mode = Observable(1)
    obs_show_ib_marker = Observable(true)
    obs_time_range = Observable(Float64.(time_range) |> Tuple)
    data_now = get_data(idx_property)
    obs_evol = [Observable(get_evolution(data_now, idx_ib, idx_rep, i)) for i in 1:2]
    obs_evol_reps = [Observable(get_evolution(data_now, idx_ib, r, i)) for r in 1:n_rep, i in 1:2]
    query_weight_kwargs = isnothing(cache.config_corr) ? NamedTuple() :
        (hasproperty(cache.config_corr, :query_weight_kwargs) ? cache.config_corr.query_weight_kwargs : NamedTuple())
    band_now = calc_ib_band(data_now, val_t, freq_query, obs_time_range[], band_modes[idx_band_mode]; query_weight_kwargs)
    obs_band_weight = [Observable(band_now.weight[:, :, i]) for i in 1:2]
    obs_band_phase_color = [Observable(color_ib_band_phase(band_now.weight[:, :, i], band_now.phase[:, :, i])) for i in 1:2]

    fig = Figure(size=(1080, 1030))
    ctrl = GridLayout(fig[1, 1:2])
    property_controls = GridLayout(ctrl[1, 1:2])
    Button(property_controls[1, 1]; label="←", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_property[] = mod1(obs_property[] - 1, length(specs))
    end
    Label(property_controls[1, 2], lift(i -> begin
        spec, variant = specs[i]
        "property\n$(spec.name) / $(variant.label)  ($i/$(length(specs)))"
    end, obs_property); tellwidth=false, halign=:center, justification=:center)
    Button(property_controls[1, 3]; label="→", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_property[] = mod1(obs_property[] + 1, length(specs))
    end
    for (col, width) in enumerate((36, 350, 36))
        colsize!(property_controls, col, Fixed(width))
    end
    Label(ctrl[2, 1], lift(r -> @sprintf("spectrum span  %.3g–%.3g ms", r...), obs_time_range))
    slider_range = IntervalSlider(ctrl[2, 2]; range=val_t, startvalues=obs_time_range[], tellwidth=false)
    on(slider_range.interval) do interval
        selected = (Float64(interval[1]), Float64(interval[2]))
        ranges_by_property[specs[obs_property[]][1].name] = selected
        obs_time_range[] = selected
    end
    on(obs_property) do i
        range_saved = ranges_by_property[specs[i][1].name]
        set_close_to!(slider_range, range_saved...)
    end

    upper = GridLayout(fig[2, 1:2])
    upper_controls = GridLayout(upper[1, 1])
    Button(upper_controls[1, 1]; label="←", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_ib[] = mod1(obs_ib[] - 1, length(val_ib))
    end
    Label(upper_controls[1, 2], lift(i -> @sprintf("IB\n%.4g A  (%d/%d)", val_ib[i], i, length(val_ib)), obs_ib);
        tellwidth=false, halign=:center, justification=:center)
    Button(upper_controls[1, 3]; label="→", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_ib[] = mod1(obs_ib[] + 1, length(val_ib))
    end
    Button(upper_controls[1, 4]; label="←", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_rep[] = mod1(obs_rep[] - 1, n_rep + 2)
    end
    Label(upper_controls[1, 5], lift(i -> "rep\n" *
        (i == n_rep + 1 ? "stack" : i == n_rep + 2 ? "all" : "$i"), obs_rep);
        tellwidth=false, halign=:center, justification=:center)
    Button(upper_controls[1, 6]; label="→", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_rep[] = mod1(obs_rep[] + 1, n_rep + 2)
    end
    for (col, width) in enumerate((36, 150, 36, 36, 90, 36))
        colsize!(upper_controls, col, Fixed(width))
    end
    obs_profile_ylabel = Observable("")
    Label(upper[2:3, 0], obs_profile_ylabel; rotation=π/2, tellheight=false)
    axs_profile = [Axis(upper[i+1, 1]; ylabel=val_istp[i], xticklabelsvisible=false,
        xticksvisible=false, xlabelvisible=false) for i in 1:2]
    ax_evol = Axis(upper[4, 1]; xlabel="t hold (ms)", ylabel=specs[idx_property][1].ylabel,
        title="per-solo fitted property")
    colsize!(upper, 1, Relative(0.88))
    rowsize!(upper, 2, Fixed(108))
    rowsize!(upper, 3, Fixed(108))
    linkxaxes!(axs_profile..., ax_evol)
    linkyaxes!(axs_profile...)
    profile_plots = Any[nothing, nothing]
    span_lo = lift(r -> r[1], obs_time_range)
    span_hi = lift(r -> r[2], obs_time_range)
    vspan!(ax_evol, span_lo, span_hi; color=(:lightgray, 0.35))
    for i in 1:2
        clr = Oklch(0.52, 0.14, hue_theme_istp[val_istp[i]])
        lines!(ax_evol, val_t, obs_evol[i]; color=clr, linewidth=1.4,
            visible=lift(r -> r != n_rep + 2, obs_rep))
        scatter!(ax_evol, val_t, obs_evol[i]; color=clr, markersize=3,
            visible=lift(r -> r != n_rep + 2, obs_rep))
        for r in 1:n_rep
            lines!(ax_evol, val_t, obs_evol_reps[r, i]; color=(clr, 0.65), linewidth=1.0,
                visible=lift(mode -> mode == n_rep + 2, obs_rep))
        end
    end
    axislegend(ax_evol,
        [LineElement(color=Oklch(0.52, 0.14, hue_theme_istp[val_istp[i]]), linewidth=1.4) for i in 1:2],
        val_istp; position=:rt, framevisible=false)
    set_trend_tick_grid!(Dict("evol-property" => ax_evol))

    lower = GridLayout(fig[3, 1:2])
    lower_controls = GridLayout(lower[2, 1])
    Button(lower_controls[1, 1]; label="←", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_band_mode[] = mod1(obs_band_mode[] - 1, length(band_modes))
    end
    Label(lower_controls[1, 2], lift(i -> "spectrum mode\n" *
        (i == 1 ? "Processed after stacked" : "Mean of rep spectra"), obs_band_mode);
        tellwidth=false, halign=:center, justification=:center)
    Button(lower_controls[1, 3]; label="→", width=36, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_band_mode[] = mod1(obs_band_mode[] + 1, length(band_modes))
    end
    Button(lower_controls[1, 4];
        label=lift(show -> show ? "IB marker: on" : "IB marker: off", obs_show_ib_marker),
        width=125, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_show_ib_marker[] = !obs_show_ib_marker[]
    end
    Button(lower_controls[1, 5];
        label=lift(i -> i == 1 ? "Color: isotope" : "Color: FT phase", obs_color_mode),
        width=140, tellwidth=false).clicks |> clicks -> on(clicks) do _
        obs_color_mode[] = mod1(obs_color_mode[] + 1, 2)
    end
    label_save = Label(lower_controls[1, 7], "○"; width=26, tellwidth=false, halign=:center, color=:gray50)
    button_reset_band = Button(lower_controls[1, 8]; label="Reset axes", width=100, tellwidth=false)
    Button(lower_controls[1, 6]; label="Save SVG").clicks |> clicks -> on(clicks) do _
        spec, variant = specs[obs_property[]]
        name = "$(tag)_ib_band_$(spec.name)_$(variant.name)_$(Dates.format(now(), "yyyymmdd_HHMMSS")).svg"
        path_save = joinpath(path_output, name)
        label_save.text = "…"
        label_save.color = :gray50
        try
            save_svg_preserve_axis_limits(path_save, fig, Axis[axs_profile..., ax_evol, axs_band...])
            label_save.text = "✓"
            label_save.color = :green4
            println("Saved $path_save")
        catch err
            label_save.text = "✕"
            label_save.color = :firebrick3
            println(stderr, "SVG save failed: $(sprint(showerror, err))")
        end
    end
    for (col, width) in enumerate((36, 220, 36, 125, 140, 100, 26, 100))
        colsize!(lower_controls, col, Fixed(width))
    end
    band_layout = GridLayout(lower[1, 1])
    axs_band = [Axis(band_layout[1, i]; xlabel="IB (A)", ylabel=i == 1 ? "freq (Hz)" : "",
        title="istp=$(val_istp[i])") for i in 1:2]
    colorrange = hasproperty(meta, :trend_spectrum_IB_plot_kwargs) ?
        get(meta.trend_spectrum_IB_plot_kwargs, :colorrange, (0.3, 1.0)) : (0.3, 1.0)
    for i in 1:2
        heatmap!(axs_band[i], val_ib, freq_query, obs_band_weight[i];
            colormap=gen_clrmap_solo(hue_theme_istp[val_istp[i]]), colorrange,
            visible=lift(mode -> mode == 1, obs_color_mode))
        heatmap!(axs_band[i], val_ib, freq_query, obs_band_phase_color[i];
            visible=lift(mode -> mode == 2, obs_color_mode))
        vlines!(axs_band[i], lift(idx -> [val_ib[idx]], obs_ib);
            color=:black, linestyle=:dash, linewidth=1.5, visible=obs_show_ib_marker)
    end
    x_pad_left = length(val_ib) > 1 ? (val_ib[2] - val_ib[1]) / 2 : 0.5
    x_pad_right = length(val_ib) > 1 ? (val_ib[end] - val_ib[end-1]) / 2 : 0.5
    y_pad_low = length(freq_query) > 1 ? (freq_query[2] - freq_query[1]) / 2 : 0.5
    y_pad_high = length(freq_query) > 1 ? (freq_query[end] - freq_query[end-1]) / 2 : 0.5
    limits_ib_band = (first(val_ib) - x_pad_left, last(val_ib) + x_pad_right)
    limits_freq_band = (first(freq_query) - y_pad_low, last(freq_query) + y_pad_high)
    for ax in axs_band
        ax.xticks = WilkinsonTicks(5; k_min=3)
    end
    linkaxes!(axs_band...)
    xlims!(first(axs_band), limits_ib_band...)
    ylims!(first(axs_band), limits_freq_band...)
    on(button_reset_band.clicks) do _
        xlims!(first(axs_band), limits_ib_band...)
        ylims!(first(axs_band), limits_freq_band...)
    end

    function update_evolution!()
        data = get_data(obs_property[])
        for i in 1:2
            obs_evol[i][] = get_evolution(data, obs_ib[], obs_rep[], i)
            for r in 1:n_rep
                obs_evol_reps[r, i][] = get_evolution(data, obs_ib[], r, i)
            end
        end
        spec, variant = specs[obs_property[]]
        ax_evol.ylabel = spec.ylabel
        ax_evol.title = "$(spec.name) / $(variant.label)"
        axis_settings = get_evolution_axis(obs_property[])
        ax_evol.yticks = axis_settings.ticks
        ylims!(ax_evol, axis_settings.limits...)
    end
    function update_profile!()
        kind, data, pos = get_profile(obs_property[])
        name_profile = Dict(
            :axial => "axial core (μm)", :radial => "radial core (μm)",
            :axial_norm => "axial core norm (μm)", :radial_norm => "radial core norm (μm)",
            :modl => "modl (μm⁻¹)",
        )[kind]
        obs_profile_ylabel[] = name_profile
        limits = kind === :modl && !isnothing(cache.config_corr) &&
            hasproperty(cache.config_corr, :vis_evol_prfl_modl) ?
            cache.config_corr.vis_evol_prfl_modl.ylims : extrema(pos)
        for i in 1:2
            profile = obs_rep[] <= n_rep ? data.reps[obs_ib[], obs_rep[], i] : data.stack[obs_ib[], i]
            size(profile) == (length(pos), length(val_t)) || throw(DimensionMismatch(
                "$(kind) profile size $(size(profile)) must be $((length(pos), length(val_t)))"))
            isnothing(profile_plots[i]) || delete!(axs_profile[i], profile_plots[i])
            profile_plots[i] = heatmap!(axs_profile[i], val_t, pos, profile';
                colormap=gen_clrmap_solo(hue_theme_istp[val_istp[i]]))
        end
        ylims!(first(axs_profile), limits...)
    end
    function update_band!()
        data = get_data(obs_property[])
        band = calc_ib_band(data, val_t, freq_query, obs_time_range[], band_modes[obs_band_mode[]]; query_weight_kwargs)
        for i in 1:2
            obs_band_weight[i][] = band.weight[:, :, i]
            obs_band_phase_color[i][] = color_ib_band_phase(band.weight[:, :, i], band.phase[:, :, i])
        end
    end
    onany(obs_property, obs_ib, obs_rep) do _...
        update_evolution!()
        update_profile!()
    end
    onany(obs_property, obs_time_range, obs_band_mode) do _...
        update_band!()
    end
    update_evolution!()
    update_profile!()
    xlims!(ax_evol, first(val_t), last(val_t))
    rowsize!(fig.layout, 2, Fixed(535))
    display(fig)
    return fig
end

fig_ib_band = run_ib_band_sel()
