using Statistics: std

## 
if false
t_stage = log_step("fitting PCA modes")
selector_t_pca_dens = @isdefined(selector_t_pca_dens) ? selector_t_pca_dens : (@isdefined(selector_t_pca) ? selector_t_pca : (t -> trues(length(t))))
selector_t_pca_modl = @isdefined(selector_t_pca_modl) ? selector_t_pca_modl : selector_t_pca_dens
n_pca_modes_prfl_modl = @isdefined(n_pca_modes_prfl_modl) ? n_pca_modes_prfl_modl : 12
modes_pca_dens2d = [
    begin
        println("  [$tag] fitting PCA IB_idx=$c")
        flush(stdout)
        essn_2d_fmt[c, :, :, :] |> es -> map(a -> a.dens2d_core |> filter_core_pca, es) |> es -> eachslice(es; dims=(1, 2)) |> m -> fit_pca_modes(n_pca_modes, m)
    end
    for c in axes(essn_2d_fmt, 1)
]
pca_spectra = [[
    begin
        mode = modes_pca_dens2d[c][m]
        spectral_weight = calc_spct_rep_evol(eachslice(mode.weight; dims=1), val_vars.t_hold, freq_query_pca; sel_evol=selector_t_pca_dens)
        peaks_prominent = spectral_weight.spct_mean_mask |> spct -> get_spectrum_peaks(freq_query_pca, spct; min_prom=0.2)
        (; spectral_weight, peaks_prominent)
    end
    for m in 1:n_pca_modes
] for c in axes(modes_pca_dens2d, 1)
]
log_done("fit density PCA modes", t_stage)

t_stage = log_step("analyzing per-shot trends")
trend_sidepeak_nvlp = [
    extr_fmt[c, r, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_spectrum, query_weight_kwargs)
    for c in axes(extr_fmt, 1), r in axes(extr_fmt, 2), i in axes(extr_fmt, 4)
]
log_done("analyzed per-shot trends", t_stage)

t_stage = log_step("fitting trend evolution properties")
fit_evol_properties = fit_evol_properties_from_trends(trend_sidepeak_nvlp, trend_property_specs)
log_done("fit trend evolution properties", t_stage)

t_stage = log_step("analyzing stacked trends")
trend_extr_stacked_over_rep = [
    extr_stacked_over_rep[c, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_spectrum, query_weight_kwargs)
    for c in axes(extr_stacked_over_rep, 1), i in axes(extr_stacked_over_rep, 3)
]
trend_stacked_over_rep = [
    trend_sidepeak_nvlp[c, :, i] |> mean_dict
    for c in axes(trend_sidepeak_nvlp, 1), i in axes(trend_sidepeak_nvlp, 3)
]
log_done("analyzed stacked trends", t_stage)
end
##
t_stage = log_step("composing FT sidepeak profile evolution")
prfl_evol = [
    [
        extr_fmt[c, r, t, i].sidepeak.prfl_norm_tailess_px
        for t in axes(extr_fmt, 3)
    ] |> prfls -> reduce(hcat, prfls)
    for c in axes(extr_fmt, 1), r in axes(extr_fmt, 2), i in axes(extr_fmt, 4)
]
prfl_evol_stacked = nothing
log_done("finished composing FT sidepeak profile evolution", t_stage)

t_stage = log_step("fitting number evolution and building validity mask")
function twostep_decay_num(t::AbstractVector, p::AbstractVector)
    D1, λ1, D2, λ2 = p
    @. D1 * exp(-t / λ1) + D2 * exp(-t / λ2)
end

function fit_num_decay(t::AbstractVector{<:Real}, nums::AbstractVector{<:Real})
    ymax = maximum(nums)
    ymax > 0 || throw(ArgumentError("number values must contain a positive maximum"))
    p_init = [ymax / 2, 50.0, ymax / 2, 200.0]
    p_lower = [0.0, 1e-6, 0.0, 1e-6]
    p_upper = [ymax, Inf, ymax, Inf]
    coef(curve_fit(twostep_decay_num, Float64.(t), Float64.(nums), p_init; lower=p_lower, upper=p_upper))
end

rng_t_num_decay = @isdefined(rng_t_num_decay) ? rng_t_num_decay : (0.0, 150.0)
mask_t_num_decay = rng_t_num_decay[1] .<= val_vars.t_hold .<= rng_t_num_decay[2]
num_fit = map(extr_fmt) do extr
    p = extr.envelope.params_asymm
    2π * prod(p.size) * p.max / (px_in_um^2)
end

function make_num_fit(ib::Int, istp::Int)
    t_hold_num = Float64.(val_vars.t_hold)
    t_fit = t_hold_num[mask_t_num_decay]
    nums = [Float64.(vec(num_fit[ib, r, :, istp])) for r in axes(num_fit, 2)]
    nums_fit = [nums[r][mask_t_num_decay] for r in axes(num_fit, 2)]
    params = fit_num_decay(repeat(t_fit, n_rep), vcat(nums_fit...))
    fitted = [twostep_decay_num(t_hold_num, params) for _ in axes(num_fit, 2)]
    errors_rel = [(nums[r] .- fitted[r]) ./ fitted[r] for r in axes(num_fit, 2)]
    σ = std(vcat([errors_rel[r][mask_t_num_decay] for r in axes(num_fit, 2)]...))
    valid = reduce(vcat, (permutedims(abs.(errors_rel[r]) .<= 2σ) for r in axes(num_fit, 2)))
    errors = [nums[r] .- fitted[r] for r in axes(num_fit, 2)]
    return (; t=t_hold_num, nums, params, fitted, errors, errors_rel, σ, valid_num=valid)
end

fits_num = [make_num_fit(ib, istp) for ib in axes(num_fit, 1), istp in axes(num_fit, 4)]
valid_num = permutedims(
    reduce((a, b) -> cat(a, b; dims=4),
        [cat([fits_num[ib, istp].valid_num for istp in axes(num_fit, 4)]...; dims=3)
         for ib in axes(num_fit, 1)]),
    (4, 1, 2, 3),
)
log_done("fitted number evolution and built validity mask", t_stage)

prfl_evol_stacked = [
    [
        begin
            vals = [extr_fmt[c, r, t, i].sidepeak.prfl_norm_tailess_px
                    for r in axes(extr_fmt, 2) if valid_num[c, r, t, i]]
            isempty(vals) ?
                fill(missing, length(extr_fmt[c, firstindex(extr_fmt, 2), t, i].sidepeak.prfl_norm_tailess_px)) :
                vec(mean(reduce(hcat, vals); dims=2))
        end
        for t in axes(extr_fmt, 3)
    ] |> prfls -> reduce(hcat, prfls)
    for c in axes(extr_fmt, 1), i in axes(extr_fmt, 4)
]

t_stage = log_step("composing core density profile evolution")
compose_core_prfl_evol(essns, field::Symbol) = [
    getproperty(essns[t].prfls_core, field)
    for t in eachindex(essns)
] |> prfls -> reduce(hcat, prfls)
prfl_axial_evol = [
    compose_core_prfl_evol(essn_2d_fmt[c, r, :, i], :axial)
    for c in axes(essn_2d_fmt, 1), r in axes(essn_2d_fmt, 2), i in axes(essn_2d_fmt, 4)
]
prfl_axial_evol_stacked = [
    compose_core_prfl_evol(essn_stacked_over_rep[c, :, i], :axial)
    for c in axes(essn_stacked_over_rep, 1), i in axes(essn_stacked_over_rep, 3)
]
prfl_radial_evol = [
    compose_core_prfl_evol(essn_2d_fmt[c, r, :, i], :radial)
    for c in axes(essn_2d_fmt, 1), r in axes(essn_2d_fmt, 2), i in axes(essn_2d_fmt, 4)
]
prfl_radial_evol_stacked = [
    compose_core_prfl_evol(essn_stacked_over_rep[c, :, i], :radial)
    for c in axes(essn_stacked_over_rep, 1), i in axes(essn_stacked_over_rep, 3)
]

function normalize_core_prfl_evol(
    prfl_evol::AbstractMatrix,
    pos::AbstractVector{<:Real},
    thres_prfl_bot_mask::AbstractVector{<:Real},
    valid_t::AbstractVector{Bool},
)
    n_pos, n_t = size(prfl_evol)
    n_t >= 4 || throw(DimensionMismatch("core profile evolution needs at least 4 t_hold profiles, got $n_t"))
    length(pos) == n_pos || throw(DimensionMismatch("profile position length $(length(pos)) does not match profile size $n_pos"))
    length(thres_prfl_bot_mask) == n_t || throw(DimensionMismatch("profile threshold length $(length(thres_prfl_bot_mask)) does not match t_hold count $n_t"))
    length(valid_t) == n_t || throw(DimensionMismatch("validity length $(length(valid_t)) does not match t_hold count $n_t"))

    mask_com = [
        valid_t[t] && isfinite(prfl_evol[p, t]) && prfl_evol[p, t] > thres_prfl_bot_mask[t]
        for p in axes(prfl_evol, 1), t in axes(prfl_evol, 2)
    ]
    t_first = firstindex(prfl_evol, 2)
    total_first = sum(
        sum((prfl_evol[p, t] for p in axes(prfl_evol, 1) if mask_com[p, t]); init=0.0)
        for t in t_first:(t_first + 3); init=0.0
    )
    isfinite(total_first) && !iszero(total_first) || return Array{Union{Missing,Float64}}(missing, size(prfl_evol))

    step_pos = mean(abs, diff(pos))
    isfinite(step_pos) && !iszero(step_pos) || throw(ArgumentError("profile positions must have nonzero finite spacing"))
    prfl_norm = Array{Union{Missing,Float64}}(missing, size(prfl_evol))
    for t in axes(prfl_evol, 2)
        idx_mask = findall(@view(mask_com[:, t]))
        isempty(idx_mask) && continue
        total_mask = sum(prfl_evol[p, t] for p in idx_mask)
        isfinite(total_mask) && !iszero(total_mask) || continue
        scale = total_first / total_mask
        com = sum(pos[p] * prfl_evol[p, t] for p in idx_mask) / total_mask
        shift = round(Int, -com / step_pos)
        for p in axes(prfl_evol, 1)
            p_shift = p + shift
            checkbounds(Bool, prfl_norm, p_shift, t) && (prfl_norm[p_shift, t] = prfl_evol[p, t] * scale)
        end
    end
    return prfl_norm
end

function calc_top_prfl_core_thresholds(
    prfl_evol::AbstractArray{<:Any,3};
    quantile_mask_prfl::Real,
    valid_num::AbstractArray{Bool},
)
    0 < quantile_mask_prfl < 1 || throw(ArgumentError("quantile_mask_prfl must lie in (0, 1), got $quantile_mask_prfl"))
    n_t = size(first(prfl_evol), 2)
    thres_prfl_top = Array{Float64}(undef, size(prfl_evol, 1), n_t, size(prfl_evol, 3))
    for c in axes(prfl_evol, 1), t in axes(first(prfl_evol), 2), i in axes(prfl_evol, 3)
        vals = Float64[]
        for r in axes(prfl_evol, 2)
            valid_num[c, r, t, i] || continue
            append!(vals, filter(isfinite, vec(prfl_evol[c, r, i][:, t])))
        end
        thres_prfl_top[c, t, i] = isempty(vals) ? NaN : quantile(vals, 1 - quantile_mask_prfl)
    end
    return thres_prfl_top
end

function stack_core_prfl_evol(prfl_evol::AbstractArray{<:Any,3}, valid_num::AbstractArray{Bool})
    n_pos, n_t = size(first(prfl_evol))
    return [
        begin
            prfl_stack = Matrix{Union{Missing,Float64}}(missing, n_pos, n_t)
            for p in axes(prfl_stack, 1), t in axes(prfl_stack, 2)
                vals = [prfl_evol[c, r, i][p, t] for r in axes(prfl_evol, 2) if valid_num[c, r, t, i] && !ismissing(prfl_evol[c, r, i][p, t])]
                isempty(vals) && (vals = [prfl_evol[c, r, i][p, t] for r in axes(prfl_evol, 2) if !ismissing(prfl_evol[c, r, i][p, t])])
                isempty(vals) || (prfl_stack[p, t] = mean(vals))
            end
            prfl_stack
        end
        for c in axes(prfl_evol, 1), i in axes(prfl_evol, 3)
    ]
end

prfl_axial_evol_stacked = stack_core_prfl_evol(prfl_axial_evol, valid_num)
prfl_radial_evol_stacked = stack_core_prfl_evol(prfl_radial_evol, valid_num)

quantile_mask_prfl = @isdefined(quantile_mask_prfl) ? quantile_mask_prfl : 0.05
thres_frac_bot_mask_prfl = @isdefined(thres_frac_bot_mask_prfl) ? thres_frac_bot_mask_prfl : 0.1
essn_ref = first(essn_2d_fmt)
pos_axial = ((-essn_ref.smwh_core[2]):essn_ref.smwh_core[2]) .* essn_ref.step_posi[2]
pos_radial = ((-essn_ref.smwh_core[1]):essn_ref.smwh_core[1]) .* essn_ref.step_posi[1]
thres_prfl_top_axial = calc_top_prfl_core_thresholds(prfl_axial_evol; quantile_mask_prfl, valid_num)
thres_prfl_top_radial = calc_top_prfl_core_thresholds(prfl_radial_evol; quantile_mask_prfl, valid_num)
thres_prfl_bot_mask_axial = thres_frac_bot_mask_prfl .* thres_prfl_top_axial
thres_prfl_bot_mask_radial = thres_frac_bot_mask_prfl .* thres_prfl_top_radial
prfl_axial_evol_norm = [
    normalize_core_prfl_evol(prfl_axial_evol[c, r, i], pos_axial, @view(thres_prfl_bot_mask_axial[c, :, i]), @view(valid_num[c, r, :, i]))
    for c in axes(prfl_axial_evol, 1), r in axes(prfl_axial_evol, 2), i in axes(prfl_axial_evol, 3)
]
prfl_radial_evol_norm = [
    normalize_core_prfl_evol(prfl_radial_evol[c, r, i], pos_radial, @view(thres_prfl_bot_mask_radial[c, :, i]), @view(valid_num[c, r, :, i]))
    for c in axes(prfl_radial_evol, 1), r in axes(prfl_radial_evol, 2), i in axes(prfl_radial_evol, 3)
]
prfl_axial_evol_norm_stacked = stack_core_prfl_evol(prfl_axial_evol_norm, valid_num)
prfl_radial_evol_norm_stacked = stack_core_prfl_evol(prfl_radial_evol_norm, valid_num)
log_done("finished composing core density profile evolution", t_stage)

##
if false
t_stage = log_step("fitting modulation profile PCA modes")
mask_y_modl_pca = (0.06 .<= y_modl .<= 0.6)
any(mask_y_modl_pca) || throw(ArgumentError("modulation profile PCA wavenum selector 0.06-0.6 selected no y_modl values."))
y_modl_pca = y_modl[mask_y_modl_pca]
idx_y_modl_pca = findall(mask_y_modl_pca)
samples_pca_prfl_modl = [
    [
        prfl_evol[c, r, i][idx_y, t]
        for i in axes(prfl_evol, 3), idx_y in idx_y_modl_pca
    ]
    for c in axes(prfl_evol, 1), t in axes(first(prfl_evol), 2), r in axes(prfl_evol, 2)
]
modes_pca_prfl_modl = fit_pca_modes(n_pca_modes_prfl_modl, samples_pca_prfl_modl)
pca_spectra_prfl_modl = [
    begin
        mode = modes_pca_prfl_modl[m]
        spectral_weight = calc_spct_rep_evol(
            [vec(mode.weight[c, :, r]) for r in axes(mode.weight, 3)],
            val_vars.t_hold,
            freq_query_pca_modl;
            sel_evol=selector_t_pca_modl,
        )
        peaks_prominent = spectral_weight.spct_mean_mask |> spct -> get_spectrum_peaks(freq_query_pca_modl, spct; min_prom=0.2)
        (; spectral_weight, peaks_prominent)
    end
    for m in 1:n_pca_modes_prfl_modl, c in axes(samples_pca_prfl_modl, 1)
]
log_done("fit modulation profile PCA modes", t_stage)
end
##
config_corr = (;
    filter_core_pca_sigma,
    n_pca_modes_prfl_modl,
    freq_query_pca_modl,
    query_weight_kwargs,
    selector_t_pca_dens_val=val_vars.t_hold[selector_t_pca_dens(val_vars.t_hold)],
    selector_t_pca_modl_val=val_vars.t_hold[selector_t_pca_modl(val_vars.t_hold)],
    selector_t_spectrum_val=NamedTuple{propertynames(selector_t_spectrum)}(
        Tuple(selector_t_spectrum[key](val_vars.t_hold) |> mask -> val_vars.t_hold[mask] for key in propertynames(selector_t_spectrum))
    ),
    trend_property_specs,
    trend_panel_per_IB_kwargs,
    trend_panel_per_prop_kwargs,
    trend_all_IB_groups,
    trend_spectrum_IB_groups,
    trend_spectrum_IB_kwargs,
    trend_spectrum_IB_plot_kwargs,
    vis_evol_prfl_modl,
    vis_evol_prfl_axial,
    vis_evol_prfl_radial,
    quantile_mask_prfl,
    thres_frac_bot_mask_prfl,
    rng_t_num_decay,
)

meta_corr = merge(
    meta_extr |> m -> merge(m, (;fmt_dens=merge(m.fmt_dens, (;dens_full_fmt=nothing,)))),
    (;
        kind="excitation_corr",
        path_output,
        config_corr,
        # n_pca_modes_prfl_modl,
        # trend_property_specs,
        # y_modl_pca,
        # freq_query_pca_modl,
        # config_corr.selector_t_pca_dens_val,
        # config_corr.selector_t_pca_modl_val,
        # config_corr.selector_t_spectrum_val,
    ),
)

t_stage = log_step("saving excitation correlation cache")
path_cache_corr = joinpath(path_output, @sprintf("%s_corr.jld2", tag))
JLD2.jldsave(
    path_cache_corr;
    meta_corr,
    # trend_sidepeak_nvlp,
    # fit_evol_properties,
    # trend_extr_stacked_over_rep,
    # trend_stacked_over_rep,
    num_fit,
    fits_num,
    valid_num,
    prfl_evol,
    prfl_evol_stacked,
    prfl_axial_evol,
    prfl_axial_evol_stacked,
    prfl_radial_evol,
    prfl_radial_evol_stacked,
    prfl_axial_evol_norm,
    prfl_axial_evol_norm_stacked,
    prfl_radial_evol_norm,
    prfl_radial_evol_norm_stacked,
    # modes_pca_dens2d,
    # pca_spectra,
    # modes_pca_prfl_modl,
    # pca_spectra_prfl_modl,
)
log_done("saved excitation correlation cache", t_stage)
