## Prepare selected SSNTFR repeat collections for later analysis.
using GLMakie
using HDF5
using ImageFiltering
using LsqFit: curve_fit, stderror
using Printf
using Statistics

GLMakie.activate!()

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
path_data = joinpath(path_root, "0204_interference", "result", "data.h5")
path_output = joinpath(path_root, "AnlzRoutine", "36.PhaseDistro")

tag = "SSNTFR"
val_istp = ["162", "164"]
val_IB_ref = [
    5.310,
    5.312,
    5.314,
    5.316,
    5.317,
    5.318,
    5.319,
    5.320,
    5.322,
    5.324,
    5.326,
    5.328,
    5.330,
    5.332,
    5.334,
    5.338,
    5.342,
]
smwh = (150, 150)
mag = 22.06
pixsz = 6.5
bin = 1
num_err = 0.6e4
range_center = 181:220
sigma_center_filter = 5
x_max_fit = 10 # μm
sigma_peak_init = 12.0
eta_peak_init = 0.5
lambda_peak_init = 5.0
phi_peak_init = 0.0
# [amp_init, sigma_peak_init, eta_peak_init, lambda_peak_init, phi_peak_init]
fit_lower_peak = [ 0.0, 10.0, 0.0, 3.0, -2pi]
fit_upper_peak = [25.0, 25.0, 1.5, 6.0, 2pi]

# live inspector selections
ib, istp, idx_rep = (5, 1, 1)
y_row = 0.0
x_col = 0.0
ylims_profile = (-1.0, 15.0)

function modl_peak_1d(x, p)
    (A, σ, η, λ, φ) = p
    @. A * exp(-(x/σ)^2) * (1 + η * cos(2π * x/λ - φ))
end

function load_density_payload(path_data::AbstractString, val_istp::AbstractVector{<:AbstractString})
    name_dataset_by_istp = Dict(
        "162" => "im64us",
        "164" => "im62us",
    )

    h5open(path_data, "r") do file
        dens_loaded = map(val_istp) do istp
            read(file[name_dataset_by_istp[istp]])
        end
        _, _, n_rep, n_IB = size(first(dens_loaded))
        dens_raw = Array{Matrix{Float64}}(undef, n_IB, length(val_istp), n_rep)
        for idx_IB in 1:n_IB, idx_istp in eachindex(val_istp), idx_rep in 1:n_rep
            dens_raw[idx_IB, idx_istp, idx_rep] = Float64.(copy(@view dens_loaded[idx_istp][:, :, idx_rep, idx_IB]))
        end
        return dens_raw
    end
end

function gaussian_offset_1d(x, p)
    return @. p[1] * exp(-((x - p[2])^2) / (2 * p[3]^2)) + p[4]
end

function draw_profile_inspector!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    dens_core::AbstractMatrix,
    fit_peak::AbstractMatrix,
    val_istp::AbstractVector;
    ib::Integer,
    istp::Integer,
    idx_rep::Integer,
    y_row::Real,
    x_col::Real,
    smidx_mean_profile::Integer,
    ylims_profile::Tuple{<:Real,<:Real},
)
    ib in axes(dens_core, 1) || throw(ArgumentError("ib must be in $(axes(dens_core, 1)), got $ib."))
    istp in axes(dens_core, 2) || throw(ArgumentError("istp must be in $(axes(dens_core, 2)), got $istp."))
    size(dens_core, 2) == length(val_istp) || throw(DimensionMismatch(
        "dens_core second dimension $(size(dens_core, 2)) must match length(val_istp) $(length(val_istp)).",
    ))
    size(fit_peak) == size(dens_core) || throw(DimensionMismatch(
        "fit_peak size $(size(fit_peak)) must match dens_core size $(size(dens_core)).",
    ))
    for idx in CartesianIndices(dens_core)
        isempty(dens_core[idx]) && continue
        length(fit_peak[idx]) == length(dens_core[idx]) || throw(DimensionMismatch(
            "fit_peak[$(Tuple(idx)...)] length $(length(fit_peak[idx])) must match dens_core length $(length(dens_core[idx])).",
        ))
        size(first(dens_core[idx])) == (length(x_dens), length(x_dens)) || throw(DimensionMismatch(
            "dens_core[$(Tuple(idx)...)] crop size $(size(first(dens_core[idx]))) must match " *
            "(length(x_dens), length(x_dens)) $((length(x_dens), length(x_dens))).",
        ))
    end

    dens_vec = dens_core[ib, istp]
    isempty(dens_vec) && throw(ArgumentError("dens_core[$ib, $istp] has no selected crops."))
    idx_rep = mod1(idx_rep, length(dens_vec))
    idx_row = argmin(abs.(x_dens .- y_row))
    idx_col = argmin(abs.(x_dens .- x_col))
    idx_center = cld(length(x_dens), 2)
    idxs_center = max(1, idx_center - smidx_mean_profile):min(length(x_dens), idx_center + smidx_mean_profile)
    dens2d = dens_vec[idx_rep]
    fit_info = fit_peak[ib, istp][idx_rep]
    dens_mean = mean(dens_vec)

    gen_theme_clr(idx_istp::Integer, alpha::Real) =
        RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
    gen_theme_clrmap(idx_istp::Integer) =
        gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)

    clr_mean = RGBAf(0.35, 0.35, 0.35, 0.62)
    clr_strip = RGBAf(0.86, 0.86, 0.86, 0.22)
    clr_fit = RGBAf(Oklch(0.60, 0.17, 145), 0.95)
    step_dens = median(diff(x_dens))
    x_strip_min = x_dens[first(idxs_center)] - step_dens / 2
    x_strip_max = x_dens[last(idxs_center)] + step_dens / 2
    y_strip_min, y_strip_max = x_strip_min, x_strip_max

    obs_idx_IB = Observable(ib)
    obs_idx_istp = Observable(istp)
    obs_idx_rep = Observable(idx_rep)
    obs_idx_row = Observable(idx_row)
    obs_idx_col = Observable(idx_col)
    obs_val_row = Observable(x_dens[idx_row])
    obs_val_col = Observable(x_dens[idx_col])
    obs_dens2d = Observable(dens2d)
    obs_dens2d_hm = lift(ds -> ds', obs_dens2d)
    obs_colorrange = Observable((0.0, maximum(dens2d)))
    obs_clrmap = Observable(gen_theme_clrmap(istp))
    obs_clr_theme = Observable(gen_theme_clr(istp, 0.3))
    obs_clr_theme_faint = Observable(gen_theme_clr(istp, 0.40))
    obs_profile_row = Observable(vec(@view dens2d[idx_row, :]))
    obs_profile_col = Observable(vec(@view dens2d[:, idx_col]))
    obs_profile_row_mean = Observable(vec(mean(@view(dens2d[idxs_center, :]); dims=1)))
    obs_profile_col_mean = Observable(vec(mean(@view(dens2d[:, idxs_center]); dims=2)))
    obs_fit_row = Observable(fit_info.fit)
    obs_fit_text = Observable(
        @sprintf(
            "A=%.3g\nσ=%.3g\nη=%.3g\nλ=%.3g\nφ=%.3g",
            fit_info.params...,
        ),
    )
    obs_title = lift(obs_idx_IB, obs_idx_istp, obs_idx_rep, obs_val_row, obs_val_col) do idx_IB_live, idx_istp_live, idx_rep_live, val_row_live, val_col_live
        @sprintf(
            "IB idx=%d, istp=%s, rep=%d/%d, y_row=%.3f μm, x_col=%.3f μm",
            idx_IB_live,
            string(val_istp[idx_istp_live]),
            idx_rep_live,
            length(dens_core[idx_IB_live, idx_istp_live]),
            val_row_live,
            val_col_live,
        )
    end
    obs_title_row = lift(obs_val_row) do val_row_live
        @sprintf("y_row=%.3f μm", val_row_live)
    end
    obs_title_col = lift(obs_val_col) do val_col_live
        @sprintf("x_col=%.3f μm", val_col_live)
    end

    Label(fig[0, 1:2]; text=obs_title, tellwidth=false, halign=:left)

    ax_hm = Axis(
        fig[1, 1];
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect(),
        xgridvisible=true,
        ygridvisible=true,
    )
    try
        deregister_interaction!(ax_hm, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    hspan!(ax_hm, y_strip_min, y_strip_max; color=clr_strip)
    vspan!(ax_hm, x_strip_min, x_strip_max; color=clr_strip)
    hm = heatmap!(ax_hm, x_dens, x_dens, obs_dens2d_hm; colormap=obs_clrmap, colorrange=obs_colorrange, rasterize=true)
    hlines!(ax_hm, lift(x -> [x], obs_val_row); color=obs_clr_theme, linewidth=0.9)
    vlines!(ax_hm, lift(x -> [x], obs_val_col); color=obs_clr_theme_faint, linewidth=0.9)

    ax_row = Axis(
        fig[2, 1];
        xlabel="x (μm)",
        ylabel="density",
        title=obs_title_row,
    )
    try
        deregister_interaction!(ax_row, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    lines!(ax_row, x_dens, obs_profile_row_mean; color=clr_mean, linewidth=2.5)
    lines!(ax_row, x_dens, obs_profile_row; color=obs_clr_theme, linewidth=1.7)
    lines!(ax_row, x_dens, obs_fit_row; color=clr_fit, linewidth=1.0)
    text!(
        ax_row,
        0.98,
        0.96;
        text=obs_fit_text,
        space=:relative,
        align=(:right, :top),
        color=clr_fit,
        fontsize=10,
    )
    xlims!(ax_row, extrema(x_dens))
    ylims!(ax_row, ylims_profile)

    ax_col = Axis(
        fig[1, 2];
        xlabel="density",
        ylabel="y (μm)",
        title=obs_title_col,
    )
    try
        deregister_interaction!(ax_col, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    lines!(ax_col, obs_profile_col_mean, x_dens; color=clr_mean, linewidth=2.5)
    lines!(ax_col, obs_profile_col, x_dens; color=obs_clr_theme, linewidth=1.7)
    xlims!(ax_col, ylims_profile)
    ylims!(ax_col, extrema(x_dens))
    ax_col.xreversed = true

    function update_profiles!()
        dens_vec_live = dens_core[obs_idx_IB[], obs_idx_istp[]]
        isempty(dens_vec_live) && return nothing
        obs_idx_rep[] = mod1(obs_idx_rep[], length(dens_vec_live))
        dens2d_live = dens_vec_live[obs_idx_rep[]]
        fit_info_live = fit_peak[obs_idx_IB[], obs_idx_istp[]][obs_idx_rep[]]
        dens_mean_live = mean(dens_vec_live)
        obs_dens2d[] = dens2d_live
        obs_colorrange[] = (0.0, maximum(dens2d_live))
        obs_clrmap[] = gen_theme_clrmap(obs_idx_istp[])
        obs_clr_theme[] = gen_theme_clr(obs_idx_istp[], 0.3)
        obs_clr_theme_faint[] = gen_theme_clr(obs_idx_istp[], 0.70)
        obs_profile_col[] = vec(@view dens2d_live[:, obs_idx_col[]])
        obs_profile_row[] = vec(@view dens2d_live[obs_idx_row[], :])
        obs_profile_row_mean[] = vec(mean(@view(dens2d_live[idxs_center, :]); dims=1))
        obs_profile_col_mean[] = vec(mean(@view(dens2d_live[:, idxs_center]); dims=2))
        obs_fit_row[] = fit_info_live.fit
        obs_fit_text[] = @sprintf(
            "A=%.3g\nσ=%.3g\nη=%.3g\nλ=%.3g\nφ=%.3g",
            fit_info_live.params...,
        )
        return nothing
    end

    function update_cut_profiles!(x_click::Real, y_click::Real)
        idx_col_live = argmin(abs.(x_dens .- x_click))
        idx_row_live = argmin(abs.(x_dens .- y_click))
        obs_idx_col[] = idx_col_live
        obs_idx_row[] = idx_row_live
        obs_val_col[] = x_dens[idx_col_live]
        obs_val_row[] = x_dens[idx_row_live]
        update_profiles!()
        return nothing
    end

    function update_data_index!(step_IB::Integer, step_istp::Integer, step_profile::Integer)
        obs_idx_IB[] = mod1(obs_idx_IB[] + step_IB, size(dens_core, 1))
        obs_idx_istp[] = mod1(obs_idx_istp[] + step_istp, size(dens_core, 2))
        dens_vec_live = dens_core[obs_idx_IB[], obs_idx_istp[]]
        isempty(dens_vec_live) && return nothing
        obs_idx_rep[] = mod1(obs_idx_rep[] + step_profile, length(dens_vec_live))
        update_profiles!()
        return nothing
    end

    click_handler = on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press && is_mouseinside(ax_hm.scene)
            xy_click = mouseposition(ax_hm)
            update_cut_profiles!(xy_click[1], xy_click[2])
        end
        return Consume(false)
    end

    gl_ctrl = GridLayout(fig[2, 2])
    labels = ("IB", "istp", "rep")
    steps = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    button_handlers = map(enumerate(labels)) do (idx_ctrl, label_ctrl)
        step = steps[idx_ctrl]
        btn_prev = Button(gl_ctrl[idx_ctrl, 1]; label="←", width=34, height=30)
        Label(gl_ctrl[idx_ctrl, 2]; text=label_ctrl, tellwidth=true, tellheight=false, halign=:center, valign=:center)
        btn_next = Button(gl_ctrl[idx_ctrl, 3]; label="→", width=34, height=30)
        (
            on(btn_prev.clicks) do _
                update_data_index!((-step[1]), (-step[2]), (-step[3]))
            end,
            on(btn_next.clicks) do _
                update_data_index!(step...)
            end,
        )
    end

    colsize!(fig.layout, 1, Fixed(360))
    colsize!(fig.layout, 2, Fixed(300))
    rowsize!(fig.layout, 1, Fixed(360))
    rowsize!(fig.layout, 2, Fixed(260))
    resize_to_layout!(fig)
    return (;
        ax_hm,
        ax_row,
        ax_col,
        hm,
        idx_IB=obs_idx_IB,
        idx_istp=obs_idx_istp,
        idx_rep=obs_idx_rep,
        idx_row=obs_idx_row,
        idx_col=obs_idx_col,
        y_row=obs_val_row,
        x_col=obs_val_col,
        click_handler,
        button_handlers,
    )
end

println("  [$tag] loading densities from $path_data")
dens_raw_fmt = load_density_payload(path_data, val_istp)
n_IB, n_istp, n_rep = size(dens_raw_fmt)
wh_raw = size(dens_raw_fmt[1, 1, 1])
println("  [$tag] formatted densities as (IB, istp, rep)=$(size(dens_raw_fmt)), image size=$wh_raw")
length(val_IB_ref) == n_IB || throw(DimensionMismatch("val_IB_ref length $(length(val_IB_ref)) must match IB count $n_IB."))
length(val_istp) == n_istp || throw(DimensionMismatch("val_istp length $(length(val_istp)) must match istp count $n_istp."))

cfg_prfl = get_prfl_modl_1d_config(smwh)
x_dens = (pixsz * bin / mag) .* collect(-smwh[2]:smwh[2])
val_IB = copy(val_IB_ref)

num = Array{Float64}(undef, n_IB, n_istp, n_rep)
xy_center = Array{Tuple{Int,Int}}(undef, n_IB, n_istp, n_rep)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp, idx_rep in 1:n_rep
    dens = dens_raw_fmt[idx_IB, idx_istp, idx_rep]
    dens_smooth = imfilter(dens, Kernel.gaussian(sigma_center_filter))

    prfl_x = vec(sum(dens_smooth; dims=1))
    x_fit = collect(1.0:length(prfl_x))
    p0_x = [maximum(prfl_x), (length(prfl_x) + 1) / 2, length(prfl_x) / 10, minimum(prfl_x)]
    x_center = curve_fit(gaussian_offset_1d, x_fit, Float64.(prfl_x), p0_x).param[2]

    prfl_y = vec(sum(dens_smooth; dims=2))
    y_fit = collect(1.0:length(prfl_y))
    p0_y = [maximum(prfl_y), (length(prfl_y) + 1) / 2, length(prfl_y) / 10, minimum(prfl_y)]
    y_center = curve_fit(gaussian_offset_1d, y_fit, Float64.(prfl_y), p0_y).param[2]

    num[idx_IB, idx_istp, idx_rep] = sum(dens)
    xy_center[idx_IB, idx_istp, idx_rep] = round.(Int, (x_center, y_center))
end

num_median = dropdims(median(num; dims=3); dims=3)
mask_valid_duet = falses(n_IB, n_rep)
for idx_IB in 1:n_IB, idx_rep in 1:n_rep
    is_valid_number = all(
        abs(num[idx_IB, idx_istp, idx_rep] - num_median[idx_IB, idx_istp]) <= num_err
        for idx_istp in 1:n_istp
    )
    is_valid_center = all(
        (xy -> xy[1] in range_center && xy[2] in range_center)(xy_center[idx_IB, idx_istp, idx_rep])
        for idx_istp in 1:n_istp
    )
    mask_valid_duet[idx_IB, idx_rep] = is_valid_number && is_valid_center
end

count_profile_shot = vec(sum(mask_valid_duet; dims=2))
println("  [$tag] valid duet counts per IB=$(count_profile_shot)")

ids_rep_valid = [findall(@view mask_valid_duet[idx_IB, :]) for idx_IB in 1:n_IB]
dens_core = Array{Vector{Matrix{Float64}}}(undef, n_IB, n_istp)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    dens_core[idx_IB, idx_istp] = [
        crop_center(dens_raw_fmt[idx_IB, idx_istp, idx_rep], xy_center[idx_IB, idx_istp, idx_rep], smwh) |> copy
        for idx_rep in 1:n_rep
        if mask_valid_duet[idx_IB, idx_rep]
    ]
end

idx_center = cld(length(x_dens), 2)
idxs_center = max(1, idx_center - cfg_prfl.smh_dens_strip):min(length(x_dens), idx_center + cfg_prfl.smh_dens_strip)
mask_fit = abs.(x_dens) .<= x_max_fit
x_fit_peak = x_dens[mask_fit]

fit_peak = Array{Vector{NamedTuple}}(undef, n_IB, n_istp)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    fit_peak[idx_IB, idx_istp] = map(enumerate(dens_core[idx_IB, idx_istp])) do (idx_rep_valid, dens2d)
        profile = vec(mean(@view(dens2d[idxs_center, :]); dims=1))
        prfl_strip_mean = Float64.(profile[mask_fit])
        amp_init = max(maximum(prfl_strip_mean), eps(Float64))
        p_init = [amp_init, sigma_peak_init, eta_peak_init, lambda_peak_init, phi_peak_init+pi]
        p_lower = copy(fit_lower_peak)
        p_upper = copy(fit_upper_peak)
        p_upper[1] = 2 * amp_init
        try
            fit = curve_fit(modl_peak_1d, x_fit_peak, prfl_strip_mean, p_init; lower=p_lower, upper=p_upper)
            param_err = try
                stderror(fit)
            catch err
                err isa SingularException || rethrow()
                fill(NaN, length(fit.param))
            end
            (;
                idx_rep=ids_rep_valid[idx_IB][idx_rep_valid],
                success=true,
                params=copy(fit.param),
                param_err,
                profile,
                fit=modl_peak_1d(x_dens, fit.param),
                resid=copy(fit.resid),
            )
        catch err
            @warn "modl_peak_1d fit failed" idx_IB idx_istp idx_rep=ids_rep_valid[idx_IB][idx_rep_valid] err
            (;
                idx_rep=ids_rep_valid[idx_IB][idx_rep_valid],
                success=false,
                params=fill(NaN, length(p_init)),
                param_err=fill(NaN, length(p_init)),
                profile,
                fit=fill(NaN, length(x_dens)),
                resid=fill(NaN, length(x_fit_peak)),
            )
        end
    end
end
count_fit = sum(sum(f.success for f in fits) for fits in fit_peak)
count_fit_err = sum(sum(f.success && any(isnan, f.param_err) for f in fits) for fits in fit_peak)
println("  [$tag] fitted modl_peak_1d for $count_fit selected crops; singular error estimates for $count_fit_err crops")

ntfr2d_mean = map(dens_core) do ds
    isempty(ds) && throw(ArgumentError("No valid densities available for a condition."))
    dropdims(mean(stack(ds); dims=3); dims=3)
end

isdir(path_output) || mkpath(path_output)

fig_live = Figure(fontsize=14)
profile_axes = draw_profile_inspector!(
    fig_live,
    x_dens,
    dens_core,
    fit_peak,
    val_istp;
    ib,
    istp,
    idx_rep,
    y_row,
    x_col,
    smidx_mean_profile=cfg_prfl.smh_dens_strip,
    ylims_profile,
)
display(fig_live)
