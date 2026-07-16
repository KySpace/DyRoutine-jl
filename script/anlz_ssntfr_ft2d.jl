using GLMakie
using FFTW: fft, fftshift, ifftshift
using HDF5
using ImageFiltering
using ImageMorphology
using JLD2
using LinearAlgebra: norm
using LsqFit
using Printf
using Statistics
import CairoMakie

GLMakie.activate!()

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
path_data = joinpath(path_root, "0204_interference", "result", "data.h5")

# commit 57d2e69f9017be9957b38674e01e6fbc3aae013d
# Form FT1D profile after FT2D
# Adjusting the selected area
path_output = joinpath(path_root, "AnlzRoutine", "46.MeanAbsl2D.Mask")
path_fit_jld2 = joinpath(path_output, "SSNTFR_ft2d_fit.jld2")

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
smwh_dens_ft = (80, 80)
alpha_tukey = (1.0, 1.0)
mag = 22.06
pixsz = 6.5
bin = 1
sigma_center_filter = 5
use_common_xy_center = :free # :free, :fixed_x, or :fixed_xy

vis_cmpx_ampl_prescaler_power = 2
vis_cmpx_ampl_prescaler_scale = 0.4
vis_cmpx_ampl_prescaler = a -> a^vis_cmpx_ampl_prescaler_power * vis_cmpx_ampl_prescaler_scale
clrrng_ft2d_absl_mean = (0, 50)
clrrng_ft2d_cmpx_mean = (0, 30)
use_mask_sidepeak = true
ky_max_modl = 0.05

# reconstructed envelope and tail-removal settings
r_tail_min_profile = 20.0
fit_stride_2d = 3
fit_maxiter_2d = 1_000
fit_threshold_2d = 1.5e-1
fit_sigma_wide_min = 15.0
fit_r_narrow_max = r_tail_min_profile
kx_max_scale_reconstr = 0.1

ib, istp, idx_rep = (5, 1, 1)

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

function calc_modl_tail_masked(
    kx_ft::AbstractVector{<:Real},
    prfl::AbstractArray{<:Real,3},
    prfl_modl_fit::AbstractArray{<:Real,3},
)
    size(prfl) == size(prfl_modl_fit) || throw(DimensionMismatch(
        "profile size $(size(prfl)) must match prfl_modl_fit size $(size(prfl_modl_fit)).",
    ))
    length(kx_ft) == size(prfl, 1) || throw(DimensionMismatch(
        "kx_ft length $(length(kx_ft)) must match profile length $(size(prfl, 1)).",
    ))
    prfl_tailess = prfl .- prfl_modl_fit
    scale_reconstr = ones(Float64, size(prfl, 2), size(prfl, 3))
    prfl_total_avg = vec(mean(prfl; dims=(2, 3)))
    tail = vec(mean(prfl_modl_fit; dims=(2, 3)))
    tail_istp_avg = dropdims(mean(prfl_modl_fit; dims=3); dims=3)
    return (; prfl_modl_fit, prfl_modl_fit_scaled=prfl_modl_fit, scale_reconstr, tail, tail_istp_avg, prfl_total_avg, prfl_tailess)
end

function gaussian_offset_1d(x, p)
    return @. p[1] * exp(-((x - p[2])^2) / (2 * p[3]^2)) + p[4]
end

function shader_cmpx(ampl, phase; l_max=0.7438, c_max=0.1255, hue_offset=0, prescale=(t -> t), alpha_base=0.1, thres_alpha=0.05)
    l = @pipe ampl |> prescale |> clamp(_, 0, 1) |> 1 - _
    alpha = 1 - l |> u -> u > thres_alpha ? 1.0 : (u / thres_alpha * (1 - alpha_base) + alpha_base)
    c = l < l_max ? (l / l_max * c_max) : ((1 - l) / (1 - l_max) * c_max)
    hue = @pipe mod(phase, 2pi)/2pi*360 + hue_offset |> mod(_, 360)
    RGBAf(Oklch(l, c, hue), alpha)
end

function calc_modl_tail(
    kx_ft::AbstractVector{<:Real},
    prfl::AbstractArray{<:Real,3},
    prfl_modl_fit::AbstractArray{<:Real,3};
    kx_max_scale_reconstr::Real,
)
    size(prfl) == size(prfl_modl_fit) || throw(DimensionMismatch(
        "profile size $(size(prfl)) must match prfl_modl_fit size $(size(prfl_modl_fit)).",
    ))
    length(kx_ft) == size(prfl, 1) || throw(DimensionMismatch(
        "kx_ft length $(length(kx_ft)) must match profile length $(size(prfl, 1)).",
    ))
    sel_center = abs.(kx_ft) .<= kx_max_scale_reconstr
    any(sel_center) || throw(ArgumentError("No kx_ft values found within ±$kx_max_scale_reconstr."))
    prfl_tailess = similar(prfl, Float64)
    scale_reconstr = Array{Float64}(undef, size(prfl, 2), size(prfl, 3))
    for idx_IB in axes(prfl, 3), idx_istp in axes(prfl, 2)
        ids = (:, idx_istp, idx_IB)
        denom = sum(@view prfl_modl_fit[ids...][sel_center])
        denom > 0 || throw(ArgumentError("Reconstructed profile has nonpositive center sum at (IB=$idx_IB, istp=$idx_istp)."))
        scale = sum(@view prfl[ids...][sel_center]) / denom
        scale_reconstr[idx_istp, idx_IB] = scale
        prfl_tailess[ids...] .= @view(prfl[ids...]) .- scale .* @view(prfl_modl_fit[ids...])
    end
    prfl_modl_fit_scaled = similar(prfl_modl_fit, Float64)
    for idx_IB in axes(prfl, 3), idx_istp in axes(prfl, 2)
        prfl_modl_fit_scaled[:, idx_istp, idx_IB] .=
            scale_reconstr[idx_istp, idx_IB] .* @view(prfl_modl_fit[:, idx_istp, idx_IB])
    end
    prfl_total_avg = vec(mean(prfl; dims=(2, 3)))
    tail = vec(mean(prfl_modl_fit_scaled; dims=(2, 3)))
    tail_istp_avg = dropdims(mean(prfl_modl_fit_scaled; dims=3); dims=3)
    return (; prfl_modl_fit, prfl_modl_fit_scaled, scale_reconstr, tail, tail_istp_avg, prfl_total_avg, prfl_tailess)
end

##
println("  [$tag] loading densities from $path_data")
dens_raw_fmt = load_density_payload(path_data, val_istp)
n_IB, n_istp, n_rep = size(dens_raw_fmt)
wh_raw = size(dens_raw_fmt[1, 1, 1])
(x_center_px0, y_center_px0) = (wh_raw .+ 1) ./ 2
println("  [$tag] formatted densities as (IB, istp, rep)=$(size(dens_raw_fmt)), image size=$wh_raw")
length(val_IB_ref) == n_IB || throw(DimensionMismatch("val_IB_ref length $(length(val_IB_ref)) must match IB count $n_IB."))
length(val_istp) == n_istp || throw(DimensionMismatch("val_istp length $(length(val_istp)) must match istp count $n_istp."))

step_dens = pixsz * bin / mag
step_ft = 1 ./ (2 .* smwh_dens_ft .* step_dens)
x_dens = step_dens .* collect(-smwh[2]:smwh[2])
y_dens = step_dens .* collect(-smwh[1]:smwh[1])
val_IB = copy(val_IB_ref)
use_common_xy_center in (:free, :fixed_x, :fixed_xy) || throw(ArgumentError(
    "use_common_xy_center must be :free, :fixed_x, or :fixed_xy, got $use_common_xy_center.",
))

num = Array{Float64}(undef, n_IB, n_istp, n_rep)
xy_center_nvlp_px = Array{Tuple{Float64, Float64}}(undef, n_IB, n_istp, n_rep)
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
    xy_center_nvlp_px[idx_IB, idx_istp, idx_rep] = (x_center, y_center)
end
xy_center_shift = Array{Tuple{Int,Int}}(undef, n_IB, n_istp, n_rep)
if use_common_xy_center == :fixed_xy
    for idx_IB in 1:n_IB, idx_istp in 1:n_istp
        x_common = round(Int, mean(first.(xy_center_nvlp_px[idx_IB, idx_istp, :])))
        y_common = round(Int, mean(last.(xy_center_nvlp_px[idx_IB, idx_istp, :])))
        xy_center_shift[idx_IB, idx_istp, :] .= Ref((x_common, y_common))
    end
    println("  [$tag] using fixed xy_center repeated over reps for each IB, istp")
elseif use_common_xy_center == :fixed_x
    for idx_IB in 1:n_IB, idx_istp in 1:n_istp
        x_common = round(Int, mean(first.(xy_center_nvlp_px[idx_IB, idx_istp, :])))
        for idx_rep in 1:n_rep
            _, y_fit = xy_center_nvlp_px[idx_IB, idx_istp, idx_rep]
            xy_center_shift[idx_IB, idx_istp, idx_rep] = (x_common, round(Int, y_fit))
        end
    end
    println("  [$tag] using fixed x center and fitted y centers")
else
    xy_center_shift = map(c -> round.(Int, c), xy_center_nvlp_px)
    println("  [$tag] using fitted xy center")
end

xy_peak_nvlp = map(c -> (c .- x_center_px0) .* step_dens, xy_center_nvlp_px)

mask_valid_duet = trues(n_IB, n_rep)

count_profile_shot = vec(sum(mask_valid_duet; dims=2))
println("  [$tag] using all duet counts per IB=$(count_profile_shot)")

ids_rep_valid = [findall(@view mask_valid_duet[idx_IB, :]) for idx_IB in 1:n_IB]
dens_core = Array{Vector{Matrix{Float64}}}(undef, n_IB, n_istp)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    dens_core[idx_IB, idx_istp] = [
        crop_center(dens_raw_fmt[idx_IB, idx_istp, idx_rep], xy_center_shift[idx_IB, idx_istp, idx_rep], smwh) |> copy
        for idx_rep in 1:n_rep
        if mask_valid_duet[idx_IB, idx_rep]
    ]
end

##
gen_freq_shifted(n::Integer, step::Real) =
    iseven(n) ? collect((-n ÷ 2):(n ÷ 2 - 1)) ./ (n * step) :
    collect((-(n ÷ 2)):(n ÷ 2)) ./ (n * step)
idx_ft_x_center = argmin(abs.(x_dens))
idx_ft_y_center = argmin(abs.(y_dens))
idxs_ft_x = (idx_ft_x_center - smwh_dens_ft[1]):(idx_ft_x_center + smwh_dens_ft[1])
idxs_ft_y = (idx_ft_y_center - smwh_dens_ft[2]):(idx_ft_y_center + smwh_dens_ft[2])
first(idxs_ft_x) >= firstindex(x_dens) && last(idxs_ft_x) <= lastindex(x_dens) || throw(ArgumentError(
    "smwh_dens_ft=$smwh_dens_ft exceeds x_dens bounds around x=0.",
))
first(idxs_ft_y) >= firstindex(y_dens) && last(idxs_ft_y) <= lastindex(y_dens) || throw(ArgumentError(
    "smwh_dens_ft=$smwh_dens_ft exceeds y_dens bounds around y=0.",
))
kx_ft_full = gen_freq_shifted(length(idxs_ft_x), step_dens)
ky_ft_full = gen_freq_shifted(length(idxs_ft_y), step_dens)
mask_ft_kx = kx_ft_full .>= 0
mask_ft_ky = trues(length(ky_ft_full))

kx_ft = kx_ft_full[mask_ft_kx]
ky_ft = ky_ft_full[mask_ft_ky]
x_dens_ft = x_dens[idxs_ft_x]
y_dens_ft = y_dens[idxs_ft_y]
range_kx_ft_plot = (0.0, 0.5)
range_ky_ft_plot = (-0.2, 0.2)
step_kx_ft = median(diff(kx_ft))
step_ky_ft = median(diff(ky_ft))
mask_kx_ft_plot = (kx_ft .>= first(range_kx_ft_plot) - step_kx_ft / 2) .& (kx_ft .<= last(range_kx_ft_plot) + step_kx_ft / 2)
mask_ky_ft_plot = (ky_ft .>= first(range_ky_ft_plot) - step_ky_ft / 2) .& (ky_ft .<= last(range_ky_ft_plot) + step_ky_ft / 2)
kx_ft_plot = kx_ft[mask_kx_ft_plot]
ky_ft_plot = ky_ft[mask_ky_ft_plot]
lims_kx_ft_plot = (first(kx_ft_plot) - step_kx_ft / 2, last(kx_ft_plot) + step_kx_ft / 2)
lims_ky_ft_plot = (first(ky_ft_plot) - step_ky_ft / 2, last(ky_ft_plot) + step_ky_ft / 2)
mask_modl_ft_ky = (ky_max_modl .>= ky_ft .>= -ky_max_modl)
tukey_dens_ft = tukey1d(smwh_dens_ft[2]; alpha=alpha_tukey[2]) * tukey1d(smwh_dens_ft[1]; alpha=alpha_tukey[1])'


ntfr2d_mean = map(dens_core) do ds
    isempty(ds) && throw(ArgumentError("No valid densities available for a condition."))
    dropdims(mean(stack(ds); dims=3); dims=3)
end

fit_envelope = map(ntfr2d_mean) do dens2d
    guess_1d = fit_two_gaussian_1d_guess(x_dens, dens2d, r_tail_min_profile, :gaussian)
    fit_two_gaussian_2d(
        x_dens,
        dens2d,
        guess_1d;
        center_bound=0.0,
        stride=fit_stride_2d,
        threshold=fit_threshold_2d,
        sigma_wide_min=fit_sigma_wide_min,
        r_narrow_max=fit_r_narrow_max,
        maxiter=fit_maxiter_2d,
        model_center=:gaussian,
    )
end
coords_envelope = reduce(hcat, ([x, y] for y in y_dens, x in x_dens))
ntfr2d_reconstr = map(fit_envelope) do fit
    reshape(
        double_gaussian_disk_2d_model_abrr(coords_envelope, fit.params),
        length(y_dens),
        length(x_dens),
    )
end

function calc_ft2d_cmpx(
    dens2d::AbstractMatrix{<:Real},
    idxs_ft_y,
    idxs_ft_x,
    tukey_dens_ft,
    mask_ft_ky,
    mask_ft_kx,
    step_ft,
)
    @pipe dens2d[idxs_ft_y, idxs_ft_x] .* tukey_dens_ft |>
        ifftshift |> fft |> fftshift |>
        _[mask_ft_ky, mask_ft_kx] |>
        ft2d -> ft2d ./ (sum(abs.(ft2d)) .* prod(step_ft) ./ 4)
end

function calc_prfl_modl1d(
    ft2d_absl::AbstractMatrix{<:Real},
    mask_modl_ft_ky;
    mask_sidepeak=nothing,
)
    isnothing(mask_sidepeak) || size(ft2d_absl) == size(mask_sidepeak) || throw(DimensionMismatch(
        "ft2d size $(size(ft2d_absl)) must match mask size $(size(mask_sidepeak)).",
    ))
    ft2d_selected = isnothing(mask_sidepeak) ? ft2d_absl : ft2d_absl .* mask_sidepeak
    prfl_selected = vec(sum(@view(ft2d_selected[mask_modl_ft_ky, :]); dims=1))
    prfl_unmasked = vec(sum(@view(ft2d_absl[mask_modl_ft_ky, :]); dims=1))
    norm_unmasked = sum(prfl_unmasked[2:end]) + prfl_unmasked[1] / 2
    return prfl_selected ./ norm_unmasked
end

ft2d_cmpx = map(dens_core) do dens_reps
    map(dens_reps) do dens2d
        calc_ft2d_cmpx(dens2d, idxs_ft_y, idxs_ft_x, tukey_dens_ft, mask_ft_ky, mask_ft_kx, step_ft)
    end
end

ft2d_cmpx_mean = map(ft2d_cmpx) do ft2d_cmpx_reps
    reduce(+, ft2d_cmpx_reps) ./ length(ft2d_cmpx_reps)
end
ft2d_absl_mean = map(ft2d_cmpx) do ft2d_cmpx_reps
    reduce(+, map(ft -> abs.(ft), ft2d_cmpx_reps)) ./ length(ft2d_cmpx_reps)
end

mask_sidepeak = begin
    arg_seed = argmin([hypot(y, x .- 0) for y in ky_ft, x in kx_ft])
    ft_mean = reduce(+, vec(ft2d_absl_mean)) ./ length(ft2d_absl_mean)
    mask = ft_mean |> ft -> ft .>= (0.023 * ft[arg_seed])
    labels = label_components(mask, Bool[
    0 1 0
    1 1 1
    0 1 0
    ])
    labels .!= labels[arg_seed]
end

prfl_modl1d_cohr_unmasked = map(ft2d_cmpx_mean) do ft2d
    calc_prfl_modl1d(abs.(ft2d), mask_modl_ft_ky)
end
prfl_modl1d_inco_unmasked = map(ft2d_absl_mean) do ft2d
    calc_prfl_modl1d(ft2d, mask_modl_ft_ky)
end
prfl_modl1d_cohr_masked = map(ft2d_cmpx_mean) do ft2d
    calc_prfl_modl1d(abs.(ft2d), mask_modl_ft_ky; mask_sidepeak=mask_sidepeak)
end
prfl_modl1d_inco_masked = map(ft2d_absl_mean) do ft2d
    calc_prfl_modl1d(ft2d, mask_modl_ft_ky; mask_sidepeak=mask_sidepeak)
end
prfl_modl1d_cohr = use_mask_sidepeak ? prfl_modl1d_cohr_masked : prfl_modl1d_cohr_unmasked
prfl_modl1d_inco = use_mask_sidepeak ? prfl_modl1d_inco_masked : prfl_modl1d_inco_unmasked

ft2d_reconstr = map(ntfr2d_reconstr) do dens2d
    calc_ft2d_cmpx(dens2d, idxs_ft_y, idxs_ft_x, tukey_dens_ft, mask_ft_ky, mask_ft_kx, step_ft)
end
prfl_modl_fit = map(ft2d_reconstr) do ft2d
    calc_prfl_modl1d(abs.(ft2d), mask_modl_ft_ky; mask_sidepeak=use_mask_sidepeak ? mask_sidepeak : nothing)
end
prfl_modl_fit_fmt = Array{Float64}(undef, length(kx_ft), n_istp, n_IB)
prfl_modl1d_cohr_fmt = similar(prfl_modl_fit_fmt)
prfl_modl1d_inco_fmt = similar(prfl_modl_fit_fmt)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    prfl_modl_fit_fmt[:, idx_istp, idx_IB] .= prfl_modl_fit[idx_IB, idx_istp]
    prfl_modl1d_cohr_fmt[:, idx_istp, idx_IB] .= prfl_modl1d_cohr[idx_IB, idx_istp]
    prfl_modl1d_inco_fmt[:, idx_istp, idx_IB] .= prfl_modl1d_inco[idx_IB, idx_istp]
end
if use_mask_sidepeak
    tail_cohr = calc_modl_tail_masked(kx_ft, prfl_modl1d_cohr_fmt, prfl_modl_fit_fmt)
    tail_inco = calc_modl_tail_masked(kx_ft, prfl_modl1d_inco_fmt, prfl_modl_fit_fmt)
else
    tail_cohr = calc_modl_tail(kx_ft, prfl_modl1d_cohr_fmt, prfl_modl_fit_fmt; kx_max_scale_reconstr)
    tail_inco = calc_modl_tail(kx_ft, prfl_modl1d_inco_fmt, prfl_modl_fit_fmt; kx_max_scale_reconstr)
end
prfl_modl1d_cohr_tailess = [vec(tail_cohr.prfl_tailess[:, idx_istp, idx_IB]) for idx_IB in 1:n_IB, idx_istp in 1:n_istp]
prfl_modl1d_inco_tailess = [vec(tail_inco.prfl_tailess[:, idx_istp, idx_IB]) for idx_IB in 1:n_IB, idx_istp in 1:n_istp]
prfl_modl_fit_scaled_cohr = [vec(tail_cohr.prfl_modl_fit_scaled[:, idx_istp, idx_IB]) for idx_IB in 1:n_IB, idx_istp in 1:n_istp]
prfl_modl_fit_scaled_inco = [vec(tail_inco.prfl_modl_fit_scaled[:, idx_istp, idx_IB]) for idx_IB in 1:n_IB, idx_istp in 1:n_istp]

## Data save and FT-only visualizations

isdir(path_output) || mkpath(path_output)
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)

ylims_profile_cohr = (-0.005, 0.045)
ylims_profile_inco = (-0.005, 0.105)
fit_config = (;
    tag, path_data, path_output, smwh, smwh_dens_ft, alpha_tukey, mag, pixsz, bin,
    sigma_center_filter, use_common_xy_center, clrrng_ft2d_cmpx_mean,
    clrrng_ft2d_absl_mean, ylims_profile_cohr, ylims_profile_inco,
    vis_cmpx_ampl_prescaler_power, vis_cmpx_ampl_prescaler_scale,
    x_center_px0, step_dens, step_ft, r_tail_min_profile, fit_stride_2d,
    fit_maxiter_2d, fit_threshold_2d, fit_sigma_wide_min, fit_r_narrow_max,
    kx_max_scale_reconstr, use_mask_sidepeak, ky_max_modl,
    range_kx_ft_plot, range_ky_ft_plot, lims_kx_ft_plot, lims_ky_ft_plot,
)
JLD2.@save path_fit_jld2 fit_config x_dens y_dens x_dens_ft y_dens_ft kx_ft ky_ft kx_ft_plot ky_ft_plot val_IB val_istp num xy_center_nvlp_px xy_center_shift mask_valid_duet ids_rep_valid dens_core ntfr2d_mean fit_envelope ntfr2d_reconstr mask_sidepeak ft2d_cmpx ft2d_cmpx_mean ft2d_absl_mean prfl_modl_fit prfl_modl1d_cohr_unmasked prfl_modl1d_inco_unmasked prfl_modl1d_cohr_masked prfl_modl1d_inco_masked prfl_modl1d_cohr prfl_modl1d_inco prfl_modl1d_cohr_tailess prfl_modl1d_inco_tailess prfl_modl_fit_scaled_cohr prfl_modl_fit_scaled_inco
println("  [$tag] saved FT data to $path_fit_jld2")

function to_masked_clr(
    dens::AbstractMatrix{<:Real},
    mask::AbstractMatrix{Bool},
    hue;
    sat_max=0.24,
    max=16,
    thres_alpha=0.1,
    l_max=1.0,
    l_min=0.0,
    alpha_base=0.1,
)
    size(dens) == size(mask) || throw(DimensionMismatch("dens size $(size(dens)) does not match mask size $(size(mask))."))
    dens_norm = clamp.(dens, 0, max) ./ max
    alpha = (n, m) -> m ? (thres_alpha <= 0 ? (n > 0 ? 1.0 : alpha_base) : (n > thres_alpha ? 1.0 : (n / thres_alpha * (1 - alpha_base) + alpha_base))) : 0.0
    shader = (n, m) -> Oklch(l_max - (l_max - l_min) * abs(n), sat_max * abs(n), hue) |> c -> RGBAf(c, alpha(n, m))
    return [shader(dens_norm[x, y], mask[x, y]) for x in axes(dens, 1), y in axes(dens, 2)]
end

gen_amp_masked_clr(dens, idx_istp; max=clrrng_ft2d_absl_mean[2]) = (
    to_masked_clr(dens, .!mask_sidepeak, 0; sat_max=0.0, max),
    to_masked_clr(dens, mask_sidepeak, hue_theme_istp[string(val_istp[idx_istp])]; max),
)

clr_theme(idx_istp; alpha=1.0) = RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
clrmap_theme(idx_istp) = gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)
gen_ft_rgba(ft; rng_amp=clrrng_ft2d_cmpx_mean) = begin
    amp = abs.(ft)
    phs = angle.(ft)
    amp_max = isnothing(rng_amp) ? max(maximum(vec(amp[mask_sidepeak])), eps(Float64)) : rng_amp[2] / 2
    phase = shader_cmpx.(amp ./ amp_max, phs; prescale=vis_cmpx_ampl_prescaler)
    phase_axis = collect(range(-pi, pi; length=256))
    amp_axis = isnothing(rng_amp) ? collect(range(0, 2amp_max; length=256)) : collect(range(rng_amp...; length=256))
    legend = [shader_cmpx.(a ./ amp_max, p; prescale=vis_cmpx_ampl_prescaler) for a in amp_axis, p in phase_axis]
    phase, (; phase_axis, amp_axis, legend)
end

function draw_ft_table!(fig)
    n_row = length(val_IB)
    col_dens, col_cmpx, col_cohr, col_inco, col_tail = 1, 2, 3, 4, 5
    Label(fig[0, 1:5]; text="$tag FT2D: mean density, mean complex FT2D, coherent and incoherent spectra", tellwidth=false, halign=:left)
    for (col, title) in ((col_dens, "mean density"), (col_cmpx, "mean complex FT2D"), (col_cohr, "|mean FT2D| + coherent profile"), (col_inco, "mean |FT2D| + incoherent profile"), (col_tail, "tailess profiles"))
        Label(fig[1, col]; text=title, tellwidth=false, halign=:center, font=:bold)
    end
    draw_theme_colorbar!(gl, row_range; colorrange=clrrng_ft2d_absl_mean) = begin
        gl_cb = GridLayout(gl[row_range, 3])
        cb_x = [0.0, 1.0]
        cb_y = collect(range(colorrange...; length=256))
        cb = [y for _ in cb_x, y in cb_y]
        for j in 1:n_istp
            ax_cb = Axis(gl_cb[1, j]; ylabel=j == 1 ? "|FT|" : "")
            heatmap!(ax_cb, cb_x, cb_y, cb; colormap=clrmap_theme(j), colorrange=colorrange, rasterize=true)
            hidexdecorations!(ax_cb; grid=false)
            hideydecorations!(ax_cb; grid=false, label=j != 1, ticklabels=j != 1, ticks=j != 1)
        end
        colsize!(gl_cb, 1, Fixed(34)); colsize!(gl_cb, 2, Fixed(34)); colgap!(gl_cb, 2)
    end
    for (idx_IB, IB) in enumerate(val_IB)
        row = idx_IB + 1
        bottom = idx_IB == n_row
        Label(fig[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, halign=:right)
        gl_dens = GridLayout(fig[row, col_dens])
        gl_cmpx = GridLayout(fig[row, col_cmpx])
        gl_cohr = GridLayout(fig[row, col_cohr])
        gl_inco = GridLayout(fig[row, col_inco])
        gl_tail = GridLayout(fig[row, col_tail])
        for (idx_istp, istp_name) in enumerate(val_istp)
            clrmap = clrmap_theme(idx_istp)
            Label(gl_dens[2, idx_istp]; text=istp_name, fontsize=9, halign=:center)
            Label(gl_cmpx[2, idx_istp]; text=istp_name, fontsize=9, halign=:center)
            Label(gl_cohr[2, idx_istp]; text=istp_name, fontsize=9, halign=:center)
            Label(gl_inco[2, idx_istp]; text=istp_name, fontsize=9, halign=:center)
            Label(gl_tail[2, idx_istp]; text=istp_name, fontsize=9, halign=:center)
            ax_dens = Axis(gl_dens[3, idx_istp]; aspect=DataAspect(), xlabel=bottom ? "x (μm)" : "", ylabel=idx_istp == 1 ? "y (μm)" : "")
            heatmap!(ax_dens, x_dens, y_dens, ntfr2d_mean[idx_IB, idx_istp]'; colormap=clrmap, rasterize=true)
            hideydecorations!(ax_dens; label=idx_istp != 1, grid=false)
            ax_cmpx = Axis(gl_cmpx[3, idx_istp]; aspect=DataAspect(), xlabel=bottom ? "kx (μm⁻¹)" : "", ylabel=idx_istp == 1 ? "ky (μm⁻¹)" : "")
            heatmap!(ax_cmpx, kx_ft_plot, ky_ft_plot, first(gen_ft_rgba(ft2d_cmpx_mean[idx_IB, idx_istp]))[mask_ky_ft_plot, mask_kx_ft_plot]'; rasterize=true)
            xlims!(ax_cmpx, lims_kx_ft_plot); ylims!(ax_cmpx, lims_ky_ft_plot)
            hideydecorations!(ax_cmpx; label=idx_istp != 1, grid=false)
            ax_cohr = Axis(gl_cohr[3, idx_istp]; aspect=DataAspect(), xlabel="", ylabel=idx_istp == 1 ? "ky (μm⁻¹)" : "")
            amp_cohr_nonmasked, amp_cohr_masked = gen_amp_masked_clr(abs.(ft2d_cmpx_mean[idx_IB, idx_istp]), idx_istp; max=clrrng_ft2d_cmpx_mean[2])
            heatmap!(ax_cohr, kx_ft_plot, ky_ft_plot, amp_cohr_nonmasked[mask_ky_ft_plot, mask_kx_ft_plot]'; rasterize=true)
            heatmap!(ax_cohr, kx_ft_plot, ky_ft_plot, amp_cohr_masked[mask_ky_ft_plot, mask_kx_ft_plot]'; rasterize=true)
            hlines!(ax_cohr, [-ky_max_modl, ky_max_modl]; color=(:black, 0.55), linewidth=0.7, linestyle=:dash)
            xlims!(ax_cohr, lims_kx_ft_plot); ylims!(ax_cohr, lims_ky_ft_plot)
            ax_cohr_profile = Axis(gl_cohr[4, idx_istp]; xlabel=bottom ? "kx (μm⁻¹)" : "", ylabel=idx_istp == 1 ? "cohr" : "", xticks=0.0:0.1:0.5)
            lines!(ax_cohr_profile, kx_ft_plot, prfl_modl_fit_scaled_cohr[idx_IB, idx_istp][mask_kx_ft_plot]; color=(:gray40, 0.55), linewidth=1.0)
            lines!(ax_cohr_profile, kx_ft_plot, prfl_modl1d_cohr[idx_IB, idx_istp][mask_kx_ft_plot]; color=clr_theme(idx_istp), linewidth=1.2)
            ylims!(ax_cohr_profile, ylims_profile_cohr)
            linkxaxes!(ax_cohr, ax_cohr_profile)
            hideydecorations!(ax_cohr; label=idx_istp != 1, grid=false)
            hideydecorations!(ax_cohr_profile; label=idx_istp != 1, grid=false)
            ax_inco = Axis(gl_inco[3, idx_istp]; aspect=DataAspect(), xlabel="", ylabel=idx_istp == 1 ? "ky (μm⁻¹)" : "")
            amp_inco_nonmasked, amp_inco_masked = gen_amp_masked_clr(ft2d_absl_mean[idx_IB, idx_istp], idx_istp; max=clrrng_ft2d_absl_mean[2])
            heatmap!(ax_inco, kx_ft_plot, ky_ft_plot, amp_inco_nonmasked[mask_ky_ft_plot, mask_kx_ft_plot]'; rasterize=true)
            heatmap!(ax_inco, kx_ft_plot, ky_ft_plot, amp_inco_masked[mask_ky_ft_plot, mask_kx_ft_plot]'; rasterize=true)
            hlines!(ax_inco, [-ky_max_modl, ky_max_modl]; color=(:black, 0.55), linewidth=0.7, linestyle=:dash)
            xlims!(ax_inco, lims_kx_ft_plot); ylims!(ax_inco, lims_ky_ft_plot)
            ax_inco_profile = Axis(gl_inco[4, idx_istp]; xlabel=bottom ? "kx (μm⁻¹)" : "", ylabel=idx_istp == 1 ? "incohr" : "", xticks=0.0:0.1:0.5)
            lines!(ax_inco_profile, kx_ft_plot, prfl_modl1d_inco[idx_IB, idx_istp][mask_kx_ft_plot]; color=clr_theme(idx_istp), linewidth=1.3)
            ylims!(ax_inco_profile, ylims_profile_inco)
            linkxaxes!(ax_inco, ax_inco_profile)
            hideydecorations!(ax_inco; label=idx_istp != 1, grid=false)
            hideydecorations!(ax_inco_profile; label=idx_istp != 1, grid=false)
            ax_tail_profile = Axis(gl_tail[3, idx_istp]; xlabel=bottom ? "kx (μm⁻¹)" : "", ylabel=idx_istp == 1 ? "tailess" : "", xticks=0.0:0.1:0.5)
            lines!(ax_tail_profile, kx_ft_plot, prfl_modl1d_cohr_tailess[idx_IB, idx_istp][mask_kx_ft_plot]; color=clr_theme(idx_istp), linewidth=1.3)
            lines!(ax_tail_profile, kx_ft_plot, prfl_modl1d_inco_tailess[idx_IB, idx_istp][mask_kx_ft_plot]; color=clr_theme(idx_istp), linewidth=1.3, linestyle=:dash)
            ylims!(ax_tail_profile, min(ylims_profile_cohr[1], ylims_profile_inco[1]), max(ylims_profile_cohr[2], ylims_profile_inco[2]))
            hideydecorations!(ax_tail_profile; label=idx_istp != 1, grid=false)
        end
        for gl in (gl_dens, gl_cmpx)
            Label(gl[1, 1:2]; text=@sprintf("IB=%.3f A", IB), fontsize=9, halign=:left)
            rowsize!(gl, 1, Fixed(24)); rowsize!(gl, 2, Fixed(24)); rowsize!(gl, 3, Fixed(366)); rowgap!(gl, 0)
        end
        for gl in (gl_cohr, gl_inco)
            Label(gl[1, 1:3]; text=@sprintf("IB=%.3f A", IB), fontsize=9, halign=:left)
            rowsize!(gl, 1, Fixed(24)); rowsize!(gl, 2, Fixed(24)); rowsize!(gl, 3, Fixed(270)); rowsize!(gl, 4, Fixed(100))
            rowgap!(gl, 0)
        end
        Label(gl_tail[1, 1:2]; text=@sprintf("IB=%.3f A", IB), fontsize=9, halign=:left)
        rowsize!(gl_tail, 1, Fixed(24)); rowsize!(gl_tail, 2, Fixed(24)); rowsize!(gl_tail, 3, Fixed(100)); rowgap!(gl_tail, 0)
        colsize!(gl_dens, 1, Fixed(330)); colsize!(gl_dens, 2, Fixed(330))
        colsize!(gl_cmpx, 1, Fixed(285)); colsize!(gl_cmpx, 2, Fixed(285))
        colsize!(gl_cohr, 1, Fixed(285)); colsize!(gl_cohr, 2, Fixed(285))
        colsize!(gl_inco, 1, Fixed(285)); colsize!(gl_inco, 2, Fixed(285))
        colsize!(gl_tail, 1, Fixed(285)); colsize!(gl_tail, 2, Fixed(285))
        ax_cmpx_cb = Axis(gl_cmpx[3, 3]; xlabel="φ", ylabel="amp", titlesize=8)
        cmpx_cb = gen_ft_rgba(ft2d_cmpx_mean[idx_IB, 1])[2]
        heatmap!(ax_cmpx_cb, cmpx_cb.phase_axis, cmpx_cb.amp_axis, cmpx_cb.legend'; rasterize=true)
        xlims!(ax_cmpx_cb, -pi, pi); ylims!(ax_cmpx_cb, 0, maximum(cmpx_cb.amp_axis))
        draw_theme_colorbar!(gl_cohr, 3:4; colorrange=clrrng_ft2d_cmpx_mean)
        draw_theme_colorbar!(gl_inco, 3:4; colorrange=clrrng_ft2d_absl_mean)
        colsize!(gl_cmpx, 3, Fixed(80))
        colsize!(gl_cohr, 3, Fixed(80))
        colsize!(gl_inco, 3, Fixed(80))
        rowsize!(fig.layout, row, Fixed(420))
        colsize!(fig.layout, col_tail, Fixed(650))
    end
    rowsize!(fig.layout, 0, Fixed(28))
    rowsize!(fig.layout, 1, Fixed(28))
    colsize!(fig.layout, 0, Fixed(48))
    colsize!(fig.layout, col_dens, Fixed(700))
    colsize!(fig.layout, col_cmpx, Fixed(750))
    colsize!(fig.layout, col_cohr, Fixed(750))
    colsize!(fig.layout, col_inco, Fixed(750))
    colsize!(fig.layout, col_tail, Fixed(650))
    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)
    resize_to_layout!(fig)
    fig
end

fig_phase_distro = Figure(fontsize=12)
draw_ft_table!(fig_phase_distro)
for ext in ("png", "pdf")
    save(joinpath(path_output, "$(tag)_ft2d_table.$ext"), fig_phase_distro; backend=CairoMakie)
end

function draw_stacked_profile_heatmaps!(
    fig::Figure,
    row::Integer,
    prfl::AbstractArray{<:Real,3},
    x_plot::AbstractVector{<:Real},
    val_IB,
    val_istp,
    label_prfl;
    range_x_plot,
    colorrange,
)
    Label(fig[row - 2, 1:length(val_istp)]; text=label_prfl, tellwidth=false, halign=:center, font=:bold)
    hm = nothing
    for (idx_istp, istp_name) in enumerate(val_istp)
        ax = Axis(fig[row, idx_istp]; xlabel="kx (μm⁻¹)", ylabel="IB (A)", yaxisposition=idx_istp == 1 ? :left : :right, title="istp=$istp_name")
        hm = heatmap!(ax, x_plot, val_IB, @view(prfl[:, idx_istp, :]); colormap=gen_clrmap_solo(hue_theme_istp[istp_name]), colorrange, rasterize=true)
        idx_istp == 1 ? xlims!(ax, reverse(range_x_plot)) : xlims!(ax, range_x_plot)
        ax.xticks = idx_istp == 1 ? (last(range_x_plot):-0.1:first(range_x_plot)) : (first(range_x_plot):0.1:last(range_x_plot))
    end
    Colorbar(fig[row, length(val_istp) + 1], hm; label="profile")
    rowsize!(fig.layout, row - 2, Fixed(24))
    rowsize!(fig.layout, row - 1, Fixed(24))
    rowsize!(fig.layout, row, Fixed(260))
    colsize!(fig.layout, 1, Fixed(520))
    colsize!(fig.layout, 2, Fixed(520))
    colsize!(fig.layout, 3, Fixed(80))
    rowgap!(fig.layout, 0)
    resize_to_layout!(fig)
    fig
end

prfl_tailed_unmasked_fmt = Array{Float64}(undef, length(kx_ft), n_istp, n_IB)
prfl_tailess_masked_fmt = similar(prfl_tailed_unmasked_fmt)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    prfl_tailed_unmasked_fmt[:, idx_istp, idx_IB] .= prfl_modl1d_inco_unmasked[idx_IB, idx_istp]
    prfl_tailess_masked_fmt[:, idx_istp, idx_IB] .= prfl_modl1d_inco_tailess[idx_IB, idx_istp]
end
prfl_tailed_unmasked_fmt_cohr = similar(prfl_tailed_unmasked_fmt)
prfl_tailess_masked_fmt_cohr = similar(prfl_tailed_unmasked_fmt)
for idx_IB in 1:n_IB, idx_istp in 1:n_istp
    prfl_tailed_unmasked_fmt_cohr[:, idx_istp, idx_IB] .= prfl_modl1d_cohr_unmasked[idx_IB, idx_istp]
    prfl_tailess_masked_fmt_cohr[:, idx_istp, idx_IB] .= prfl_modl1d_cohr_tailess[idx_IB, idx_istp]
end
idx_x_plot = findall(mask_kx_ft_plot)
range_profile_stack = lims_kx_ft_plot
profile_tailess_min = minimum(prfl_tailess_masked_fmt[idx_x_plot, :, :])
profile_tailess_max = maximum(prfl_tailess_masked_fmt[idx_x_plot, :, :])
profile_range_inco = (min(0.0, profile_tailess_min), max(profile_tailess_max, eps(Float64)))
profile_tailess_cohr_min = minimum(prfl_tailess_masked_fmt_cohr[idx_x_plot, :, :])
profile_tailess_cohr_max = maximum(prfl_tailess_masked_fmt_cohr[idx_x_plot, :, :])
profile_range_cohr = (min(0.0, profile_tailess_cohr_min), max(profile_tailess_cohr_max, eps(Float64)))
fig_profiles_tailed = Figure(fontsize=12)
draw_stacked_profile_heatmaps!(fig_profiles_tailed, 3, prfl_tailed_unmasked_fmt[:, :, :], kx_ft, val_IB, val_istp, "tailed, unmasked incoherent profiles"; range_x_plot=range_profile_stack, colorrange=profile_range_inco)
draw_stacked_profile_heatmaps!(fig_profiles_tailed, 6, prfl_tailed_unmasked_fmt_cohr[:, :, :], kx_ft, val_IB, val_istp, "tailed, unmasked coherent profiles"; range_x_plot=range_profile_stack, colorrange=profile_range_cohr)
fig_profiles_tailess = Figure(fontsize=12)
draw_stacked_profile_heatmaps!(fig_profiles_tailess, 3, prfl_tailess_masked_fmt[:, :, :], kx_ft, val_IB, val_istp, "tailess, masked incoherent profiles"; range_x_plot=range_profile_stack, colorrange=profile_range_inco)
draw_stacked_profile_heatmaps!(fig_profiles_tailess, 6, prfl_tailess_masked_fmt_cohr[:, :, :], kx_ft, val_IB, val_istp, "tailess, masked coherent profiles"; range_x_plot=range_profile_stack, colorrange=profile_range_cohr)
for (fig_profiles, filename) in ((fig_profiles_tailed, "$(tag)_profiles_tailed_unmasked"), (fig_profiles_tailess, "$(tag)_profiles_masked_tailess"))
    for ext in ("png", "pdf")
        save(joinpath(path_output, "$filename.$ext"), fig_profiles; backend=CairoMakie)
    end
end

function draw_ft_live!(fig)
    obs_ib = Observable(ib)
    obs_istp = Observable(istp)
    obs_rep = Observable(idx_rep)
    dens_live() = dens_core[obs_ib[], obs_istp[]][obs_rep[]]
    ft_live() = ft2d_cmpx[obs_ib[], obs_istp[]][obs_rep[]]
    ft_mean_live() = ft2d_cmpx_mean[obs_ib[], obs_istp[]]
    inco_mean_live() = ft2d_absl_mean[obs_ib[], obs_istp[]]
    obs_density = Observable(dens_live()')
    obs_ft_shot = Observable(first(gen_ft_rgba(ft_live(); rng_amp=nothing))[mask_ky_ft_plot, mask_kx_ft_plot]')
    obs_ft_mean = Observable(first(gen_ft_rgba(ft_mean_live()))[mask_ky_ft_plot, mask_kx_ft_plot]')
    mean_amp_nonmasked, mean_amp_masked = gen_amp_masked_clr(abs.(ft_mean_live()), istp; max=clrrng_ft2d_cmpx_mean[2])
    inco_amp_nonmasked, inco_amp_masked = gen_amp_masked_clr(inco_mean_live(), istp)
    obs_ft_mean_amp_nonmasked = Observable(mean_amp_nonmasked[mask_ky_ft_plot, mask_kx_ft_plot]')
    obs_ft_mean_amp_masked = Observable(mean_amp_masked[mask_ky_ft_plot, mask_kx_ft_plot]')
    obs_ft_inco_nonmasked = Observable(inco_amp_nonmasked[mask_ky_ft_plot, mask_kx_ft_plot]')
    obs_ft_inco_masked = Observable(inco_amp_masked[mask_ky_ft_plot, mask_kx_ft_plot]')
    obs_cohr_profile = Observable(prfl_modl1d_cohr[ib, istp])
    obs_inco_profile = Observable(prfl_modl1d_inco[ib, istp])
    obs_cohr_reconstr = Observable(prfl_modl_fit_scaled_cohr[ib, istp])
    obs_inco_reconstr = Observable(prfl_modl_fit_scaled_inco[ib, istp])
    obs_density_range = Observable((0.0, maximum(dens_live())))
    obs_theme = Observable(clrmap_theme(istp))
    obs_title = lift(obs_ib, obs_istp, obs_rep) do i, j, r
        @sprintf("IB=%.3f A, istp=%s, rep=%d/%d", val_IB[i], val_istp[j], r, length(dens_core[i, j]))
    end
    Label(fig[0, 1:2]; text=obs_title, tellwidth=false, halign=:left)
    ax_density = Axis(fig[1, 1]; xlabel="x (μm)", ylabel="y (μm)", aspect=DataAspect())
    heatmap!(ax_density, x_dens, y_dens, obs_density; colormap=obs_theme, colorrange=obs_density_range, rasterize=true)
    ax_shot = Axis(fig[1, 2]; xlabel="kx (μm⁻¹)", ylabel="ky (μm⁻¹)", title="shot complex FT2D", aspect=DataAspect())
    heatmap!(ax_shot, kx_ft_plot, ky_ft_plot, obs_ft_shot; rasterize=true)
    xlims!(ax_shot, lims_kx_ft_plot); ylims!(ax_shot, lims_ky_ft_plot)

    gl_mean = GridLayout(fig[2, 1])
    ax_mean = Axis(gl_mean[1, 1]; xlabel="kx (μm⁻¹)", ylabel="ky (μm⁻¹)", title="mean complex FT2D", aspect=DataAspect())
    heatmap!(ax_mean, kx_ft_plot, ky_ft_plot, obs_ft_mean; rasterize=true)
    xlims!(ax_mean, lims_kx_ft_plot); ylims!(ax_mean, lims_ky_ft_plot)
    mean_rgba = gen_ft_rgba(ft_mean_live())[2]
    ax_mean_cb = Axis(gl_mean[1, 2]; xlabel="φ", ylabel="amp", title="FT shader")
    heatmap!(ax_mean_cb, mean_rgba.phase_axis, mean_rgba.amp_axis, mean_rgba.legend'; rasterize=true)
    xlims!(ax_mean_cb, -pi, pi); ylims!(ax_mean_cb, 0, maximum(mean_rgba.amp_axis))
    ax_ctrl = GridLayout(fig[2, 2])
    btn_ib_prev = Button(ax_ctrl[1, 1]; label="IB ←")
    btn_ib_next = Button(ax_ctrl[1, 2]; label="IB →")
    btn_istp_prev = Button(ax_ctrl[2, 1]; label="istp ←")
    btn_istp_next = Button(ax_ctrl[2, 2]; label="istp →")
    btn_rep_prev = Button(ax_ctrl[3, 1]; label="rep ←")
    btn_rep_next = Button(ax_ctrl[3, 2]; label="rep →")
    btn_reset = Button(ax_ctrl[4, 1:2]; label="reset")

    gl_cohr = GridLayout(fig[3, 1])
    ax_cohr = Axis(gl_cohr[1, 1]; xlabel="kx (μm⁻¹)", ylabel="ky (μm⁻¹)", title="|mean complex FT2D|", aspect=DataAspect())
    heatmap!(ax_cohr, kx_ft_plot, ky_ft_plot, obs_ft_mean_amp_nonmasked; rasterize=true)
    heatmap!(ax_cohr, kx_ft_plot, ky_ft_plot, obs_ft_mean_amp_masked; rasterize=true)
    hlines!(ax_cohr, [-ky_max_modl, ky_max_modl]; color=(:black, 0.55), linewidth=0.7, linestyle=:dash)
    xlims!(ax_cohr, lims_kx_ft_plot); ylims!(ax_cohr, lims_ky_ft_plot)
    ax_cohr_prfl = Axis(gl_cohr[2, 1]; xlabel="kx (μm⁻¹)", ylabel="cohr")
    lines!(ax_cohr_prfl, kx_ft_plot, lift(v -> v[mask_kx_ft_plot], obs_cohr_reconstr); color=(:gray40, 0.55), linewidth=1.0)
    lines!(ax_cohr_prfl, kx_ft_plot, lift(v -> v[mask_kx_ft_plot], obs_cohr_profile); color=lift(clr_theme, obs_istp), linewidth=1.3)
    ylims!(ax_cohr_prfl, ylims_profile_cohr); linkxaxes!(ax_cohr, ax_cohr_prfl)
    ax_cohr_cb = Axis(gl_cohr[1:2, 2]; ylabel="|FT|")
    heatmap!(ax_cohr_cb, [0.0, 1.0], collect(range(clrrng_ft2d_cmpx_mean...; length=256)), [y for _ in 1:2, y in range(clrrng_ft2d_cmpx_mean...; length=256)]; colormap=obs_theme, colorrange=clrrng_ft2d_cmpx_mean)
    hidexdecorations!(ax_cohr_cb; grid=false)

    gl_inco = GridLayout(fig[3, 2])
    ax_inco = Axis(gl_inco[1, 1]; xlabel="kx (μm⁻¹)", ylabel="ky (μm⁻¹)", title="mean |FT2D|", aspect=DataAspect())
    heatmap!(ax_inco, kx_ft_plot, ky_ft_plot, obs_ft_inco_nonmasked; rasterize=true)
    heatmap!(ax_inco, kx_ft_plot, ky_ft_plot, obs_ft_inco_masked; rasterize=true)
    hlines!(ax_inco, [-ky_max_modl, ky_max_modl]; color=(:black, 0.55), linewidth=0.7, linestyle=:dash)
    xlims!(ax_inco, lims_kx_ft_plot); ylims!(ax_inco, lims_ky_ft_plot)
    ax_inco_prfl = Axis(gl_inco[2, 1]; xlabel="kx (μm⁻¹)", ylabel="incohr")
    lines!(ax_inco_prfl, kx_ft_plot, lift(v -> v[mask_kx_ft_plot], obs_inco_reconstr); color=(:gray40, 0.55), linewidth=1.0)
    lines!(ax_inco_prfl, kx_ft_plot, lift(v -> v[mask_kx_ft_plot], obs_inco_profile); color=lift(clr_theme, obs_istp), linewidth=1.3)
    ylims!(ax_inco_prfl, ylims_profile_inco); linkxaxes!(ax_inco, ax_inco_prfl)
    ax_inco_cb = Axis(gl_inco[1:2, 2]; ylabel="|FT|")
    heatmap!(ax_inco_cb, [0.0, 1.0], collect(range(clrrng_ft2d_absl_mean...; length=256)), [y for _ in 1:2, y in range(clrrng_ft2d_absl_mean...; length=256)]; colormap=obs_theme, colorrange=clrrng_ft2d_absl_mean)
    hidexdecorations!(ax_inco_cb; grid=false)

    function update_live!()
        n = length(dens_core[obs_ib[], obs_istp[]])
        obs_rep[] = mod1(obs_rep[], n)
        obs_density[] = dens_live()'
        obs_ft_shot[] = first(gen_ft_rgba(ft_live(); rng_amp=nothing))[mask_ky_ft_plot, mask_kx_ft_plot]'
        obs_ft_mean[] = first(gen_ft_rgba(ft_mean_live()))[mask_ky_ft_plot, mask_kx_ft_plot]'
        mean_amp_nonmasked_live, mean_amp_masked_live = gen_amp_masked_clr(abs.(ft_mean_live()), obs_istp[]; max=clrrng_ft2d_cmpx_mean[2])
        inco_amp_nonmasked_live, inco_amp_masked_live = gen_amp_masked_clr(inco_mean_live(), obs_istp[])
        obs_ft_mean_amp_nonmasked[] = mean_amp_nonmasked_live[mask_ky_ft_plot, mask_kx_ft_plot]'
        obs_ft_mean_amp_masked[] = mean_amp_masked_live[mask_ky_ft_plot, mask_kx_ft_plot]'
        obs_ft_inco_nonmasked[] = inco_amp_nonmasked_live[mask_ky_ft_plot, mask_kx_ft_plot]'
        obs_ft_inco_masked[] = inco_amp_masked_live[mask_ky_ft_plot, mask_kx_ft_plot]'
        obs_cohr_profile[] = prfl_modl1d_cohr[obs_ib[], obs_istp[]]
        obs_inco_profile[] = prfl_modl1d_inco[obs_ib[], obs_istp[]]
        obs_cohr_reconstr[] = prfl_modl_fit_scaled_cohr[obs_ib[], obs_istp[]]
        obs_inco_reconstr[] = prfl_modl_fit_scaled_inco[obs_ib[], obs_istp[]]
        obs_density_range[] = (0.0, maximum(dens_live()))
        obs_theme[] = clrmap_theme(obs_istp[])
        notify(obs_title)
    end
    on(btn_ib_prev.clicks) do _; obs_ib[] = mod1(obs_ib[] - 1, n_IB); update_live!(); end
    on(btn_ib_next.clicks) do _; obs_ib[] = mod1(obs_ib[] + 1, n_IB); update_live!(); end
    on(btn_istp_prev.clicks) do _; obs_istp[] = mod1(obs_istp[] - 1, n_istp); update_live!(); end
    on(btn_istp_next.clicks) do _; obs_istp[] = mod1(obs_istp[] + 1, n_istp); update_live!(); end
    on(btn_rep_prev.clicks) do _; obs_rep[] -= 1; update_live!(); end
    on(btn_rep_next.clicks) do _; obs_rep[] += 1; update_live!(); end
    on(btn_reset.clicks) do _; obs_ib[] = ib; obs_istp[] = istp; obs_rep[] = idx_rep; update_live!(); end
    colsize!(fig.layout, 1, Fixed(400)); colsize!(fig.layout, 2, Fixed(400))
    rowsize!(fig.layout, 1, Fixed(300)); rowsize!(fig.layout, 2, Fixed(300)); rowsize!(fig.layout, 3, Fixed(400))
    colsize!(gl_cohr, 1, 240); colsize!(gl_cohr, 2, 60); rowsize!(gl_cohr, 1, 240); rowsize!(gl_cohr, 2, 120)
    colsize!(gl_inco, 1, 240); colsize!(gl_inco, 2, 60); rowsize!(gl_inco, 1, 240); rowsize!(gl_inco, 2, 120)
    colsize!(gl_mean, 1, 240); colsize!(gl_mean, 2, 60)
    resize_to_layout!(fig)
    fig
end

fig_live = Figure(fontsize=14)
draw_ft_live!(fig_live)
display(fig_live)


fig_cohr = Figure()
axs_cohr = Axis(fig_cohr[1, 1]; title="coherent 2D modulation spectral weight within mask")
axs_modl = Axis(fig_cohr[2, 1]; title="averaged 2D modulation spectral weight within mask")
filename_plot_cohr_modl = "$(tag)_cohr_modl"
for (i, istp) in enumerate(["162", "164"])
    lines!(axs_cohr, val_IB_ref, sum.(prfl_modl1d_cohr_tailess)[:,i]; color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.8))
    lines!(axs_modl, val_IB_ref, sum.(prfl_modl1d_inco_tailess)[:,i]; color=RGBAf(Oklch(0.52, 0.14, hue_theme_istp[istp]), 0.8))
end
for ext in ("png", "pdf")
    save(joinpath(path_output, "$filename_plot_cohr_modl.$ext"), fig_cohr; backend=CairoMakie)
end
println("  [$tag] saved phase distro table to $(joinpath(path_output, "$filename_plot_cohr_modl.png"))")

