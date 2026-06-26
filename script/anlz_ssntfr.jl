using CairoMakie
using DelimitedFiles
using FFTW
using HDF5
using ImageFiltering
using LsqFit
using Printf
using Statistics

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
title_anlz = "01.Reproduce"
path_data = joinpath(path_root, "0204_interference", "result", "data.h5")
path_sidepeak_ref = joinpath(path_root, "0204_interference", "result", "sidepeak_B.csv")
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

tag = "SSNTFR"
IBfesh = [
    5.310, 5.312, 5.314, 5.316, 5.317, 5.318, 5.319, 5.320, 5.322,
    5.324, 5.326, 5.328, 5.330, 5.332, 5.334, 5.338, 5.342,
]
# The original notebook labels the raw HDF5 isotope datasets as 62/64, but its
# plotting convention treats those labels as inverted.
istp = ["162", "164"]

avg = 100
nerr = 0.6e4
smwh_roi = (150, 150)
center_window = 181:220
mag = 22.06
pixsz = 6.5
pad = 0
fsz = 20
fszy = 120
fit_range_n = (0.10, 0.40)
fit_range_phi = (0.10, 0.40)

idx_IB_axis = 1
idx_rep_axis = 2
idx_t_hold_axis = 3
idx_istp_axis = 4
idx_162 = findfirst(==("162"), istp)
idx_164 = findfirst(==("164"), istp)

function gaussian_center_1d(prfl::AbstractVector{<:Real})
    x = collect(1.0:length(prfl))
    y = Float64.(prfl)
    amp = maximum(y) - minimum(y)
    p0 = [amp, Float64(argmax(y)), max(length(y) / 10, 2.0), minimum(y)]
    lower = [0.0, 1.0, 1.0, -Inf]
    upper = [Inf, length(y), length(y), Inf]
    model_center(x, p) = @. p[1] * exp(-((x - p[2])^2) / (2p[3]^2)) + p[4]
    fit = curve_fit(model_center, x, y, p0; lower, upper)
    return round(Int, coef(fit)[2])
end

function find_center(img::AbstractMatrix{<:Real}; sigma::Real=10)
    smoothed = imfilter(Float64.(img), Kernel.gaussian(sigma))
    cy = gaussian_center_1d(vec(sum(smoothed; dims=2)))
    cx = gaussian_center_1d(vec(sum(smoothed; dims=1)))
    return cy, cx
end

function crop_center_yx(img::AbstractMatrix{<:Real}, cy::Integer, cx::Integer, smwh::Tuple{<:Integer,<:Integer})
    smx, smy = smwh
    return @view img[cy-smy:cy+smy, cx-smx:cx+smx]
end

function format_dens_ssntfr_h5_ib(f::HDF5.File, idx_ib::Integer, val_IB::Real, vals_rep, vals_istp)
    dens62 = f["im62us"][:, :, :, idx_ib]
    dens64 = f["im64us"][:, :, :, idx_ib]
    size(dens62, 3) == length(vals_rep) || throw(DimensionMismatch("rep count mismatch for IB=$val_IB."))

    val_vars = (; IB=[val_IB], rep=collect(vals_rep), t_hold=[0.0], istp=collect(vals_istp))
    name_dims = propertynames(val_vars)
    n_dim_vars = Tuple(length(getproperty(val_vars, name)) for name in name_dims)
    dens_full_fmt = Array{Matrix{Float64}}(undef, n_dim_vars...)
    for idx_rep in eachindex(vals_rep)
        dens_full_fmt[1, idx_rep, 1, idx_164] = Float64.(dens62[:, :, idx_rep])
        dens_full_fmt[1, idx_rep, 1, idx_162] = Float64.(dens64[:, :, idx_rep])
    end
    return (; val_vars, dens_full_fmt, name_dims, n_dim_vars, wh_dens=reverse(size(dens62)[1:2]))
end

function calc_postselection(dens_full_fmt, val_IB::Real)
    dens62 = @view dens_full_fmt[1, :, 1, idx_164]
    dens64 = @view dens_full_fmt[1, :, 1, idx_162]
    n62_all = sum.(dens62)
    n64_all = sum.(dens64)
    n62_med = median(n62_all)
    n64_med = median(n64_all)

    mask_rep = falses(length(dens62))
    dens_roi_fmt = Array{Union{Missing,Matrix{Float64}}}(missing, 1, length(dens62), 1, length(istp))
    num_fmt = Array{Float64}(undef, 1, length(dens62), 1, length(istp))
    center_fmt = Array{Tuple{Int,Int}}(undef, 1, length(dens62), 1, length(istp))

    for idx_rep in eachindex(dens62)
        img62 = dens62[idx_rep]
        img64 = dens64[idx_rep]
        n62 = n62_all[idx_rep]
        n64 = n64_all[idx_rep]
        cy62, cx62 = find_center(img62)
        cy64, cx64 = find_center(img64)
        keep = abs(n62 - n62_med) <= nerr &&
               abs(n64 - n64_med) <= nerr &&
               all(in(center_window), (cy62, cx62, cy64, cx64))
        num_fmt[1, idx_rep, 1, idx_164] = n62
        num_fmt[1, idx_rep, 1, idx_162] = n64
        center_fmt[1, idx_rep, 1, idx_164] = (cx62, cy62)
        center_fmt[1, idx_rep, 1, idx_162] = (cx64, cy64)
        if keep
            mask_rep[idx_rep] = true
            dens_roi_fmt[1, idx_rep, 1, idx_164] = copy(crop_center_yx(img62, cy62, cx62, smwh_roi))
            dens_roi_fmt[1, idx_rep, 1, idx_162] = copy(crop_center_yx(img64, cy64, cx64, smwh_roi))
        end
    end

    any(mask_rep) || throw(ArgumentError("post-selection kept no shots for IB=$val_IB"))
    return (; mask_rep, dens_roi_fmt, num_fmt, center_fmt)
end

function mean_selected_images(dens_roi_fmt, mask_rep::AbstractVector{Bool}, idx_istp::Integer)
    images = skipmissing(vec(dens_roi_fmt[1, mask_rep, 1, idx_istp])) |> collect
    isempty(images) && throw(ArgumentError("cannot average an empty selected image list."))
    return mean(images)
end

function tukey_window(n::Integer, alpha::Real=0.2)
    m = n - 1
    edge = floor(Int, alpha * m / 2)
    win = ones(Float64, n)
    edge == 0 && return win
    for idx in 1:n
        k = idx - 1
        if k < edge
            win[idx] = 0.5 * (1 - cos(2pi * k / (alpha * m)))
        elseif k > m - edge
            win[idx] = 0.5 * (1 - cos(2pi * (m - k) / (alpha * m)))
        end
    end
    return win
end

fftshift1(v::AbstractVector) = circshift(v, -((length(v) + 1) ÷ 2))

function calc_fft_profile(images::AbstractVector{<:AbstractMatrix{<:Real}})
    ly, lx = size(first(images))
    szy = (ly + 1) ÷ 2
    szx = (lx + 1) ÷ 2
    fftsz = (ly - 1) / 2 + pad
    fftx = [mag / (fftsz * 2 * pixsz) * n for n in -fftsz:fftsz]
    rxf = (szx - fsz):(szx + fsz)
    ryf = (szy + pad - fszy):(szy + pad + fszy)
    win = tukey_window(ly, 0.2)
    dk = fftx[2] - fftx[1]

    fft_reps = map(images) do img
        prfl = vec(sum(imfilter(Float64.(@view img[:, rxf]), Kernel.gaussian(5)); dims=2)) ./ length(rxf)
        ft = fftshift1(fft(prfl .* win))[ryf]
        ft ./ (sum(abs.(ft)) * dk / 2)
    end

    return (;
        fftx=fftx[ryf],
        fn=mean(abs.(stack(fft_reps)); dims=2) |> vec,
        fphi=abs.(mean(stack(fft_reps); dims=2) |> vec),
        fft_reps,
    )
end

function calc_fft_profile_fmt(dens_roi_fmt, mask_rep::AbstractVector{Bool})
    fft_fmt = Array{NamedTuple}(undef, 1, 1, length(istp))
    for idx_istp in eachindex(istp)
        images = skipmissing(vec(dens_roi_fmt[1, mask_rep, 1, idx_istp])) |> collect
        fft_fmt[1, 1, idx_istp] = calc_fft_profile(images)
    end
    return fft_fmt
end

model_n(x, p) = @. p[5] + p[4] * x + p[1]^2 * exp(-((x - p[2])^2) / (2p[3]^2))
model_phi(x, p) = @. p[8] + p[7] * x + p[1]^2 * exp(-((x - p[2])^2) / (2p[3]^2)) +
                     p[4]^2 * exp(-((x - p[6])^2) / (2p[5]^2))

function relative_error_a0(fit)
    cov = try
        estimate_covar(fit)
    catch
        fill(NaN, length(coef(fit)), length(coef(fit)))
    end
    se_a0 = sqrt(abs(cov[1, 1]))
    a0 = abs(coef(fit)[1])
    return iszero(a0) ? NaN : se_a0 / a0
end

function fit_sidepeak(fftx::AbstractVector, prfl::AbstractVector, kind::Symbol)
    lo, hi = kind == :phi ? fit_range_phi : fit_range_n
    mask = (lo .<= fftx) .& (fftx .<= hi)
    x = Float64.(fftx[mask])
    y = Float64.(prfl[mask])
    amp_hint = sqrt(max(maximum(y) - median(y), 1e-6))
    fits = Any[]
    if kind == :n
        lower = [1e-8, 0.10, 0.03, -Inf, -Inf]
        upper = [Inf, 0.27, 0.07, 0.0, Inf]
        for x0 in (0.16, 0.18, 0.21, 0.24), sigma0 in (0.04, 0.05, 0.06)
            p0 = [amp_hint, x0, sigma0, -1.0, minimum(y)]
            try
                push!(fits, curve_fit(model_n, x, y, p0; lower, upper, maxIter=10_000))
            catch
            end
        end
    else
        lower = [1e-8, 0.16, 0.01, 1e-8, 0.01, -Inf, -Inf, -Inf]
        upper = [Inf, 0.27, 0.07, Inf, 0.20, 0.10, 0.0, Inf]
        for x0 in (0.18, 0.21, 0.24), sigma0 in (0.02, 0.03, 0.05), x1 in (0.0, 0.05, 0.09)
            p0 = [amp_hint, x0, sigma0, amp_hint / 2, 0.10, x1, -1.0, minimum(y)]
            try
                push!(fits, curve_fit(model_phi, x, y, p0; lower, upper, maxIter=10_000))
            catch
            end
        end
    end
    isempty(fits) && throw(ArgumentError("all $kind sidepeak fit attempts failed."))
    fit = fits[argmin([sum(abs2, residuals(f)) for f in fits])]
    return (; params=coef(fit), amp=abs(coef(fit)[1]^2), relerr=relative_error_a0(fit), x, y)
end

function fit_np_from_fft_fmt(fft_fmt)
    out = Array{NamedTuple}(undef, 1, length(istp), 2)
    for idx_istp in eachindex(istp)
        fft = fft_fmt[1, 1, idx_istp]
        out[1, idx_istp, 1] = fit_sidepeak(fft.fftx, fft.fn, :n)
        out[1, idx_istp, 2] = fit_sidepeak(fft.fftx, fft.fphi, :phi)
    end
    return out
end

function load_sidepeak_reference_fmt(path_sidepeak::AbstractString)
    data = readdlm(path_sidepeak, ',', Float64)
    size(data, 2) >= 10 || throw(DimensionMismatch("Expected at least 10 columns in $path_sidepeak, got $(size(data, 2))."))
    amp_fmt = Array{Float64}(undef, size(data, 1), length(istp), 2)
    err_fmt = similar(amp_fmt)
    amp_fmt[:, idx_162, 1] .= data[:, 3]
    err_fmt[:, idx_162, 1] .= data[:, 3] .* data[:, 4]
    amp_fmt[:, idx_164, 1] .= data[:, 5]
    err_fmt[:, idx_164, 1] .= data[:, 5] .* data[:, 6]
    amp_fmt[:, idx_162, 2] .= data[:, 7]
    err_fmt[:, idx_162, 2] .= data[:, 7] .* data[:, 8]
    amp_fmt[:, idx_164, 2] .= data[:, 9]
    err_fmt[:, idx_164, 2] .= data[:, 9] .* data[:, 10]
    return (; IB=data[:, 1], amp_fmt, err_fmt)
end

function draw_reference_style_plot(path_save::AbstractString, x, amp_fmt, err_fmt, idx_prop::Integer; ylabel::AbstractString, ylim)
    fig = Figure(size=(810, 510), backgroundcolor=:white, fontsize=24)
    ax = Axis(
        fig[1, 1];
        xlabel="xlabel",
        ylabel,
        xlabelsize=32,
        ylabelsize=32,
        xticklabelsize=24,
        yticklabelsize=24,
        spinewidth=1.2,
        xgridvisible=false,
        ygridvisible=false,
        xticks=5.310:0.005:5.340,
        xminorticksvisible=true,
        yminorticksvisible=true,
        xticksvisible=true,
        yticksvisible=true,
        xminorgridvisible=false,
        yminorgridvisible=false,
    )
    xlims!(ax, minimum(x) - 0.0015, maximum(x) + 0.0015)
    ylims!(ax, ylim)
    ax.xtickalign = 1
    ax.ytickalign = 1
    ax.xminortickalign = 1
    ax.yminortickalign = 1

    y164 = amp_fmt[:, idx_164, idx_prop]
    e164 = err_fmt[:, idx_164, idx_prop]
    y162 = amp_fmt[:, idx_162, idx_prop]
    e162 = err_fmt[:, idx_162, idx_prop]
    errorbars!(ax, x, y164, e164; color=:blue, whiskerwidth=0, linewidth=2)
    ln164 = scatterlines!(ax, x, y164; color=:blue, marker=:circle, markersize=10, linewidth=3, linestyle=:dash)
    errorbars!(ax, x, y162, e162; color=:red, whiskerwidth=0, linewidth=2)
    ln162 = scatterlines!(ax, x, y162; color=:red, marker=:rect, markersize=10, linewidth=3, linestyle=:dash)
    axislegend(ax, [ln164, ln162], ["164Dy", "162Dy"]; position=(0.73, 0.82), framevisible=false, labelsize=24, patchsize=(70, 18))
    save(path_save, fig)
    return fig
end

println("Loading $path_data")
val_vars = (; IB=IBfesh, rep=collect(1:avg), t_hold=[0.0], istp)
name_dims = propertynames(val_vars)
n_dim_vars = Tuple(length(getproperty(val_vars, name)) for name in name_dims)
post_mask_fmt = falses(length(IBfesh), avg)
num_stat_fmt = Array{Tuple{Float64,Float64}}(undef, length(IBfesh), length(istp))
dens_avg_fmt = Array{Matrix{Float64}}(undef, length(IBfesh), 1, length(istp))
fft_prfl_fmt = Array{NamedTuple}(undef, length(IBfesh), 1, length(istp))
fit_np_fmt = Array{NamedTuple}(undef, length(IBfesh), length(istp), 2)

h5open(path_data, "r") do f
    size(f["im62us"], 3) == avg || @warn "Expected $avg reps, got $(size(f["im62us"], 3))."
    size(f["im62us"], 4) == length(IBfesh) || throw(DimensionMismatch("IB axis mismatch."))
    println("Formatted axes $name_dims dims $n_dim_vars")
    for idx_ib in eachindex(IBfesh)
        println("  [$tag] formatting/post-selecting IB=$(IBfesh[idx_ib]) ($(idx_ib)/$(length(IBfesh)))")
        fmt = format_dens_ssntfr_h5_ib(f, idx_ib, IBfesh[idx_ib], val_vars.rep, val_vars.istp)
        post = calc_postselection(fmt.dens_full_fmt, IBfesh[idx_ib])
        post_mask_fmt[idx_ib, :] .= post.mask_rep
        for idx_istp in eachindex(istp)
            nums = vec(post.num_fmt[1, post.mask_rep, 1, idx_istp])
            num_stat_fmt[idx_ib, idx_istp] = (mean(nums), std(nums))
            dens_avg_fmt[idx_ib, 1, idx_istp] = mean_selected_images(post.dens_roi_fmt, post.mask_rep, idx_istp)
        end
        fft_ib_fmt = calc_fft_profile_fmt(post.dens_roi_fmt, post.mask_rep)
        fit_ib_fmt = fit_np_from_fft_fmt(fft_ib_fmt)
        fft_prfl_fmt[idx_ib:idx_ib, :, :] .= fft_ib_fmt
        fit_np_fmt[idx_ib:idx_ib, :, :] .= fit_ib_fmt
    end
end

n_kept = vec(sum(post_mask_fmt; dims=2))
println("Post-selected shots per IB: ", n_kept)

np_fit_amp_fmt = Array{Float64}(undef, length(IBfesh), length(istp), 2)
np_fit_err_fmt = similar(np_fit_amp_fmt)
for idx_ib in eachindex(IBfesh), idx_istp in eachindex(istp), idx_prop in 1:2
    np_fit_amp_fmt[idx_ib, idx_istp, idx_prop] = fit_np_fmt[idx_ib, idx_istp, idx_prop].amp
    np_fit_err_fmt[idx_ib, idx_istp, idx_prop] = fit_np_fmt[idx_ib, idx_istp, idx_prop].amp * fit_np_fmt[idx_ib, idx_istp, idx_prop].relerr
end

draw_reference_style_plot(
    joinpath(path_output, "ssntfr_np_fit_An.png"),
    IBfesh,
    np_fit_amp_fmt,
    np_fit_err_fmt,
    1;
    ylabel="An(a.u.)",
    ylim=(0, max(2.25, maximum(np_fit_amp_fmt[:, :, 1]) * 1.15)),
)
draw_reference_style_plot(
    joinpath(path_output, "ssntfr_np_fit_Aphi.png"),
    IBfesh,
    np_fit_amp_fmt,
    np_fit_err_fmt,
    2;
    ylabel="Aϕ(a.u.)",
    ylim=(0, max(0.7, maximum(np_fit_amp_fmt[:, :, 2]) * 1.15)),
)

if isfile(path_sidepeak_ref)
    ref = load_sidepeak_reference_fmt(path_sidepeak_ref)
    draw_reference_style_plot(
        joinpath(path_output, "Ref,An.reproduce.png"),
        ref.IB,
        ref.amp_fmt,
        ref.err_fmt,
        1;
        ylabel="An(a.u.)",
        ylim=(0, 2.25),
    )
    draw_reference_style_plot(
        joinpath(path_output, "Ref.Aphi.reproduce.png"),
        ref.IB,
        ref.amp_fmt,
        ref.err_fmt,
        2;
        ylabel="Aϕ(a.u.)",
        ylim=(0, 0.7),
    )
end

open(joinpath(path_output, "ssntfr_np_fit.csv"), "w") do io
    println(io, "IB,np64n,err64n,np62n,err62n,np64phi,err64phi,np62phi,err62phi,n_kept")
    for idx_ib in eachindex(IBfesh)
        @printf(io, "%.6f,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%d\n",
            IBfesh[idx_ib],
            np_fit_amp_fmt[idx_ib, idx_164, 1], np_fit_err_fmt[idx_ib, idx_164, 1],
            np_fit_amp_fmt[idx_ib, idx_162, 1], np_fit_err_fmt[idx_ib, idx_162, 1],
            np_fit_amp_fmt[idx_ib, idx_164, 2], np_fit_err_fmt[idx_ib, idx_164, 2],
            np_fit_amp_fmt[idx_ib, idx_162, 2], np_fit_err_fmt[idx_ib, idx_162, 2],
            n_kept[idx_ib])
    end
end

cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
println("Saved SSNTFR reproduction outputs to $path_output")
