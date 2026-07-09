using CairoMakie
using HDF5
using ImageFiltering
using Printf
using Statistics

isdefined(Main, :gen_clrmap_solo) || include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
isdefined(Main, :crop_center) || include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
isdefined(Main, :calc_prfl_modl_1d) || include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

path_root = isdefined(@__MODULE__, :path_root) ? path_root : raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\DualSS"
path_data_src = joinpath(path_root, "0204_interference", "result", "data.h5")
path_prfl_ref = joinpath(path_root, "0204_interference", "result", "prfl.h5")

tag = isdefined(@__MODULE__, :tag) ? tag : "SSNTFR"
val_istp = isdefined(@__MODULE__, :val_istp) ? val_istp : ["162", "164"]
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
smwh_src = (150, 150)
xy_fixed_src = (201, 201)
mag_src = 22.06
pixsz_src = 6.5
bin_src = 1
pad_src = 0
save_src_profiles = isdefined(@__MODULE__, :save_src_profiles) ? save_src_profiles : false
path_output_src = isdefined(@__MODULE__, :path_output) ? path_output : joinpath(path_root, "AnlzRoutine", "28.IncoCohrModlNtfr.[WL-migration]")
num_err_src = 0.6e4
range_center_src = 181:220
smwh_center_fit_src = (50, 50)
sigma_center_filter_src = 5

function calc_x_modl_src(w::Integer, smw_modl::Integer; mag::Real, pixsz::Real, bin::Integer=1)
    smw = (w - 1) ÷ 2
    x_modl_full = [mag / (smw * 2 * pixsz * bin) * n for n in -smw:smw]
    idx_center = (length(x_modl_full) + 1) ÷ 2
    return x_modl_full[idx_center-smw_modl:idx_center+smw_modl]
end

function load_src_density_payload(path_data::AbstractString, val_istp::AbstractVector{<:AbstractString})
    name_dataset_by_istp = Dict(
        "162" => "im64us",
        "164" => "im62us",
    )

    h5open(path_data, "r") do file
        unknown_istp = setdiff(val_istp, collect(keys(name_dataset_by_istp)))
        isempty(unknown_istp) || throw(ArgumentError("No source dataset mapping for istp values $unknown_istp."))

        dens_loaded = map(val_istp) do istp
            name_dataset = name_dataset_by_istp[istp]
            haskey(file, name_dataset) || throw(ArgumentError("Missing source dataset $name_dataset for istp=$istp."))
            read(file[name_dataset])
        end
        size_ref = size(first(dens_loaded))
        all(size(d) == size_ref for d in dens_loaded) || throw(DimensionMismatch(
            "All source density datasets must have the same size, got $(size.(dens_loaded)).",
        ))
        ndims(first(dens_loaded)) == 4 || throw(DimensionMismatch(
            "Expected raw source images with dimensions (y, x, rep, IB), got $size_ref.",
        ))

        _, _, n_rep, n_IB = size_ref
        dens_src = Array{Matrix{Float64}}(undef, n_IB, length(val_istp), n_rep)
        for idx_IB in 1:n_IB, idx_istp in eachindex(val_istp), idx_rep in 1:n_rep
            dens_src[idx_IB, idx_istp, idx_rep] = Float64.(copy(@view dens_loaded[idx_istp][:, :, idx_rep, idx_IB]))
        end
        return dens_src
    end
end

function check_reference_axes(
    path_prfl::AbstractString,
    x_dens::AbstractVector{<:Real},
    x_modl::AbstractVector{<:Real};
    atol::Real=1e-12,
)
    isfile(path_prfl) || return nothing
    x_dens_ref, x_modl_ref = h5open(path_prfl, "r") do file
        read(file["x_dens"]), read(file["x_modl"])
    end
    length(x_dens_ref) == length(x_dens) || throw(DimensionMismatch(
        "computed x_dens length $(length(x_dens)) must match prfl.h5 length $(length(x_dens_ref)).",
    ))
    length(x_modl_ref) == length(x_modl) || throw(DimensionMismatch(
        "computed x_modl length $(length(x_modl)) must match prfl.h5 length $(length(x_modl_ref)).",
    ))
    max_diff_x_dens = maximum(abs.(x_dens .- x_dens_ref))
    max_diff_x_modl = maximum(abs.(x_modl .- x_modl_ref))
    println("  [$tag] x_dens computed vs prfl.h5 max_abs_diff=$max_diff_x_dens")
    println("  [$tag] x_modl computed vs prfl.h5 max_abs_diff=$max_diff_x_modl")
    max_diff_x_dens <= atol || throw(DimensionMismatch("computed x_dens differs from prfl.h5 by max_abs_diff=$max_diff_x_dens."))
    if max_diff_x_modl <= atol
        println("  [$tag] using x_modl matching prfl.h5")
    else
        println("  [$tag] using recomputed symmetric x_modl; prfl.h5 x_modl is shifted by one bin")
    end
    return nothing
end

function gaussian_offset_1d(x, p)
    return @. p[1] * exp(-((x - p[2])^2) / (2 * p[3]^2)) + p[4]
end

function fit_gaussian_offset_center_1d(prfl::AbstractVector{<:Real})
    n = length(prfl)
    y = Float64.(prfl)
    x = collect(1.0:n)
    p0 = [maximum(y), (n + 1) / 2, n / 10, minimum(y)]
    fit = curve_fit(gaussian_offset_1d, x, y, p0)
    return fit.param[2]
end

function calc_center_peak_src(dens::AbstractMatrix{<:Real})
    dens_smooth = imfilter(dens, Kernel.gaussian(sigma_center_filter_src))
    x_center = fit_gaussian_offset_center_1d(vec(sum(dens_smooth; dims=1)))
    y_center = fit_gaussian_offset_center_1d(vec(sum(dens_smooth; dims=2)))
    return round.(Int, (x_center, y_center))
end

function calc_valid_duet_mask_src(
    dens_src_raw_fmt::AbstractArray{<:AbstractMatrix,3};
    num_err::Real=num_err_src,
    range_center=range_center_src,
)
    n_IB, n_istp, n_rep = size(dens_src_raw_fmt)
    num_src = Array{Float64}(undef, n_IB, n_istp, n_rep)
    xy_center_src = Array{Tuple{Int,Int}}(undef, n_IB, n_istp, n_rep)

    for idx_IB in 1:n_IB, idx_istp in 1:n_istp, idx_rep in 1:n_rep
        dens = dens_src_raw_fmt[idx_IB, idx_istp, idx_rep]
        num_src[idx_IB, idx_istp, idx_rep] = sum(dens)
        xy_center_src[idx_IB, idx_istp, idx_rep] = calc_center_peak_src(dens)
    end

    num_median_src = dropdims(median(num_src; dims=3); dims=3)
    mask_valid_duet = falses(n_IB, n_rep)
    for idx_IB in 1:n_IB, idx_rep in 1:n_rep
        is_valid_number = all(
            abs(num_src[idx_IB, idx_istp, idx_rep] - num_median_src[idx_IB, idx_istp]) <= num_err
            for idx_istp in 1:n_istp
        )
        is_valid_center = all(
            (xy -> xy[1] in range_center && xy[2] in range_center)(xy_center_src[idx_IB, idx_istp, idx_rep])
            for idx_istp in 1:n_istp
        )
        mask_valid_duet[idx_IB, idx_rep] = is_valid_number && is_valid_center
    end

    return (; mask_valid_duet, xy_center_src, num_src, num_median_src)
end

function crop_valid_source_densities(
    dens_src_raw_fmt::AbstractArray{<:AbstractMatrix,3},
    mask_valid_duet::AbstractMatrix{Bool},
    xy_center_src::AbstractArray{<:Tuple{Int,Int},3},
    smwh::Tuple{<:Integer,<:Integer},
)
    n_IB, n_istp, n_rep = size(dens_src_raw_fmt)
    size(mask_valid_duet) == (n_IB, n_rep) || throw(DimensionMismatch(
        "mask_valid_duet size $(size(mask_valid_duet)) must match (IB, rep) $((n_IB, n_rep)).",
    ))
    size(xy_center_src) == size(dens_src_raw_fmt) || throw(DimensionMismatch(
        "xy_center_src size $(size(xy_center_src)) must match density source size $(size(dens_src_raw_fmt)).",
    ))

    dens_src_core = Array{Vector{Matrix{Float64}}}(undef, n_IB, n_istp)
    for idx_IB in 1:n_IB, idx_istp in 1:n_istp
        dens_src_core[idx_IB, idx_istp] = [
            crop_center(dens_src_raw_fmt[idx_IB, idx_istp, idx_rep], xy_center_src[idx_IB, idx_istp, idx_rep], smwh) |> copy
            for idx_rep in 1:n_rep
            if mask_valid_duet[idx_IB, idx_rep]
        ]
    end
    return dens_src_core
end

function save_src_profiles_h5(
    path_output::AbstractString,
    x_dens::AbstractVector{<:Real},
    x_modl::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector,
    ntfr2d_mean::AbstractMatrix{<:AbstractMatrix},
    prfl_inco::AbstractArray{<:Real,3},
    prfl_cohr::AbstractArray{<:Real,3},
    count_profile_shot::AbstractVector{<:Integer},
)
    isdir(path_output) || mkpath(path_output)
    path_prfl_src = joinpath(path_output, "$(tag)_prfl_src.h5")
    ntfr2d_mean_packed = Array{Float64}(undef, length(x_dens), length(x_dens), length(val_istp), length(val_IB))
    for idx_IB in eachindex(val_IB), idx_istp in eachindex(val_istp)
        ntfr2d_mean_packed[:, :, idx_istp, idx_IB] .= ntfr2d_mean[idx_IB, idx_istp]
    end
    h5open(path_prfl_src, "w") do file
        write(file, "x_dens", x_dens)
        write(file, "x_modl", x_modl)
        write(file, "val_IB", val_IB)
        write(file, "val_istp", String.(val_istp))
        write(file, "ntfr2d_mean", ntfr2d_mean_packed)
        write(file, "prfl_inco", prfl_inco)
        write(file, "prfl_cohr", prfl_cohr)
        write(file, "count_profile_shot", collect(count_profile_shot))
    end
    println("saved $path_prfl_src")
    return path_prfl_src
end

function draw_ntfr2d_mean_src_table!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    ntfr2d_mean::AbstractMatrix{<:AbstractMatrix};
    colorrange,
    smh_dens_strip::Integer,
)
    for (idx_istp, istp) in enumerate(val_istp)
        Label(fig[1, idx_istp]; text="istp=$istp", tellwidth=false, tellheight=true, halign=:center, font=:bold)
    end

    for (idx_IB, IB) in enumerate(val_IB)
        row = idx_IB + 1
        Label(fig[row, 0]; text=@sprintf("%.3f", IB), tellwidth=true, tellheight=false, halign=:right)
        for (idx_istp, istp) in enumerate(val_istp)
            ax = Axis(
                fig[row, idx_istp];
                xlabel=idx_IB == length(val_IB) ? "position (μm)" : "",
                ylabel=idx_istp == 1 ? "position (μm)" : "",
                aspect=DataAspect(),
            )
            heatmap!(
                ax,
                x_dens,
                x_dens,
                ntfr2d_mean[idx_IB, idx_istp]';
                colormap=gen_clrmap_solo(hue_theme_istp[istp]),
                colorrange,
                rasterize=true,
            )
            step_x = median(diff(x_dens))
            y_strip = (smh_dens_strip + 0.5) * step_x
            x_rect = [
                first(x_dens) - step_x / 2,
                last(x_dens) + step_x / 2,
                last(x_dens) + step_x / 2,
                first(x_dens) - step_x / 2,
                first(x_dens) - step_x / 2,
            ]
            y_rect = [
                -y_strip,
                -y_strip,
                y_strip,
                y_strip,
                -y_strip,
            ]
            lines!(ax, x_rect, y_rect; color=(Oklch(0.4, 0.05, 320), 0.40), linewidth=0.8)
            hidexdecorations!(ax; label=idx_IB == length(val_IB) ? false : true, ticklabels=idx_IB == length(val_IB) ? false : true, ticks=idx_IB == length(val_IB) ? false : true, grid=false)
            hideydecorations!(ax; label=idx_istp == 1 ? false : true, ticklabels=false, ticks=false, grid=false)
        end
        rowsize!(fig.layout, row, Fixed(105))
    end
    return fig
end

function save_ntfr2d_mean_src_table(
    path_output::AbstractString,
    x_dens::AbstractVector{<:Real},
    val_IB::AbstractVector{<:Real},
    val_istp::AbstractVector{<:AbstractString},
    ntfr2d_mean::AbstractMatrix{<:AbstractMatrix},
    smh_dens_strip::Integer,
)

end

println("  [$tag] loading source densities from $path_data_src")
dens_src_raw_fmt = load_src_density_payload(path_data_src, val_istp)
n_IB_src, n_istp_src, n_rep_src = size(dens_src_raw_fmt)
wh_raw_src = size(dens_src_raw_fmt[1, 1, 1])
println("  [$tag] formatted source densities as (IB, istp, rep)=$(size(dens_src_raw_fmt)), image size=$wh_raw_src")
length(val_IB_ref) == n_IB_src || throw(DimensionMismatch("val_IB_ref length $(length(val_IB_ref)) must match source IB count $n_IB_src."))
length(val_istp) == n_istp_src || throw(DimensionMismatch("val_istp length $(length(val_istp)) must match source istp count $n_istp_src."))

cfg_prfl_modl = get_prfl_modl_1d_config(smwh_src)
x_dens = (pixsz_src * bin_src / mag_src) .* collect(-smwh_src[2]:smwh_src[2])
x_modl = calc_x_modl_src(2 * smwh_src[1] + 1, cfg_prfl_modl.smw_modl; mag=mag_src, pixsz=pixsz_src, bin=bin_src)
check_reference_axes(path_prfl_ref, x_dens, x_modl)
step_modl = median(diff(x_modl))
val_IB = copy(val_IB_ref)

valid_src = calc_valid_duet_mask_src(dens_src_raw_fmt)
count_profile_shot = vec(sum(valid_src.mask_valid_duet; dims=2))
println("  [$tag] valid source duet counts per IB=$(count_profile_shot)")
dens_src_core = crop_valid_source_densities(dens_src_raw_fmt, valid_src.mask_valid_duet, valid_src.xy_center_src, smwh_src)
ntfr2d_mean = map(dens_src_core) do ds
    isempty(ds) && throw(ArgumentError("No valid source densities available for a condition."))
    dropdims(mean(stack(ds); dims=3); dims=3)
end

prfl_modl_src = map(ds -> calc_prfl_modl_1d(ds, smwh_src; step_modl), dens_src_core)
    
prfl_inco = Array{Float64}(undef, length(x_modl), length(val_istp), n_IB_src)
prfl_cohr = similar(prfl_inco)
for idx_IB in axes(prfl_modl_src, 1), idx_istp in axes(prfl_modl_src, 2)
    prfl = prfl_modl_src[idx_IB, idx_istp]
    length(prfl.prfl_inco) == length(x_modl) || throw(DimensionMismatch("prfl_inco length $(length(prfl.prfl_inco)) must match x_modl length $(length(x_modl))."))
    prfl_inco[:, idx_istp, idx_IB] .= prfl.prfl_inco
    prfl_cohr[:, idx_istp, idx_IB] .= prfl.prfl_cohr
end

isdir(path_output_src) || mkpath(path_output_src)
colorrange = (0.0, maximum(maximum, ntfr2d_mean))
fig = Figure(fontsize=14)
Label(
    fig[0, 1:length(val_istp)];
    text=@sprintf("%s source mean densities from number/center-selected peak crops, common colorrange %.3g..%.3g", tag, colorrange...),
    tellwidth=false,
    tellheight=true,
    halign=:left,
)
draw_ntfr2d_mean_src_table!(fig, x_dens, val_IB, val_istp, ntfr2d_mean; colorrange, cfg_prfl_modl.smh_dens_strip)
colsize!(fig.layout, 1, Fixed(105))
colsize!(fig.layout, 2, Fixed(105))
colgap!(fig.layout, 1, 8)
rowgap!(fig.layout, 1, 4)
resize_to_layout!(fig)
path_plot = joinpath(path_output_src, "$(tag)_ntfr2d_mean_src_table.png")
save(path_plot, fig; backend=CairoMakie)
println("saved $path_plot")

if save_src_profiles
    save_src_profiles_h5(path_output_src, x_dens, x_modl, val_IB, val_istp, ntfr2d_mean, prfl_inco, prfl_cohr, count_profile_shot)
end
