using YAML
using MAT
using CairoMakie
using Printf

include(joinpath(@__DIR__, "dualmotcommons.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
CairoMakie.activate!()
function mat_field_value(container, name::AbstractString, label::AbstractString)
    while container isa AbstractArray && !(container isa MAT.MatlabStructArray)
        length(container) == 1 ||
            throw(ArgumentError("$label: expected one MATLAB struct value, got $(length(container))"))
        container = only(container)
    end
    if container isa MAT.MatlabStructArray
        return container[name]
    elseif container isa AbstractDict && haskey(container, name)
        return container[name]
    end
    throw(ArgumentError("$label: missing MATLAB struct field '$name'"))
end

function read_insitu_xscale(path_liferes::AbstractString)
    liferes = matread(path_liferes)["liferes"]
    cres = mat_field_value(liferes, "cres", path_liferes)
    entry = if cres isa MAT.MatlabStructArray
        length(cres) > 0 || throw(ArgumentError("$path_liferes: cres is empty"))
        cres[1]
    elseif cres isa AbstractDict
        cres
    elseif cres isa AbstractArray && !isempty(cres)
        first(vec(cres))
    else
        throw(ArgumentError("$path_liferes: cres must contain at least one struct"))
    end
    ep = mat_field_value(entry, "ep", path_liferes * " cres[1]")
    cam = mat_field_value(ep, "cam", path_liferes * " cres[1].ep")
    xscale = mat_field_value(cam, "xscale", path_liferes * " cres[1].ep.cam")
    while xscale isa AbstractArray
        length(xscale) == 1 ||
            throw(ArgumentError("$path_liferes: cres[1].ep.cam.xscale must be a scalar"))
        xscale = only(xscale)
    end
    xscale isa Real && isfinite(xscale) && xscale > 0 ||
        throw(ArgumentError("$path_liferes: cres[1].ep.cam.xscale must be finite and positive, got $xscale"))
    Float64(xscale)
end

function draw_insitu_scale_marker!(ax, image::AbstractMatrix, length_px::Real)
    height, width = size(image)
    margin_x = max(2.0, 0.05 * width)
    length_px <= width - 2margin_x ||
        throw(ArgumentError("scale marker length $(length_px) px does not fit image width $width px"))
    x_start = margin_x
    x_stop = x_start + length_px
    y_marker = height - max(2.0, 0.06 * height)
    lines!(ax, [x_start, x_stop], [y_marker, y_marker];
        color=:black, linewidth=5, inspectable=false)
    lines!(ax, [x_start, x_stop], [y_marker, y_marker];
        color=:white, linewidth=2, inspectable=false)
    nothing
end

function gaussian_blur(image::AbstractMatrix{<:Real}, sigma_px::Real)
    isfinite(sigma_px) && sigma_px >= 0 ||
        throw(ArgumentError("blur sigma must be finite and nonnegative, got $sigma_px"))
    sigma_px == 0 && return Float64.(image)
    radius = max(1, ceil(Int, 3sigma_px))
    offsets = collect(-radius:radius)
    kernel = exp.(-0.5 .* (offsets ./ sigma_px).^2)
    kernel ./= sum(kernel)

    height, width = size(image)
    source = Float64.(image)
    horizontal = similar(source)
    blurred = similar(source)
    for y in 1:height, x in 1:width
        horizontal[y, x] = sum(kernel[k + radius + 1] * source[y, clamp(x + k, 1, width)] for k in -radius:radius)
    end
    for y in 1:height, x in 1:width
        blurred[y, x] = sum(kernel[k + radius + 1] * horizontal[clamp(y + k, 1, height), x] for k in -radius:radius)
    end
    blurred
end

function make_insitu_sample_figures(tag::AbstractString, block, path_folder::AbstractString, path_output::AbstractString, blur_sigma_px::Real)
    name_acq = Symbol.(block["vars"])
    Set(name_acq) == Set((:f_626_MZM, :rep, :loadcfg, :istp)) && length(name_acq) == 4 ||
        throw(ArgumentError("$tag: vars must contain f_626_MZM, rep, loadcfg, and istp once each"))
    val_vars = Dict(name => block[string(name)] for name in name_acq)
    all(values -> values isa AbstractVector && !isempty(values), values(val_vars)) ||
        throw(ArgumentError("$tag: every configured variable must have a nonempty value list"))
    length(val_vars[:loadcfg]) == 1 ||
        throw(ArgumentError("$tag: expected one load configuration"))
    string.(val_vars[:istp]) == ["162", "164"] ||
        throw(ArgumentError("$tag: isotopes must be ordered 162, 164"))

    sources = block["source"]
    sources isa AbstractVector && length(sources) == 1 ||
        throw(ArgumentError("$tag: expected one liferes source file"))
    source_liferes = only(sources)
    source_allimg = replace(source_liferes, r"^liferes(?= )" => "allimg")
    source_allimg != source_liferes ||
        throw(ArgumentError("$tag: source must be named 'liferes <run>.mat'"))
    path_liferes = joinpath(path_folder, source_liferes)
    xscale = read_insitu_xscale(path_liferes)
    tag_upper = uppercase(tag)
    unit_length_m = occursin("CMOT", tag_upper) ? 1e-3 :
        occursin("MOT", tag_upper) ? 5e-3 :
        throw(ArgumentError("$tag: tag must identify MOT or CMOT for its scale marker"))
    scale_length_px = unit_length_m / xscale

    path_allimg = joinpath(path_folder, source_allimg)
    images = vec(matread(path_allimg)["allimg"])
    n_dims_acq = length.(getindex.(Ref(val_vars), name_acq))
    length(images) == prod(n_dims_acq) ||
        throw(DimensionMismatch("$tag: $(length(images)) images, expected $(prod(n_dims_acq)) from config.yaml"))
    all(image -> image isa AbstractMatrix{<:Real}, images) ||
        throw(ArgumentError("$tag: allimg cells must be numeric matrices"))
    size_image = size(first(images))
    all(image -> size(image) == size_image && all(isfinite, image), images) ||
        throw(ArgumentError("$tag: allimg cells must have equal dimensions and finite values"))

    # Acquisition variables run slowest to fastest in config.yaml.
    images_acq = reshape(images, Tuple(reverse(n_dims_acq))) |>
        data -> permutedims(data, reverse(1:length(n_dims_acq)))
    name = (:f_626_MZM, :rep, :loadcfg, :istp)
    images_fmt = permutedims(images_acq, indexin(collect(name), collect(name_acq)))
    val_f = val_vars[:f_626_MZM]
    val_rep = val_vars[:rep]
    max_image = maximum(maximum, images)
    max_image > 0 || throw(ArgumentError("$tag: all images are nonpositive"))

    fig = Figure(size=(1450, 2800), fontsize=17, backgroundcolor=:white)
    Label(fig[1, 1], "f_626_MZM")
    for (idx_rep, rep) in enumerate(val_rep)
        col_left = 2idx_rep
        Label(fig[1, col_left:col_left+1], "rep = $rep", fontsize=20)
        Label(fig[2, col_left], "162", fontsize=17)
        Label(fig[2, col_left+1], "164", fontsize=17)
    end

    hue_istp = (21, 259)
    alpha_base = 0.0
    thres_alpha = 0.05
    clrmaps = map(hue -> gen_clrmap_solo(hue; alpha_base, thres_alpha, chroma_max=0.18), hue_istp)
    for (idx_f, f) in enumerate(val_f)
        row = idx_f + 2
        Label(fig[row, 1], @sprintf("%.4f", f), fontsize=16)
        for idx_rep in eachindex(val_rep), idx_istp in 1:2
            col = 2idx_rep + idx_istp - 1
            ax = Axis(fig[row, col]; aspect=DataAspect(), yreversed=true)
            hidedecorations!(ax)
            hidespines!(ax)
            image = images_fmt[idx_f, idx_rep, 1, idx_istp]'
            heatmap!(ax, image;
                colormap=clrmaps[idx_istp], colorrange=(0, max_image), rasterize=true)
            draw_insitu_scale_marker!(ax, image, scale_length_px)
        end
    end
    for idx_rep in 1:(length(val_rep) - 1)
        colgap!(fig.layout, 2idx_rep + 1, 18)
    end
    colsize!(fig.layout, 1, Fixed(125))
    for col in 2:2length(val_rep)+1
        colsize!(fig.layout, col, Fixed(180))
    end
    rowsize!(fig.layout, 1, Fixed(34))
    rowsize!(fig.layout, 2, Fixed(30))
    for row in 3:length(val_f)+2
        rowsize!(fig.layout, row, Fixed(180))
    end
    path_plot = joinpath(path_output, "$tag samples.svg")
    save(path_plot, fig)
    println("Saved $path_plot ($(length(images)) images; shared color range 0–$max_image; scale marker $(unit_length_m * 1000) mm)")

    fig_misc = Figure(size=(850, 2800), fontsize=17, backgroundcolor=:white)
    Label(fig_misc[1, 1], "f_626_MZM")
    for (idx_rep, rep) in enumerate(val_rep)
        col = idx_rep + 1
        Label(fig_misc[1, col], "rep = $rep", fontsize=20)
        Label(fig_misc[2, col], "162 + 164", fontsize=17)
    end
    for (idx_f, f) in enumerate(val_f)
        row = idx_f + 2
        Label(fig_misc[row, 1], @sprintf("%.4f", f), fontsize=16)
        for idx_rep in eachindex(val_rep)
            ax = Axis(fig_misc[row, idx_rep+1]; aspect=DataAspect(), yreversed=true)
            hidedecorations!(ax)
            hidespines!(ax)
            clr_misc = to_miscibility_clr(
                images_fmt[idx_f, idx_rep, 1, 1],
                images_fmt[idx_f, idx_rep, 1, 2],
                hue_istp...;
                max=max_image,
                alpha_base,
                thres_alpha,
            )
            heatmap!(ax, clr_misc'; rasterize=true)
            draw_insitu_scale_marker!(ax, clr_misc', scale_length_px)
        end
    end
    colsize!(fig_misc.layout, 1, Fixed(125))
    for col in 2:length(val_rep)+1
        colsize!(fig_misc.layout, col, Fixed(180))
    end
    rowsize!(fig_misc.layout, 1, Fixed(34))
    rowsize!(fig_misc.layout, 2, Fixed(30))
    for row in 3:length(val_f)+2
        rowsize!(fig_misc.layout, row, Fixed(180))
    end
    path_plot_misc = joinpath(path_output, "$tag samples miscibility.svg")
    save(path_plot_misc, fig_misc)
    println("Saved $path_plot_misc ($(length(val_f) * length(val_rep)) paired images)")

    idx_f_selected = findfirst(f -> isapprox(f, 481.625; atol=1e-10, rtol=0), val_f)
    idx_rep_selected = findfirst(rep -> rep == 1, val_rep)
    idx_f_selected === nothing && throw(ArgumentError("$tag: f_626_MZM = 481.625 is not configured"))
    idx_rep_selected === nothing && throw(ArgumentError("$tag: rep = 1 is not configured"))

    fig_selected = Figure(size=(520, 330), fontsize=17, backgroundcolor=:white)
    Label(fig_selected[1, 1:2], "$tag — f_626_MZM = 481.625, rep = 1 (Gaussian blur σ = $(blur_sigma_px) px)", fontsize=18)
    for idx_istp in 1:2
        ax = Axis(fig_selected[3, idx_istp]; aspect=DataAspect(), yreversed=true)
        hidedecorations!(ax)
        hidespines!(ax)
        Label(fig_selected[2, idx_istp], string(val_vars[:istp][idx_istp]), fontsize=16)
        image_raw = images_fmt[idx_f_selected, idx_rep_selected, 1, idx_istp]
        image = gaussian_blur(image_raw, blur_sigma_px)'
        heatmap!(ax, image;
            colormap=clrmaps[idx_istp], colorrange=(0, max_image), rasterize=true)
        width = size(image, 2)
        height = size(image, 1)
        x_center = (width + 1) / 2
        lines!(ax, [x_center, x_center], [0.5, height + 0.5];
            color=:black, linewidth=2, inspectable=false)
        draw_insitu_scale_marker!(ax, image, scale_length_px)
    end
    colsize!(fig_selected.layout, 1, Fixed(230))
    colsize!(fig_selected.layout, 2, Fixed(230))
    path_plot_selected = joinpath(path_output, "$tag sample f_626_MZM 481.625 rep 1.svg")
    save(path_plot_selected, fig_selected)
    println("Saved $path_plot_selected (Gaussian blur σ = $(blur_sigma_px) px; rasterization enabled)")
end

path_folder = length(ARGS) >= 1 ? ARGS[1] : joinpath(path_root, "insitu samples", "162-164")
path_output = length(ARGS) >= 2 ? ARGS[2] : path_folder
config = YAML.load_file(joinpath(path_folder, "config.yaml"))
config isa AbstractVector && !isempty(config) ||
    throw(ArgumentError("expected one or more processing entries in config.yaml"))
tags = getindex.(config, "tag")
all(tag -> tag isa AbstractString && !isempty(strip(tag)), tags) ||
    throw(ArgumentError("every processing entry must have a nonempty tag"))
length(unique(tags)) == length(tags) ||
    throw(ArgumentError("processing tags in config.yaml must be unique"))
blur_sigma_px = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0
isfinite(blur_sigma_px) && blur_sigma_px >= 0 ||
    throw(ArgumentError("blur sigma must be finite and nonnegative, got $blur_sigma_px"))
mkpath(path_output)

for processing in config
    blocks = get(processing, "data", nothing)
    blocks isa AbstractVector && length(blocks) == 1 ||
        throw(ArgumentError("$(processing["tag"]): expected one rectangular data block"))
    make_insitu_sample_figures(processing["tag"], only(blocks), path_folder, path_output, blur_sigma_px)
end
