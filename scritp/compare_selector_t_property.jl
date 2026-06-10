## Run this after data is generated
for t_min in 0:10:100, t_max in 20:20:240
    if t_max - t_min <= 30
        continue
    end
    tag_range = @sprintf("%d-%d", t_min, t_max)
    selector_t_common = t -> t_min .< t .< t_max
    selector_t_spectrum = (;
        number=selector_t_common,
        sp_weight=selector_t_common,
        sp_height=selector_t_common,
        sp_width=selector_t_common,
        sp_wavenum=selector_t_common,
        nvlp=selector_t_common,
    )
    t_stage = log_step("analyzing per-shot trends")
    trend_sidepeak_nvlp = [
        extr_fmt[c, r, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_spectrum, query_weight_kwargs)
        for c in axes(extr_fmt, 1), r in axes(extr_fmt, 2), i in axes(extr_fmt, 4)
    ]
    log_done("analyzed per-shot trends", t_stage)

    t_stage = log_step("analyzing stacked trends")
    trend_extr_stacked_over_rep = [
        extr_stacked_over_rep[c, :, i] |> e -> anlz_trend_from_extr(val_vars.t_hold, e, freq_query; selector_t_spectrum, query_weight_kwargs)
        for c in axes(extr_stacked_over_rep, 1), i in axes(extr_stacked_over_rep, 3)
    ]
    log_done("analyzed stacked trends", t_stage)

    trend_stacked_over_rep = [
        trend_sidepeak_nvlp[c, :, i] |> mean_dict
        for c in axes(trend_sidepeak_nvlp, 1), i in axes(trend_sidepeak_nvlp, 3)
    ]

    t_stage = log_step("building and saving spectrum-vs-IB figures")
    trend_spectrum_groups = (;
        stacked=trend_extr_stacked_over_rep,
        all=trend_stacked_over_rep,
    )
    for spec in trend_property_specs
        fig_spectrum, axs_spectrum = set_axis_spectrum_property_IB_istp!(
            val_vars.IB,
            val_vars.istp,
            spec,
            "$tag | spectrum vs IB | $(spec.name)";
            groups=trend_spectrum_IB_groups,
            trend_spectrum_IB_kwargs...,
        )
        plot_spectrum_property_IB_istp!(
            axs_spectrum,
            trend_spectrum_groups,
            val_vars.IB,
            val_vars.istp,
            spec;
            trend_spectrum_IB_plot_kwargs...,
        )
        resize_to_layout!(fig_spectrum)
        fig_spectrum |> f -> save(
            joinpath(path_output, tag, @sprintf("%s_%s_spectrum_IB.png", spec.name, tag_range)),
            f;
            backend=CairoMakie,
            force=true,
        )
    end
    log_done("saved spectrum-vs-IB figures", t_stage)
end

## Build selector-range summary grids directly from the saved spectrum images.
using CairoMakie

function parse_selector_t_image(path_image::AbstractString)
    match_image = match(
        r"^(?<property>.+)_(?<t_min>\d+)-(?<t_max>\d+)_spectrum_IB\.png$",
        basename(path_image),
    )
    isnothing(match_image) && return nothing
    return (
        property=match_image[:property],
        t_min=parse(Int, match_image[:t_min]),
        t_max=parse(Int, match_image[:t_max]),
        path=path_image,
    )
end

function collect_selector_t_images(path_group::AbstractString)
    isdir(path_group) || throw(ArgumentError("selector image folder does not exist: $path_group"))
    records = filter(!isnothing, parse_selector_t_image.(readdir(path_group; join=true)))
    isempty(records) && throw(ArgumentError("no selector spectrum images found in $path_group"))
    return records
end

function save_selector_t_summary(
    path_group::AbstractString,
    property::AbstractString,
    records::AbstractVector;
    scale_image::Real=0.25,
)
    0 < scale_image <= 1 || throw(ArgumentError(
        "scale_image must be in (0, 1], got $scale_image",
    ))
    records_property = filter(record -> record.property == property, records)
    isempty(records_property) && throw(ArgumentError("no images found for property $property in $path_group"))

    val_t_min = sort!(unique(record.t_min for record in records_property))
    val_t_max = sort!(unique(record.t_max for record in records_property))
    path_by_range = Dict(
        (record.t_min, record.t_max) => record.path
        for record in records_property
    )
    size_image = first(records_property).path |> Makie.FileIO.load |> size
    length(size_image) == 2 || throw(DimensionMismatch(
        "expected a 2D image for $property, got size $size_image",
    ))
    image_height, image_width = size_image
    cell_height = round(Int, scale_image * image_height)
    cell_width = round(Int, scale_image * image_width)
    name_group = basename(path_group)
    println(
        "[$name_group/$property] building $(length(val_t_min)) x $(length(val_t_max)) grid " *
        "from $(length(records_property)) images; source $(image_width) x $(image_height) px, " *
        "grid cell $(cell_width) x $(cell_height) px ($(scale_image)x)",
    )
    size_title = 70
    size_axis_label = 55
    gap_cell = 4
    width_label = 140
    height_title = 110
    height_col_label = 90

    fig = Figure()
    Label(
        fig[1, 2:length(val_t_max) + 1],
        "$(basename(path_group)) | $property | selector range (ms)";
        fontsize=size_title,
        tellwidth=false,
    )
    Label(fig[2, 1], "t_min"; fontsize=size_axis_label)
    Label(fig[1, 1], "t_max"; fontsize=size_axis_label)

    for (idx_col, t_max) in enumerate(val_t_max)
        Label(fig[2, idx_col + 1], string(t_max); fontsize=size_axis_label)
    end
    for (idx_row, t_min) in enumerate(val_t_min)
        println("[$name_group/$property] row $idx_row/$(length(val_t_min)): t_min=$t_min ms")
        Label(fig[idx_row + 2, 1], string(t_min); fontsize=size_axis_label)
        for (idx_col, t_max) in enumerate(val_t_max)
            ax = Axis(
                fig[idx_row + 2, idx_col + 1];
                aspect=DataAspect(),
                yreversed=true,
            )
            hidedecorations!(ax)
            hidespines!(ax)
            path_image = get(path_by_range, (t_min, t_max), nothing)
            isnothing(path_image) && continue
            image_cell = Makie.FileIO.load(path_image)
            size(image_cell) == size_image || throw(DimensionMismatch(
                "property $property images must share one size; " *
                "$(basename(path_image)) has size $(size(image_cell)), expected $size_image",
            ))
            image!(ax, permutedims(image_cell, (2, 1)))
        end
    end

    colsize!(fig.layout, 1, Fixed(width_label))
    for idx_col in eachindex(val_t_max)
        colsize!(fig.layout, idx_col + 1, Fixed(cell_width))
    end
    rowsize!(fig.layout, 1, Fixed(height_title))
    rowsize!(fig.layout, 2, Fixed(height_col_label))
    for idx_row in eachindex(val_t_min)
        rowsize!(fig.layout, idx_row + 2, Fixed(cell_height))
    end
    colgap!(fig.layout, gap_cell)
    rowgap!(fig.layout, gap_cell)
    println("[$name_group/$property] resizing figure to layout")
    resize_to_layout!(fig)
    path_output = joinpath(path_group, "$(property)_selector_t_summary.pdf")
    println("[$name_group/$property] saving $path_output")
    save(path_output, fig; backend=CairoMakie, force=true)
    println("[$name_group/$property] saved")
    return path_output
end

function save_selector_t_summaries(path_root::AbstractString)
    paths_output = String[]
    for name_group in ("CFNM", "NTRC")
        path_group = joinpath(path_root, name_group)
        println("[$name_group] scanning $path_group")
        records = collect_selector_t_images(path_group)
        properties = sort!(unique(record.property for record in records))
        println("[$name_group] found $(length(records)) images across $(length(properties)) properties")
        append!(
            paths_output,
            save_selector_t_summary(path_group, property, records)
            for property in properties
        )
        println("[$name_group] completed $(length(properties)) property summaries")
    end
    println("completed $(length(paths_output)) selector summary PDFs")
    return paths_output
end

path_selector_t_root = path_output
paths_selector_t_summary = save_selector_t_summaries(path_selector_t_root)
println.(paths_selector_t_summary)
