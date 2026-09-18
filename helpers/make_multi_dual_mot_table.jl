using Base64

const DEFAULT_ROOT = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data"
const DEFAULT_OUTPUT_NAME = "multi_dual_mot_table.svg"
const INKSCAPE_CANDIDATES = [
    raw"C:\Program Files\Inkscape\bin\inkscape.exe",
    raw"C:\Program Files\Inkscape\inkscape.exe",
]

# Only these folders are included, in this outer-grid order.
const COLUMN_PARENT_NAMES = [
    "MOT loading 421",
    "MOT loading 626",
    "MOT lifetime",
    "CMOT lifetime",
    "ODT BField",
]
const ROW_PAIR_NAMES = [
    "162-164",
    "160-162",
    "161-162",
    "161-164",
    "163-164",
    "162-163",
    "161-163",
]

const INNER_GAP = 2.0
const OUTER_GAP = 6.0
const ROW_LABEL_WIDTH = 150.0
const COLUMN_LABEL_HEIGHT = 34.0
const PAGE_PADDING = 8.0
const PNG_SCALE = 2.0
const PNG_RENDER_TIMEOUT = 60.0

struct SvgEntry
    path::String
    test_tag::String
    style_tag::String
    subvariant_tag::String
    width::Float64
    height::Float64
end

function escape_xml(text::AbstractString)
    replace(text, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;", '\'' => "&apos;")
end

function parse_svg_length(value::AbstractString, path::AbstractString)
    matched = match(r"^\s*([0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:px|pt|pc|mm|cm|in)?\s*$"i, value)
    isnothing(matched) && throw(ArgumentError("unsupported SVG length '$value' in $path"))
    parse(Float64, matched.captures[1])
end

function read_svg_size(path::AbstractString)
    text = read(path, String)
    svg_tag = match(r"<svg\b[^>]*>"is, text)
    isnothing(svg_tag) && throw(ArgumentError("no <svg> root element found in $path"))
    tag = svg_tag.match

    width_match = match(r"\bwidth\s*=\s*['\"]([^'\"]+)['\"]"i, tag)
    height_match = match(r"\bheight\s*=\s*['\"]([^'\"]+)['\"]"i, tag)
    if !isnothing(width_match) && !isnothing(height_match)
        return (
            parse_svg_length(width_match.captures[1], path),
            parse_svg_length(height_match.captures[1], path),
        )
    end

    viewbox_match = match(r"\bviewBox\s*=\s*['\"]([^'\"]+)['\"]"i, tag)
    isnothing(viewbox_match) && throw(ArgumentError("SVG needs width/height or viewBox attributes: $path"))
    values = split(strip(viewbox_match.captures[1]), r"[\s,]+") .|> x -> parse(Float64, x)
    length(values) == 4 || throw(ArgumentError("viewBox must contain four numbers in $path"))
    values[3] > 0 && values[4] > 0 || throw(ArgumentError("viewBox size must be positive in $path"))
    (values[3], values[4])
end

function parse_svg_entry(path::AbstractString)
    matched = match(r"^\[(.+?)\]\.\[(.+?)\]\.\[(.+?)\]\.svg$"i, basename(path))
    isnothing(matched) && throw(ArgumentError(
        "SVG filename must be [test tag].[style tag].[subvariant tag].svg: $path",
    ))
    width, height = read_svg_size(path)
    SvgEntry(path, matched.captures..., width, height)
end

function selected_parent_dirs(root::AbstractString,
                              parent_names::AbstractVector{<:AbstractString})
    available = Dict(basename(path) => path for path in readdir(root; join=true) if isdir(path))
    missing = filter(name -> !haskey(available, name), parent_names)
    isempty(missing) || throw(ArgumentError(
        "configured parent folders do not exist in $root: $(join(missing, ", "))",
    ))
    [available[name] for name in parent_names]
end

function collect_entries(parent_dirs::AbstractVector{<:AbstractString},
                         pair_names::AbstractVector{<:AbstractString})
    allowed_pairs = Set(pair_names)
    entries = Dict{Tuple{String,String},Vector{SvgEntry}}()

    for parent in parent_dirs
        parent_name = basename(parent)
        for pair_dir in filter(isdir, readdir(parent; join=true))
            pair_name = basename(pair_dir)
            pair_name in allowed_pairs || continue
            svg_paths = filter(path -> isfile(path) && endswith(lowercase(path), ".svg"),
                               readdir(pair_dir; join=true))
            cell_entries = parse_svg_entry.(sort(svg_paths; by=lowercase))
            isempty(cell_entries) || (entries[(parent_name, pair_name)] = cell_entries)
        end
    end
    entries
end

function cell_layout(cell_entries::AbstractVector{SvgEntry})
    isempty(cell_entries) && return (
        styles=String[], subvariants=String[], widths=Float64[], heights=Float64[],
        width=0.0, height=0.0,
    )

    styles = sort!(unique(entry.style_tag for entry in cell_entries); by=lowercase)
    subvariants = sort!(unique(entry.subvariant_tag for entry in cell_entries); by=lowercase)
    widths = [maximum((entry.width for entry in cell_entries if entry.style_tag == tag);
                      init=0.0) for tag in styles]
    heights = [maximum((entry.height for entry in cell_entries if entry.subvariant_tag == tag);
                       init=0.0) for tag in subvariants]
    width = sum(widths) + INNER_GAP * max(length(widths) - 1, 0)
    height = sum(heights) + INNER_GAP * max(length(heights) - 1, 0)
    (; styles, subvariants, widths, heights, width, height)
end

function offsets(sizes::AbstractVector{<:Real}, gap::Real)
    result = zeros(Float64, length(sizes))
    for index in 2:length(sizes)
        result[index] = result[index - 1] + sizes[index - 1] + gap
    end
    result
end

function svg_image_element(io::IO, entry::SvgEntry, x::Real, y::Real)
    encoded = base64encode(read(entry.path))
    title = escape_xml("$(entry.test_tag) / $(entry.style_tag) / $(entry.subvariant_tag)")
    println(io, "<image x=\"$x\" y=\"$y\" width=\"$(entry.width)\" height=\"$(entry.height)\"",
            " preserveAspectRatio=\"xMidYMid meet\" href=\"data:image/svg+xml;base64,$encoded\">",
            "<title>$title</title></image>")
end

function read_png_size(path::AbstractString)
    header = open(path, "r") do io
        read(io, 24)
    end
    length(header) == 24 || throw(ArgumentError("truncated PNG header: $path"))
    header[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a] ||
        throw(ArgumentError("renderer output is not a PNG: $path"))
    decode_u32(bytes) = foldl((value, byte) -> (value << 8) | UInt32(byte), bytes;
                              init=UInt32(0))
    (Int(decode_u32(@view header[17:20])), Int(decode_u32(@view header[21:24])))
end


function render_svg_png(svg_path::AbstractString, png_path::AbstractString,
                        width::Real, height::Real; scale::Real=PNG_SCALE,
                        max_attempts::Integer=3)
    isfinite(scale) && scale > 0 || throw(ArgumentError(
        "PNG scale must be positive and finite; received $scale",
    ))
    max_attempts > 0 || throw(ArgumentError("max_attempts must be positive; received $max_attempts"))
    renderer_index = findfirst(isfile, INKSCAPE_CANDIDATES)
    isnothing(renderer_index) && throw(ArgumentError(
        "PNG output requires Inkscape; no Inkscape executable was found",
    ))
    renderer = INKSCAPE_CANDIDATES[renderer_index]
    expected_size = (ceil(Int, scale * width), ceil(Int, scale * height))
    destination = abspath(png_path)
    mkpath(dirname(destination))
    diagnostics = ""
    last_error = nothing
    succeeded = false

    for attempt in 1:max_attempts
        try
            mktempdir() do temp_dir
                temp_png = joinpath(temp_dir, "render.png")
                log_path = joinpath(temp_dir, "renderer.log")
                command = Cmd([
                    renderer,
                    abspath(svg_path),
                    "--export-type=png",
                    "--export-area-page",
                    "--export-filename=$temp_png",
                    "--export-width=$(expected_size[1])",
                ])
                try
                    open(log_path, "w") do log_io
                        process = run(pipeline(command; stdout=log_io, stderr=log_io);
                                      wait=false)
                        wait_status = timedwait(() -> process_exited(process),
                                                PNG_RENDER_TIMEOUT; pollint=0.25)
                        if wait_status == :timed_out
                            kill(process)
                            timedwait(() -> process_exited(process), 5.0; pollint=0.1)
                            error("renderer timed out after $(PNG_RENDER_TIMEOUT) seconds")
                        end
                        success(process) || error("renderer exited unsuccessfully")
                    end
                    isfile(temp_png) && filesize(temp_png) > 0 ||
                        error("renderer did not create a nonempty PNG")
                    actual_size = read_png_size(temp_png)
                    actual_size == expected_size || error(
                        "renderer created $(actual_size[1]) × $(actual_size[2]), expected " *
                        "$(expected_size[1]) × $(expected_size[2])",
                    )
                    mv(temp_png, destination; force=true)
                    succeeded = true
                catch error_render
                    last_error = error_render
                    diagnostics = isfile(log_path) ? read(log_path, String) : ""
                end
            end
        catch error_temp
            last_error = error_temp
        end
        succeeded && break
        if attempt < max_attempts
            reason = isnothing(last_error) ? "unknown renderer error" : sprint(showerror, last_error)
            println(stderr, "PNG render attempt $attempt failed ($reason); " *
                            "retrying ($max_attempts attempts maximum).")
        end
    end

    if succeeded
        println("Wrote $(expected_size[1]) × $(expected_size[2]) PNG rendering to:")
        println(destination)
        return destination
    end

    detail = sprint(showerror, last_error)
    isempty(diagnostics) || (detail *= "\nRenderer output:\n" * last(diagnostics, min(2000, length(diagnostics))))
    error("PNG rendering failed after $max_attempts attempts: $detail")
end

function make_multi_dual_mot_table(root::AbstractString, output_path::AbstractString;
                                   parent_names::AbstractVector{<:AbstractString}=COLUMN_PARENT_NAMES,
                                   pair_names::AbstractVector{<:AbstractString}=ROW_PAIR_NAMES,
                                   png_path::Union{Nothing,AbstractString}=nothing)
    isdir(root) || throw(ArgumentError("data root does not exist: $root"))
    isempty(parent_names) && throw(ArgumentError("parent_names cannot be empty"))
    isempty(pair_names) && throw(ArgumentError("pair_names cannot be empty"))
    allunique(parent_names) || throw(ArgumentError("parent_names must be unique: $parent_names"))
    allunique(pair_names) || throw(ArgumentError("pair_names must be unique: $pair_names"))
    parent_dirs = selected_parent_dirs(root, parent_names)
    entries = collect_entries(parent_dirs, pair_names)
    isempty(entries) && throw(ArgumentError("no SVG images matched the configured rows and columns"))

    layouts = Dict(
        (parent_name, pair_name) => cell_layout(
            get(entries, (parent_name, pair_name), SvgEntry[]),
        )
        for parent_name in parent_names for pair_name in pair_names
    )
    column_widths = [maximum(
        (layouts[(parent_name, pair_name)].width for pair_name in pair_names); init=0.0,
    ) for parent_name in parent_names]
    row_heights = [maximum(
        (layouts[(parent_name, pair_name)].height for parent_name in parent_names); init=0.0,
    ) for pair_name in pair_names]
    column_offsets = offsets(column_widths, OUTER_GAP)
    row_offsets = offsets(row_heights, OUTER_GAP)
    content_width = sum(column_widths) + OUTER_GAP * max(length(column_widths) - 1, 0)
    content_height = sum(row_heights) + OUTER_GAP * max(length(row_heights) - 1, 0)
    canvas_width = 2PAGE_PADDING + ROW_LABEL_WIDTH + content_width
    canvas_height = 2PAGE_PADDING + COLUMN_LABEL_HEIGHT + content_height
    x_content = PAGE_PADDING + ROW_LABEL_WIDTH
    y_content = PAGE_PADDING + COLUMN_LABEL_HEIGHT

    mkpath(dirname(abspath(output_path)))
    open(output_path, "w") do io
        println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$canvas_width\" height=\"$canvas_height\" viewBox=\"0 0 $canvas_width $canvas_height\">")
        println(io, "<title>Dual-isotope MOT figure table</title>")
        println(io, "<desc>Columns are analysis folders, rows are isotope-pair folders, and each cell uses style-tag columns and subvariant-tag rows.</desc>")
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
        println(io, "<g font-family=\"Arial, Helvetica, sans-serif\" fill=\"#202020\">")

        for (column_index, parent_name) in enumerate(parent_names)
            x = x_content + column_offsets[column_index]
            center_x = x + column_widths[column_index] / 2
            println(io, "<text x=\"$center_x\" y=\"$(PAGE_PADDING + 20)\" text-anchor=\"middle\" font-size=\"18\" font-weight=\"bold\">$(escape_xml(parent_name))</text>")
        end

        for (row_index, pair_name) in enumerate(pair_names)
            y = y_content + row_offsets[row_index]
            center_y = y + row_heights[row_index] / 2
            println(io, "<text x=\"$(x_content - 10)\" y=\"$center_y\" text-anchor=\"end\" dominant-baseline=\"middle\" font-size=\"18\" font-weight=\"bold\">$(escape_xml(pair_name))</text>")

            for (column_index, parent_name) in enumerate(parent_names)
                layout = layouts[(parent_name, pair_name)]
                x_cell = x_content + column_offsets[column_index]
                y_cell = y
                println(io, "<rect x=\"$x_cell\" y=\"$y\" width=\"$(column_widths[column_index])\" height=\"$(row_heights[row_index])\" fill=\"none\" stroke=\"#b8b8b8\" stroke-width=\"1\"/>")

                style_offsets = offsets(layout.widths, INNER_GAP)
                subvariant_offsets = offsets(layout.heights, INNER_GAP)
                cell_entries = get(entries, (parent_name, pair_name), SvgEntry[])
                occupied = Set{Tuple{String,String}}()
                for entry in cell_entries
                    position = (entry.style_tag, entry.subvariant_tag)
                    position in occupied && throw(ArgumentError(
                        "duplicate style/subvariant position $(position) in $parent_name/$pair_name",
                    ))
                    push!(occupied, position)
                    style_index = findfirst(==(entry.style_tag), layout.styles)
                    subvariant_index = findfirst(==(entry.subvariant_tag), layout.subvariants)
                    slot_width = layout.widths[style_index]
                    slot_height = layout.heights[subvariant_index]
                    image_x = x_cell + style_offsets[style_index] +
                              (slot_width - entry.width) / 2
                    image_y = y_cell + subvariant_offsets[subvariant_index] +
                              (slot_height - entry.height) / 2
                    svg_image_element(io, entry, image_x, image_y)
                end
            end
        end

        println(io, "</g>")
        println(io, "</svg>")
    end

    println("Wrote $(length(parent_names))-column × $(length(pair_names))-row SVG table to:")
    println(abspath(output_path))
    isnothing(png_path) || render_svg_png(output_path, png_path, canvas_width, canvas_height)
    abspath(output_path)
end

function main(args::AbstractVector{<:AbstractString}=ARGS)
    length(args) <= 3 || throw(ArgumentError(
        "usage: julia helpers/make_multi_dual_mot_table.jl [data_root] [output.svg] [output.png]",
    ))
    root = isempty(args) ? DEFAULT_ROOT : abspath(args[1])
    output_path = length(args) < 2 ? joinpath(root, DEFAULT_OUTPUT_NAME) : abspath(args[2])
    png_path = length(args) < 3 ? splitext(output_path)[1] * ".png" : abspath(args[3])
    make_multi_dual_mot_table(root, output_path; png_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__) || isinteractive()
    main()
end
