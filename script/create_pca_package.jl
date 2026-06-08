using JLD2
using Printf

include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
title_anlz = "[06.08].64.PCATests.[t=0-100ms]"
tag = "CFNM"

path_anlz = joinpath(path_root, "AnlzRoutine", title_anlz)
path_pca = joinpath(path_anlz, "PCA modes")
path_pca_data = joinpath(path_anlz, @sprintf("%s_pca_modes.jld2", tag))
path_package = joinpath(path_anlz, "PCA modes PicTaggerCache")
path_cache = joinpath(path_package, "cache.yaml")

@load path_pca_data runinfo val_vars n_pca_modes pca_spectra

get_bind_date(runinfo, idx_bind) = hasproperty(runinfo, :date_runid) ? first(runinfo.date_runid[idx_bind]) : runinfo.date
get_bind_runid(runinfo, idx_bind) =
    if hasproperty(runinfo, :date_runid)
        last(runinfo.date_runid[idx_bind])
    elseif hasproperty(runinfo, :runids)
        as_vector(runinfo.runids)[idx_bind]
    else
        as_vector(runinfo.runid)[idx_bind]
    end
get_bind_runinfo(runinfo, val_vars, idx_bind) = merge(
    runinfo,
    (;
        date=get_bind_date(runinfo, idx_bind),
        runid=get_bind_runid(runinfo, idx_bind),
        IB=val_vars.IB[idx_bind],
    ),
)

const TAG_RECORDS = [
    (; name="overall", hue=0.0),
    (; name="dipole x", hue=39.0),
    (; name="dipole y", hue=100.0),
    (; name="quad", hue=128.0),
    (; name="caterpillar", hue=153.0),
    (; name="breathing", hue=188.0),
    (; name="strip modulation", hue=206.0),
    (; name="modulation", hue=268.0),
    (; name="unkown", hue=346.0),
    (; name="phase/crystal mode", hue=303.0),
    (; name="array dipole x", hue=69.0),
    (; name="dipole inclined", hue=23.0),
    (; name="complex 2d pattern", hue=242.0),
]

yaml_value(::Nothing) = "null"
yaml_value(x::AbstractString) = x
yaml_value(x::Real) = string(float(x))

function write_tag_records(io, tag_records)
    println(io, "tags:")
    for tag_record in tag_records
        println(io, "- name: $(yaml_value(tag_record.name))")
        println(io, "  hue: $(yaml_value(tag_record.hue))")
    end
end

function write_freq_weight_pairs(io, peaks)
    println(io, "  freq_weight_pairs:")
    if isempty(peaks)
        println(io, "  - frequency: null")
        println(io, "    weight: null")
    else
        for peak in peaks
            println(io, "  - frequency: $(yaml_value(peak.freq))")
            println(io, "    weight: $(yaml_value(peak.value_reduced))")
        end
    end
end

function get_peaks_prominent(pca_spectrum)
    if hasproperty(pca_spectrum, :peaks_prominent)
        return pca_spectrum.peaks_prominent
    elseif pca_spectrum isa Tuple && length(pca_spectrum) >= 2
        return pca_spectrum[2]
    else
        throw(ArgumentError("PCA spectrum entry must contain peaks_prominent or be a tuple with peaks as its second entry; got $(typeof(pca_spectrum))."))
    end
end

function find_mode_image(path_pca, tag_IB::AbstractString, idx_mode::Integer)
    path_src = joinpath(path_pca, tag_IB, @sprintf("%s_%d.png", tag_IB, idx_mode))
    isfile(path_src) && return path_src

    path_src_flat = joinpath(path_pca, @sprintf("%s_%d.png", tag_IB, idx_mode))
    isfile(path_src_flat) && return path_src_flat

    throw(ArgumentError("Missing PCA mode image for $tag_IB mode $idx_mode. Tried $path_src and $path_src_flat."))
end

function write_image_record(io, record)
    println(io, "- id: null")
    println(io, "  image_path: $(record.image_path)")
    println(io, "  ib: $(yaml_value(record.ib))")
    println(io, "  source: $(record.source)")
    println(io, "  source_tag: $(record.source_tag)")
    println(io, "  tags: []")
    println(io, "  index: $(record.index)")
    write_freq_weight_pairs(io, record.peaks)
    println(io, "  created_at: null")
    println(io, "  updated_at: null")
end

function create_pca_package(; path_pca, path_package, path_cache, runinfo, val_vars, n_pca_modes, pca_spectra)
    isdir(path_pca) || throw(ArgumentError("PCA image folder does not exist: $path_pca"))
    mkpath(path_package)

    records = []
    idx_image = 0
    for (idx_ib, IB) in enumerate(val_vars.IB)
        tag_IB = gen_run_tag(get_bind_runinfo(runinfo, val_vars, idx_ib))
        path_dst_dir = joinpath(path_package, tag_IB)
        mkpath(path_dst_dir)

        for idx_mode in 1:n_pca_modes
            idx_image += 1
            filename_dst = @sprintf("PCA.Mode.%02d.png", idx_mode)
            path_src = find_mode_image(path_pca, tag_IB, idx_mode)
            path_dst = joinpath(path_dst_dir, filename_dst)
            cp(path_src, path_dst; force=true)

            push!(
                records,
                (;
                    image_path=joinpath(tag_IB, filename_dst) |> p -> replace(p, "\\" => "/"),
                    ib=IB,
                    source=filename_dst,
                    source_tag=tag_IB,
                    index=idx_mode,
                    peaks=get_peaks_prominent(pca_spectra[idx_ib][idx_mode]),
                ),
            )
        end
    end

    open(path_cache, "w") do io
        println(io, "version: 1")
        write_tag_records(io, TAG_RECORDS)
        println(io, "images:")
        for record in records
            write_image_record(io, record)
        end
    end

    return (; path_package, path_cache, n_images=idx_image, n_IB=length(val_vars.IB), n_pca_modes)
end

result = create_pca_package(; path_pca, path_package, path_cache, runinfo, val_vars, n_pca_modes, pca_spectra)
println("Created PicTagger PCA cache package at $(result.path_package)")
println("Wrote $(result.n_images) images across $(result.n_IB) IB folders with $(result.n_pca_modes) modes each.")
