using HDF5
using MAT
using Printf
using Statistics

as_vector(x) = x isa AbstractArray ? collect(x) : [x]
format_vars(vars::NamedTuple) = map(as_vector, vars)
is_selector_tuple(selector) = selector isa NamedTuple && any(n -> n in propertynames(selector), (:index, :val))
get_selector_val(selector) = is_selector_tuple(selector) ? (hasproperty(selector, :val) ? selector.val : Colon()) : selector
get_selector_index(selector) = is_selector_tuple(selector) && hasproperty(selector, :index) ? selector.index : Colon()
format_runids(runids::Integer) = @sprintf("%02d", runids)
format_runids(runids::AbstractRange{<:Integer}) = runids |> a -> "$(a)" |> s -> replace(s, ":" => "-")
format_runids(runids::AbstractVector{<:Integer}) = join(format_runids.(runids), "-")
format_runids(runids) = runids |> a -> "$(a)" |> s -> replace(s, ":" => "-")
get_date(date_runid::Tuple) = first(date_runid)
get_runid_from_date_runid(date_runid::Tuple) = last(date_runid)

function select_by_selector(vals::AbstractVector, selector, name_var::Symbol)
    idx_all = collect(eachindex(vals))
    selector_idx = get_selector_index(selector)
    selector_val = get_selector_val(selector)

    idx_from_index = if selector_idx isa Colon
        idx_all
    elseif selector_idx isa Function
        [i for i in idx_all if selector_idx(i)]
    elseif selector_idx isa AbstractVector{Bool}
        length(selector_idx) == length(vals) || throw(DimensionMismatch("sel_vars.$name_var.index has length $(length(selector_idx)); expected $(length(vals))."))
        idx_all[selector_idx]
    elseif selector_idx isa AbstractArray{<:Integer} || selector_idx isa AbstractRange{<:Integer}
        collect(selector_idx)
    elseif selector_idx isa Integer
        [selector_idx]
    else
        throw(ArgumentError("sel_vars.$name_var.index must be an index, index range, boolean mask, or predicate function; got $(typeof(selector_idx))."))
    end

    all((1 .<= idx_from_index) .& (idx_from_index .<= length(vals))) || throw(BoundsError(vals, idx_from_index))
    idx_sel = if selector_val isa Colon
        idx_from_index
    elseif selector_val isa Function
        [i for i in idx_from_index if selector_val(vals[i])]
    elseif selector_val isa AbstractVector{Bool}
        length(selector_val) == length(vals) || throw(DimensionMismatch("sel_vars.$name_var has length $(length(selector_val)); expected $(length(vals))."))
        [i for i in idx_from_index if selector_val[i]]
    elseif selector_val isa AbstractArray || selector_val isa AbstractRange
        [i for i in idx_from_index if vals[i] in selector_val]
    else
        [i for i in idx_from_index if vals[i] == selector_val]
    end

    isempty(idx_sel) && throw(ArgumentError("sel_vars.$name_var selected no values from $(vals)."))
    return idx_sel
end

function select_val_vars(val_vars::NamedTuple, sel_vars::NamedTuple)
    names_vars = propertynames(val_vars)
    names_sel = propertynames(sel_vars)
    unknown = setdiff(names_sel, names_vars)
    isempty(unknown) || throw(ArgumentError("sel_vars contains unknown variables $(unknown); expected variables are $(names_vars)."))

    idx_vars = map(names_vars) do name_var
        hasproperty(sel_vars, name_var) ? select_by_selector(getproperty(val_vars, name_var), getproperty(sel_vars, name_var), name_var) : collect(eachindex(getproperty(val_vars, name_var)))
    end
    val_vars_sel = NamedTuple{names_vars}(map((vals, idx) -> vals[idx], Tuple(val_vars), idx_vars))
    return val_vars_sel, Tuple(idx_vars)
end

select_val_vars(val_vars::NamedTuple, ::Nothing) = (val_vars, Tuple(map(vals -> collect(eachindex(vals)), Tuple(val_vars))))

function select_runinfo(runinfo, val_vars_sel::NamedTuple, val_vars_full::NamedTuple, idx_vars::Tuple)
    pairs = Any[:vars => val_vars_sel, :vars_full => val_vars_full]
    if hasproperty(runinfo, :bind_id)
        name_bound = runinfo.bind_id
        pos_bound = findfirst(==(name_bound), propertynames(val_vars_sel))
        if pos_bound !== nothing
            idx_bound = idx_vars[pos_bound]
            if hasproperty(runinfo, :date_runid)
                push!(pairs, :date_runid => as_vector(runinfo.date_runid)[idx_bound])
            elseif hasproperty(runinfo, :runids)
                push!(pairs, :runids => as_vector(runinfo.runids)[idx_bound])
            elseif hasproperty(runinfo, :runid) && length(as_vector(runinfo.runid)) > 1
                push!(pairs, :runid => as_vector(runinfo.runid)[idx_bound])
            end
        end
    end
    return merge(runinfo, (; pairs...))
end

function get_n_runids(runinfo)
    if hasproperty(runinfo, :runids)
        return length(as_vector(runinfo.runids))
    elseif hasproperty(runinfo, :runid)
        return length(as_vector(runinfo.runid))
    elseif hasproperty(runinfo, :date_runid)
        return length(as_vector(runinfo.date_runid))
    else
        error("runinfo must contain runid, runids, or date_runid.")
    end
end

function validate_bind_id(runinfo, val_vars)
    hasproperty(runinfo, :bind_id) || return nothing
    name_bound = runinfo.bind_id
    hasproperty(val_vars, name_bound) || throw(ArgumentError("runinfo.bind_id is $name_bound, but runinfo.vars does not contain that variable."))

    n_runids = get_n_runids(runinfo)
    n_bound = length(getproperty(val_vars, name_bound))
    n_runids == n_bound || throw(DimensionMismatch("runinfo.bind_id=$name_bound requires the same number of runids and $(name_bound) values; got $n_runids runids and $n_bound $(name_bound) values."))
    return nothing
end

function gen_run_tag(runinfo)
    runids = if hasproperty(runinfo, :runids)
        runinfo.runids
    elseif hasproperty(runinfo, :runid)
        runinfo.runid
    elseif hasproperty(runinfo, :date_runid)
        get_runid_from_date_runid.(as_vector(runinfo.date_runid))
    else
        error("runinfo must contain runid, runids, or date_runid.")
    end
    str_runids = format_runids(runids)
    if hasproperty(runinfo, :IB)
        return @sprintf("%s_%.3f_r%s", runinfo.tag_head, runinfo.IB, str_runids)
    end
    if hasproperty(runinfo, :vars) && hasproperty(runinfo.vars, :IB)
        val_ib = as_vector(runinfo.vars.IB)
        if length(val_ib) == 1
            return @sprintf("%s_%.3f_r%s", runinfo.tag_head, only(val_ib), str_runids)
        end
    end
    return @sprintf("%s_run%s", runinfo.tag_head, str_runids)
end

function load_dens_run(
    date::AbstractString,
    runid::Integer;
    path_root::AbstractString,
    year_test::Integer,
    wh_corner::Tuple{<:Integer,<:Integer},
)
    dir_data = gen_date_path(date, year_test)
    file_data = gen_h5name(date, runid)
    path_input_routine = joinpath(path_root, dir_data, @sprintf("run%02d", runid), file_data)
    path_input_flat = joinpath(path_root, @sprintf("run%02d", runid), file_data)
    path_input = isfile(path_input_routine) ? path_input_routine : path_input_flat

    h5open(path_input, "r") do f
        dens_run = f["/od"] |>
                   read |>
                   x_vec -> permutedims(x_vec, (3, 2, 1)) |>
                            x_vec -> stack(
                       map(d -> subtract_corner_mean(d, wh_corner), eachslice(x_vec; dims=1));
                       dims=1,
                   )
        ndims(dens_run) == 3 || error("Expected /od in $path_input to have 3 dimensions after formatting, got $(ndims(dens_run)).")
        return dens_run
    end
end

function format_image_array(dens_crop::AbstractArray{<:Real,3}, n_dim_vars)
    n_shot = size(dens_crop, 1)
    n_shot == prod(n_dim_vars) || throw(DimensionMismatch("Loaded $n_shot cropped images but expected $(prod(n_dim_vars)) from dimensions $n_dim_vars."))

    imgs = [copy(@view dens_crop[idx, :, :]) for idx in axes(dens_crop, 1)]
    img_fmt_rev = reshape(imgs, reverse(n_dim_vars)...)
    order_vars = Tuple(reverse(1:length(n_dim_vars)))
    return permutedims(img_fmt_rev, order_vars)
end

function format_dens_runinfo(
    runinfo;
    path_root::AbstractString,
    year_test::Integer,
    wh_corner::Tuple{<:Integer,<:Integer},
    smwh_roi::Tuple{<:Integer,<:Integer},
    len_avg_peak::Integer=10,
    sel_vars=nothing,
)
    val_vars = format_vars(runinfo.vars)
    name_dims = propertynames(val_vars)
    n_dim_vars = Tuple(map(length, val_vars))
    n_variation = prod(n_dim_vars)
    validate_bind_id(runinfo, val_vars)
    dates, runids = if hasproperty(runinfo, :date_runid)
        date_runids = as_vector(runinfo.date_runid)
        length(date_runids) == first(n_dim_vars) || throw(DimensionMismatch("date_runid has length $(length(date_runids)); expected $(first(n_dim_vars)) to match the first variable axis $(first(name_dims))."))
        get_date.(date_runids), get_runid_from_date_runid.(date_runids)
    elseif hasproperty(runinfo, :runids)
        fill(runinfo.date, length(as_vector(runinfo.runids))), as_vector(runinfo.runids)
    elseif hasproperty(runinfo, :runid)
        [runinfo.date], as_vector(runinfo.runid)
    else
        error("runinfo must contain runid, runids, or date_runid.")
    end

    dens = map((date, runid) -> load_dens_run(date, runid; path_root, year_test, wh_corner), dates, runids) |>
           ds -> cat(ds...; dims=1)
    n_shot, h_dens, w_dens = size(dens)
    n_shot == n_variation || throw(DimensionMismatch("Loaded $n_shot shots for $(gen_run_tag(runinfo)), but expected $n_variation from variables $name_dims with dimensions $n_dim_vars."))

    dens_mean = dropdims(mean(dens; dims=1); dims=1)
    xy_peak_px = find_positive_cluster_center(dens_mean, smwh_roi; len_avg=len_avg_peak) |> cent -> round.(Int, cent)
    dens_crop = mapslices(d -> crop_center(d, xy_peak_px, smwh_roi), dens; dims=(2, 3))
    dens_full_fmt = format_image_array(dens_crop, n_dim_vars)
    val_vars_sel, idx_vars = select_val_vars(val_vars, sel_vars)
    dens_full_fmt_sel = dens_full_fmt[idx_vars...]
    n_dim_vars_sel = Tuple(map(length, val_vars_sel))
    runinfo_sel = select_runinfo(runinfo, val_vars_sel, val_vars, idx_vars)

    return (; runinfo=runinfo_sel, val_vars=val_vars_sel, val_vars_full=val_vars, dens_full_fmt=dens_full_fmt_sel, wh_dens=(w_dens, h_dens), xy_peak_px, n_dim_vars=n_dim_vars_sel, name_dims)
end

function read_simulation_grid_units(path_input::AbstractString, filename_xy::AbstractString)
    path_xy = joinpath(path_input, filename_xy)
    isfile(path_xy) || throw(ArgumentError("Simulation grid file not found: $path_xy"))
    file_xy = matopen(path_xy)
    try
        x_mat = read(file_xy, "Y")
        y_mat = read(file_xy, "X")
        unit_x_vals = x_mat[:, 1] |> diff |> unique
        unit_y_vals = y_mat[1, :] |> diff |> unique
        length(unit_x_vals) == 1 || throw(ArgumentError("Expected one x grid interval in $path_xy, got $(unit_x_vals)."))
        length(unit_y_vals) == 1 || throw(ArgumentError("Expected one y grid interval in $path_xy, got $(unit_y_vals)."))
        return (unit_x_vals[1], unit_y_vals[1])
    finally
        close(file_xy)
    end
end

function load_dens_simulation_folder(
    path_input::AbstractString;
    names_istp::AbstractVector{<:AbstractString},
    pattern_filename_data::Regex=Regex(raw"nxy_t=(?<idx_time>\d+).mat"),
    names_density::AbstractVector{<:AbstractString}=["n1", "n2"],
)
    length(names_density) == length(names_istp) ||
        throw(DimensionMismatch("names_density has length $(length(names_density)); expected $(length(names_istp)) to match names_istp."))
    isdir(path_input) || throw(ArgumentError("Simulation input folder not found: $path_input"))

    filenames_data = readdir(path_input) |> fs -> filter(f -> occursin(pattern_filename_data, f), fs)
    isempty(filenames_data) && throw(ArgumentError("No simulation files matching $pattern_filename_data found in $path_input."))

    ids_t = Vector{Int}(undef, length(filenames_data))
    dens_raw = nothing
    for (idx_file, filename) in enumerate(filenames_data)
        ids_t[idx_file] = parse(Int, match(pattern_filename_data, filename)["idx_time"])
        file = matopen(joinpath(path_input, filename))
        try
            dens_first = read(file, first(names_density))
            if dens_raw === nothing
                wh = size(dens_first)
                dens_raw = Array{Float64}(undef, length(filenames_data), length(names_density), wh...)
            end
            dens_raw[idx_file, 1, :, :] .= dens_first
            for idx_istp in 2:length(names_density)
                dens_raw[idx_file, idx_istp, :, :] .= read(file, names_density[idx_istp])
            end
        finally
            close(file)
        end
    end

    perm_t = sortperm(ids_t)
    return ids_t[perm_t], dens_raw[perm_t, :, :, :]
end

function format_dens_simulation_runinfo(
    runinfo;
    path_root::AbstractString,
    smwh_roi::Tuple{<:Integer,<:Integer},
    xy_peak_px::Tuple{<:Integer,<:Integer},
    unit_t::Real,
    unit_in_um::Real,
    filename_xy::AbstractString="XY.mat",
    pattern_filename_data::Regex=Regex(raw"nxy_t=(?<idx_time>\d+).mat"),
    names_density::AbstractVector{<:AbstractString}=["n1", "n2"],
    sel_vars=nothing,
)
    val_vars_base = format_vars(runinfo.vars)
    hasproperty(val_vars_base, :istp) || throw(ArgumentError("simulation runinfo.vars must include istp."))
    names_istp = String.(val_vars_base.istp)
    path_input = joinpath(path_root, runinfo.dir)
    ids_t, dens_raw = load_dens_simulation_folder(path_input; names_istp, pattern_filename_data, names_density)
    unit_x, unit_y = read_simulation_grid_units(path_input, filename_xy)
    px_in_um = (unit_x, unit_y) .* unit_in_um

    for name_required in (:IB, :rep, :istp)
        hasproperty(val_vars_base, name_required) ||
            throw(ArgumentError("simulation runinfo.vars must include $name_required."))
    end
    val_vars_full = (;
        IB=val_vars_base.IB,
        rep=val_vars_base.rep,
        t_hold=ids_t .* unit_t,
        istp=val_vars_base.istp,
    )
    name_dims = propertynames(val_vars_full)
    n_dim_vars = Tuple(map(length, val_vars_full))
    expected = prod(n_dim_vars)
    n_loaded = length(val_vars_full.IB) * length(val_vars_full.rep) * length(ids_t) * length(names_istp)
    expected == n_loaded ||
        throw(DimensionMismatch("Simulation variables describe $expected shots but loader produced $n_loaded from singleton IB/rep, $(length(ids_t)) times, and $(length(names_istp)) istp values."))

    h_dens, w_dens = size(dens_raw, 3), size(dens_raw, 4)
    dens_full_fmt = [
        copy(@view dens_raw[t, i, :, :]) |> transpose |> d -> crop_center(d, xy_peak_px, smwh_roi) |> copy
        for c in eachindex(val_vars_full.IB),
            r in eachindex(val_vars_full.rep),
            t in eachindex(val_vars_full.t_hold),
            i in eachindex(val_vars_full.istp)
    ]

    val_vars_sel, idx_vars = select_val_vars(val_vars_full, sel_vars)
    dens_full_fmt_sel = dens_full_fmt[idx_vars...]
    n_dim_vars_sel = Tuple(map(length, val_vars_sel))
    runinfo_sel = select_runinfo(runinfo, val_vars_sel, val_vars_full, idx_vars)

    return (;
        runinfo=runinfo_sel,
        val_vars=val_vars_sel,
        val_vars_full,
        dens_full_fmt=dens_full_fmt_sel,
        wh_dens=(w_dens, h_dens),
        xy_peak_px,
        n_dim_vars=n_dim_vars_sel,
        name_dims,
        unit_x,
        unit_y,
        unit_t,
        unit_in_um,
        px_in_um,
        ids_t_full=ids_t,
        ids_t=ids_t[idx_vars[3]],
    )
end
