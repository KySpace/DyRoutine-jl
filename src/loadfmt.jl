using HDF5
using Printf
using Statistics

as_vector(x) = x isa AbstractArray ? collect(x) : [x]
format_vars(vars::NamedTuple) = map(as_vector, vars)
format_runids(runids::Integer) = @sprintf("%02d", runids)
format_runids(runids) = runids |> a -> "$(a)" |> s -> replace(s, ":" => "-")

function gen_run_tag(runinfo)
    runids = hasproperty(runinfo, :runids) ? runinfo.runids : runinfo.runid
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
)
    runids = hasproperty(runinfo, :runids) ? as_vector(runinfo.runids) : as_vector(runinfo.runid)
    val = format_vars(runinfo.vars)
    name_dims = propertynames(val)
    n_dim_vars = Tuple(map(length, val))
    n_variation = prod(n_dim_vars)

    dens = map(runid -> load_dens_run(runinfo.date, runid; path_root, year_test, wh_corner), runids) |>
           ds -> cat(ds...; dims=1)
    n_shot, h_dens, w_dens = size(dens)
    n_shot == n_variation || throw(DimensionMismatch("Loaded $n_shot shots for $(gen_run_tag(runinfo)), but expected $n_variation from variables $name_dims with dimensions $n_dim_vars."))

    dens_mean = dropdims(mean(dens; dims=1); dims=1)
    xy_peak_px = find_positive_cluster_center(dens_mean, smwh_roi; len_avg=len_avg_peak) |> cent -> round.(Int, cent)
    dens_crop = mapslices(d -> crop_center(d, xy_peak_px, smwh_roi), dens; dims=(2, 3))
    dens_full_fmt = format_image_array(dens_crop, n_dim_vars)

    return (; runinfo, val, dens_full_fmt, wh_dens=(w_dens, h_dens), xy_peak_px, n_dim_vars, name_dims)
end
