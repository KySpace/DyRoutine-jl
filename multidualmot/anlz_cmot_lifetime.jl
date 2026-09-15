# Included by anlz_cmot_lifetime_runner.jl; variable axes are slowest to fastest.
function read_cmot_nums(paths::AbstractVector{<:AbstractString})
    nums = map(paths) do path
        matopen(path) do file
            cres = read(file, "liferes")["cres"]
            # MatlabStructArray length counts fields; atomnum's shape counts shots.
            entries = if cres isa MAT.MatlabStructArray
                [cres[i] for i in 1:length(cres["atomnum"])]
            elseif cres isa AbstractDict
                [cres]
            else
                vec(cres)
            end
            map(entries) do entry
                value = entry["atomnum"]
                value = value isa AbstractArray ? only(value) : value
                value isa Real && isfinite(value) ||
                    throw(ArgumentError("$path: atomnum must be a finite numeric scalar"))
                Float64(value)
            end
        end
    end
    vcat(nums...)
end

num_data = read_cmot_nums(runinfo.files)
len_data = length(num_data)
n_dim_vars = length.((val_tbiasmot, val_t_hold, val_loadcfg, val_istp))
n_variation = prod(n_dim_vars)
len_data > 0 && rem(len_data, n_variation) == 0 ||
    throw(DimensionMismatch("$tag_head: $len_data atom numbers must be a positive integer multiple of $n_variation variations"))
n_rep = div(len_data, n_variation)
val_rep = 1:n_rep
vars = merge(runinfo.vars, (; rep=val_rep))
name, val_vars = keys(vars), values(vars)
n_dims_fmt = (n_rep, n_dim_vars...)
num_fmt = reshape(num_data, reverse(n_dims_fmt)) |> data -> permutedims(data, (5, 4, 3, 2, 1))
# Statistics retain (β_MOT, t_hold, loadcfg, istp); std is the sample standard deviation.
num_stat = dropdims(mean(num_fmt; dims=1); dims=1)
std_num_stat = dropdims(std(num_fmt; dims=1); dims=1)
n_rep == 1 && @warn "$tag_head: one repetition; sample standard deviations are undefined (NaN)"
println("$tag_head: $len_data samples / $n_variation variations = $n_rep repetitions")

# Replot from here using the existing statistics.
figs_cmot = Dict{Tuple{Float64, Symbol}, Figure}()
for (idx_bias, bias) in enumerate(val_tbiasmot), scale in (:lin, :log)
    balance = iszero(bias) ? "t" : "n"
    fig = Figure(size=size_figure)
    ax = Axis(fig[1, 1]; xlabel="t_hold (ms)", ylabel="CMOT atom number",
        title="$tag_head · β_MOT = $bias · $balance-balanced",
        xscale=scale == :log ? log10 : identity,
        yscale=scale == :log ? log10 : identity)
    for loadcfg in (:DDM, :DIS), (idx_istp, istp) in enumerate(val_istp)
        idx_loadcfg = findfirst(==(loadcfg), val_loadcfg)
        nums = @view num_stat[idx_bias, :, idx_loadcfg, idx_istp]
        stds = @view std_num_stat[idx_bias, :, idx_loadcfg, idx_istp]
        hue = hue_istp[istp]
        stroke = RGB(Oklch(lightness_stroke, chroma_stroke, hue))
        face = RGB(Oklch(lightness_face, chroma_face, hue))
        # Undefined log points become gaps; never modify the underlying statistics.
        mask = scale == :log ? nums .> 0 : trues(length(nums))
        if !all(mask)
            @warn "$tag_head: omitting nonpositive log points" bias loadcfg istp count=count(!, mask)
        end
        nums_plot = ifelse.(mask, nums, NaN)
        lines!(ax, val_t_hold, nums_plot; color=stroke, linewidth=2)
        scatter!(ax, val_t_hold, nums_plot; color=face, strokecolor=stroke,
            strokewidth=1.5, markersize=11, marker=marker_loadcfg[loadcfg],
            label="$(string(istp)) $loadcfg")
        mask_error = mask .& isfinite.(stds)
        if scale == :log
            # A mean - std <= 0 has no logarithm: omit that error bar and report it.
            mask_error .&= nums .- stds .> 0
            count_crossing = count(mask .& isfinite.(stds) .& (nums .- stds .<= 0))
            count_crossing > 0 && @warn "$tag_head: log error bars crossing zero omitted" bias loadcfg istp count_crossing
        end
        errorbars!(ax, val_t_hold[mask_error], nums[mask_error], stds[mask_error];
            color=stroke, whiskerwidth=7, linewidth=1.2)
    end
    axislegend(ax; position=:rt)
    path_svg = joinpath(path_output, "[CMOT.lifetime].[$scale].[$balance-balanced].svg")
    save(path_svg, fig)
    figs_cmot[(bias, scale)] = fig
end
