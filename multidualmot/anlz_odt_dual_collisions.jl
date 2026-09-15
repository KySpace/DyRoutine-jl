using MAT
using CairoMakie: Axis, Figure, GridLayout, Label, lines!, save, scatter!, text!
using LsqFit: coef, curve_fit
using Pipe: @pipe
using Printf: @sprintf
path_root = raw"\\10.16.18.204\Experimental data\2026-09\0908"
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\IsotopePairs\2026-09\0908"
val_istp = [:161, :163]
val_type = [:MX, :Sw]
val_626_163 = [:detuned, :on_resn]
val_p_odt = [0.0, 3.0, 5.4]
val_t_hold = [25, 50, 100, 200, 500, 1500, 4000]
n_dims_fmt = @pipe (val_626_163, val_p_odt, val_t_hold, val_type) |> map(length, _)

function read_nums(dirs_test::AbstractVector{<:AbstractString})
    nums_raw = map(dirs_test) do dir_test
        file_liferes = matopen(joinpath(path_root, dir_test, "liferes.mat"))
        cres = read(file_liferes, "liferes")["cres"]
        [cres[i]["atomnum"] for i in 1:42] |> vec
    end
    nums_raw = vcat(nums_raw...)
    @pipe nums_raw |> reshape(_, Tuple(reverse(n_dims_fmt))) |> permutedims(_, Tuple(reverse(1:length(n_dims_fmt))))
end

dirs_test = [raw"run58", raw"run59"]
nums = read_nums(dirs_test)

decays = [nums[idx_626_163, p, :, ty] for idx_626_163 in axes(nums, 1), p in axes(nums, 2), ty in axes(nums, 4)]

function exp_decay(t::AbstractVector, params::AbstractVector)
    num_initial, tau_decay, background = params
    @. num_initial * exp(-t / tau_decay) + background
end

mask_t_fit = @pipe [0, 1, 1, 1, 1, 1, 1] |> Bool.(_)

function fit_exp_decay(t::AbstractVector{<:Real}, decay::AbstractVector{<:Real}, mask_t::AbstractVector{<:Bool}=trues(size(t)))
    t_fit = Float64.(t) ./ 1e3
    decay_fit = Float64.(decay)
    decay_fit[begin] > 0 || throw(ArgumentError("the first decay value must be positive"))
    all(decay_fit .>= 0) || throw(ArgumentError("decay values must be nonnegative"))

    background_initial = 0 #minimum(decay_fit)
    num_initial_initial = max(decay_fit[begin] - background_initial, eps(Float64))
    tau_initial = decay_fit[end] > background_initial && t_fit[end] > t_fit[begin] ?
        max(eps(Float64), (t_fit[end] - t_fit[begin]) / log(num_initial_initial / (decay_fit[end] - background_initial))) : 1.0
    params_initial = [num_initial_initial, tau_initial, background_initial]
    fit = curve_fit(exp_decay,
                    t_fit[mask_t], decay_fit[mask_t], params_initial;
                    lower=[0.0, eps(Float64), 0],
                    # upper=[maximum(decay_fit), Inf, 0]
                    )
    coef(fit)
end

decay_fits = map(decays) do decay
    fit_exp_decay(val_t_hold, decay)
end
num_initial = getindex.(decay_fits, 1)
tau_decay = getindex.(decay_fits, 2)

clr_type = Dict(:MX => :darkviolet, :Sw => :darkblue)
idx_626_163_plot = [findfirst(==(value), val_626_163) for value in (:on_resn, :detuned)]
t_hold_sec = val_t_hold ./ 1e3
t_fit_plot = range(first(t_hold_sec), last(t_hold_sec); length=200)

fig_decay = Figure(size=(1200, 900))
for (idx_col, idx_626_163) in enumerate(idx_626_163_plot)
    Label(fig_decay[1, idx_col], string(val_626_163[idx_626_163]); tellwidth=false)
    for (idx_row, idx_p) in enumerate(axes(nums, 2))
        gl_panel = GridLayout(fig_decay[idx_row + 1, idx_col])
        for (idx_inner_col, idx_type) in enumerate(axes(nums, 4))
            ax = Axis(
                gl_panel[1, idx_inner_col];
                title=string(val_type[idx_type]),
                xlabel=idx_row == length(axes(nums, 2)) ? "t_hold (sec)" : "",
                ylabel=idx_inner_col == 1 ? "atom number (1e6)" : "",
                limits=(nothing, (-0.2, 5.2)),
            )
            decay = nums[idx_626_163, idx_p, :, idx_type]
            params = decay_fits[idx_626_163, idx_p, idx_type]
            color = clr_type[val_type[idx_type]]
            scatter!(ax, t_hold_sec, decay ./ 1e6; color, markersize=8)
            lines!(ax, t_fit_plot, exp_decay(t_fit_plot, params) ./ 1e6; color, linewidth=2)
            ratio_num = nums[idx_626_163, idx_p, 4, idx_type] / nums[idx_626_163, idx_p, 2, idx_type]
            text!(
                ax,
                0.97,
                0.97;
                text="τ_decay = $(@sprintf("%.2f", params[2])) sec\n" *
                     "num_initial = $(@sprintf("%.3g", params[1]))\n" *
                     "N(200 ms) / N(50 ms) = $(@sprintf("%.2f", ratio_num))",
                space=:relative,
                align=(:right, :top),
                fontsize=11,
            )
        end
    end
end

dir_output = joinpath(path_root, last(dirs_test))
path_svg = joinpath(dir_output, "anlz_odt_dual_collissions.svg")
save(path_svg, fig_decay)
cp(abspath(@__FILE__), joinpath(dir_output, basename(@__FILE__)); force=true)

## commonly used comparisons
# 200 ms / 50 ms
ratio_200_50 = nums[:, :, 4, :] ./ nums[:, :, 2, :]
ratio_decay_simu = @pipe tau_decay |> exp.(-0.15./_) |> round.(_; digits=2)
@pipe ratio_decay_simu |> vcat(_[2,3,:] ./ _[2,1,:], _[1,3,:] ./ _[1,1,:])'
rate_decay = @pipe tau_decay |> 1 ./ _

display_table = q -> @pipe [q[a,:,:] for a in 1:2] |> hcat(_...) |> round.(_; digits=2)

tau_decay_compose = rate_decay |> (rate -> [
        if a == 1 || p == 1
            NaN
        else
            rate[1, p, ty] + rate[a, 1, ty] - rate[1, 1, ty]
        end
        for a in axes(rate, 1), p in axes(rate, 2), ty in axes(rate, 3)]) |>
            a -> 1 ./ a
