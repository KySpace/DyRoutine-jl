using CairoMakie
include("../multidualmot/dualmotcommons.jl")
path_run17 = raw"\\10.16.18.204\Experimental data\2026-08\0831\run17\liferes.mat";
path_run18 = raw"\\10.16.18.204\Experimental data\2026-08\0831\run18\liferes.mat";
nums_optim_160 = read_cres_fields([path_run17], ("atomnum", )).atomnum |> n -> reshape(n, (9, 3))
nums_optim_162 = read_cres_fields([path_run18], ("atomnum", )).atomnum |> n -> reshape(n, (9, 2))
val_t = [1, 2, 3, 5, 8, 14, 20, 30, 40]
fig = Figure()
ax = Axis(fig[1,1]; xlabel="loading time (s)", ylabel="CMOT N (1e7)",
    yticks=0:2:10)
style = dualmot_curve_style((;loadcfg = :DCS, istp = Symbol("162")))
for optim in 1:2
    nums = (optim == 1 ? nums_optim_160 : nums_optim_162)[:,:] |> n -> mean(n, dims=2) / 1e7 |> vec
    lines!(ax, val_t, nums;
        linestyle= (optim == 1) ? :dash : :solid,
        style.line_options...,
        label="DIS 162 when optimized for $((optim == 1) ? "160" : "162")")
    scatter!(ax, val_t, nums;
        style.marker_options...,
        marker = style.marker,
        markersize=8,
    )
end
axislegend(ax; position=:rb, framevisible=false, labelsize=14)
fig |> display
