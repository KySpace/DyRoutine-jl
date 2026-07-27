GLMakie.activate!()


path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
title_load = "[07.24].100.Extr"
cache_extr = JLD2.load(joinpath(path_root, "AnlzRoutine", title_load, "NTRC_essn_extr.jld2"))
prfl_axial_evol = JLD2.load(joinpath(path_root, "AnlzRoutine", title_load, "NTRC_corr.jld2"))["prfl_axial_evol"]
meta_extr = cache_extr["meta_extr"]
val_vars = meta_extr.val_vars
px_in_um = meta_extr.px_in_um
coor = (-60:1:60)*px_in_um

num_fit = map(cache_extr["extr_fmt"]) do extr
    p = extr.envelope.params_asymm
    2π * prod(p.size) * p.max
end


istp = 2
ib = findfirst(val_vars.IB .== 5.316)
ids = (ib, 1, :, 2)
rep = 2
sample_axial = prfl_axial_evol[ib, rep, istp]
function twostep_decay(t, p)
    D1, λ1, D2, λ2 = p
    @. D1 * exp(- t / λ1) + D2 * exp(- t / λ2)
end

fig_live = Figure(); 
axs_prfl = Axis(fig_live[1,1]);
axs_num = Axis(fig_live[2,1]);
heatmap!(axs_prfl, val_vars.t_hold, coor, sample_axial'; colormap=Makie.Reverse(:grays), colorrange=(0, 2.5))
for r = 1 : 3
    num_sample = num_fit[ib, r, :, istp] |> vec
    scatter!(axs_num, val_vars.t_hold, num_sample)
end
linkxaxes!([axs_num, axs_prfl])
fig_live |> display