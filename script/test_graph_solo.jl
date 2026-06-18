include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
fig = Figure()
gl = GridLayout()
fig[1, 1] = gl

c = 2
t = 10
i = 2
r = 1

axs_live = set_panel_solo_modl!(gl)
info = info_fmt[c, r, t, i]
draw_solo_modl!(axs_live, extr_fmt[c, r, t, i], info)
resize_to_layout!(fig)
fig |> display

# essn_sample = calc_solo_essn_2d(dens_fmt_sample, smwh_peak .+ 1, smwh_peak, 10, 6.5 / 22.);
# info_sample = Dict("istp" => val_vars[3][1], "t_hold" => val_vars[2][3], "repeat" => val_vars[1][1])
# axs_live = set_panel_solo_essn_2d!(gl);
# draw_solo_essn_2d!(axs_live, essn_sample, info_sample);
# fig |> display
