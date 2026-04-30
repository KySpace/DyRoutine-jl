include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
fig = Figure()
gl = GridLayout()
fig[1, 1] = gl
dens_fmt_sample = dens_full_fmt[1, 10, 2, :, :];
essn_sample = calc_solo_essn_2d(dens_fmt_sample, smwh_peak .+ 1, smwh_peak, 10, 6.5 / 22.);
info_sample = Dict("istp" => val[3][2], "t_hold" => val[2][10], "repeat" => val[1][1])
axs = set_panel_solo_essn_2d!(gl);
draw_solo_essn_2d!(axs, essn_sample, info_sample);
fig |> display
