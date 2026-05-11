## Run this after the data has been processed
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))

extr_fmt = [
    essn_2d_fmt[r, t, i] |> e -> calc_solo_extr(e, fit_prfl_modl_over_rep_t_1d[i])
    for r in axes(essn_2d_fmt, 1), t in axes(essn_2d_fmt, 2), i in axes(essn_2d_fmt, 3)
]

GLMakie.activate!()
fig_live = Figure()
gl = GridLayout()
fig_live[1, 1] = gl
axs_live = set_panel_solo_modl!(gl);
ids_demo = (3, 10, 2)
info_demo = info_fmt[ids_demo...]
extr_demo = extr_fmt[ids_demo...]
draw_solo_modl!(axs_live, extr_demo, info_demo)
fig_live |> display
