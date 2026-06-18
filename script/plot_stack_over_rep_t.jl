include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
fig_live = Figure()
gl_live = GridLayout()
fig_live[1,1] = gl_live
axs_live = set_panel_solo_modl!(gl_live)


c = 3
istp = 1
essn = essn_stacked_over_rep_t[c,istp]
finess = 10
x, y = essn.smwh_core |> s -> map(u -> (-u:(1/finess):u), s)
x_modl, y_modl = (x, y) .* essn.step_modl .* finess


fit_result = fit_prfl_modl_over_rep_t_1d[c, istp]
draw_solo_essn_2d!(
    axs_live, essn,
    Dict(
        "istp" => val_vars.istp[istp],
        "t_hold" => 0,
        "repeat" => 0,
        "runid" => get_bind_runid(runinfo, c),
        "IB" => val_vars.IB[c],
    ))

band!(axs_live["upright"], y_modl, 0, fit_result.tail(y_modl), color=(:indianred, 0.7)) |> b -> translate!(b, 0, 0, -1)
band!(axs_live["upright"], y_modl, fit_result.tail(y_modl), fit_result.model(y_modl, fit_result.params), color=(:darkseagreen1, 0.5)) |> b -> translate!(b, 0, 0, -1)
xlims!(axs_live["upright"], 0, 1.2)

axs_live["sideway"] |> hidedecorations!
fig_live |> resize_to_layout!
fig_live |> display
fig_live |> f -> save(joinpath(@__DIR__, "..", "tests", "stack_over_rep_t_162.svg"), f; backend=CairoMakie)
