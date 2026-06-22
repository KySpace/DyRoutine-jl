## extraction figures from solo or stack

draw_solo_modl_kwargs = @isdefined(draw_solo_modl_kwargs) ? draw_solo_modl_kwargs : NamedTuple()
fig_full, axs_solo, axs_stacked = set_axis_full(n_dim_vars_per_IB, set_panel_solo_modl!)
println("Full axes ready: dimensions $(n_dim_vars_per_IB)")
for c in 1:n_dim_vars[1]
    tag_IB = tag_IBs[c]
    local t_stage = log_step("Now plotting full modulation table for $tag_IB.")
    for t in 1:n_dim_vars[3], i in 1:n_dim_vars[4]
        for r in 1:n_dim_vars[2]
            info = info_fmt[c, r, t, i]
            print("\r\033[2Kplotting for runid $(info["runid"]), rep $r, $(info["t_hold"]) ms, $(info["istp"])")
            draw_solo_modl!(axs_solo[r, t, i], extr_fmt[c, r, t, i], info; draw_solo_modl_kwargs...)
        end
        info = info_fmt[c, 1, t, i] |> d -> merge(d, Dict("repeat" => "stacked"))
        print("\r\033[2Kplotting for stacked runid $(info["runid"]), $(info["t_hold"]) ms, $(info["istp"])")
        draw_solo_modl!(axs_stacked[t, i], extr_stacked_over_rep[c, t, i], info; draw_solo_modl_kwargs...)
    end
    println("")
    println("Full modulation table drawn.")
    resize_to_layout!(fig_full)
    fig_full |> f -> save(joinpath(path_output, @sprintf("solo_table_[%s].pdf", tag_IB)), f; backend=CairoMakie)
    log_done("Full modulation plot saved for $tag_IB.", t_stage)
end
