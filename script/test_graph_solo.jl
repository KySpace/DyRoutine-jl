include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
fig = Figure()
gl = GridLayout()
fig[1, 1] = gl
dens_fmt_sample = dens_full_fmt[1, 3, 1, :, :];
essn_sample = calc_solo_essn_2d(dens_fmt_sample, smwh_peak .+ 1, smwh_peak, 10, 6.5 / 22.);
info_sample = Dict("istp" => val[3][1], "t_hold" => val[2][3], "repeat" => val[1][1])
axs_live = set_panel_solo_essn_2d!(gl);
draw_solo_essn_2d!(axs_live, essn_sample, info_sample);
fig |> display

essn_2d_fmt = dens_full_fmt |> ds -> mapslices(d -> calc_solo_essn_2d(d, smwh_peak .+ 1, smwh_peak, 10, 6.5 / 22.), ds; dims=(4, 5))  |> e -> dropdims(e; dims=(4,5));
info_fmt = [Dict("istp" => val[3][i], "t_hold" => val[2][t], "repeat" => val[1][r]) for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]]
fig_full, axs_solo = set_axis_full(n_dim_vars)
for r in 1:n_dim_vars[1], t in 1:n_dim_vars[2], i in 1:n_dim_vars[3]
    info = info_fmt[r, t, i]
    print("\rplotting for rep $i, $(info["t_hold"]) ms, $(info["istp"])")
    draw_solo_essn_2d!(axs_solo[r, t, i], essn_2d_fmt[r, t, i], info)
    draw_solo_essn_2d!(axs_live, essn_2d_fmt[r, t, i], info)
end
resize_to_layout!(fig_full)

fig_full |> f -> save(joinpath(@__DIR__, "full_essn.pdf"), f; backend=CairoMakie)
