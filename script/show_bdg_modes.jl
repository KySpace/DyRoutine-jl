name_grn = @sprintf("psi_as12=%s.mat", a_s)
name_uv = @sprintf("spectrum_as12=%s.mat", a_s)
smlx, smly, smlz = (100, 30, 40)
x_vec = range(-smlx, smlx, n_x);
y_vec = range(-smly, smly, n_y);
z_vec = range(-smlz, smlz, n_z);
ids_cut = (pos, dims) -> ntuple(i -> i == pos ? cld(dims[pos], 2) : Colon(), length(dims))
int_y = ds -> sum(ds; dims=2) |> ds -> dropdims(ds; dims=2)
int_z = ds -> sum(ds; dims=3) |> ds -> dropdims(ds; dims=3)
cut_y = ds -> ds[ids_cut(2, size(ds))...]
cut_z = ds -> ds[ids_cut(3, size(ds))...]
file_gnd = matopen(joinpath(path_root, dir_test, name_grn))
file_uv = matopen(joinpath(path_root, dir_test, name_uv))
ψ1 = read(file_gnd, "psi1")
ψ2 = read(file_gnd, "psi2")
ψ_raw = [ψ1, ψ2] |> ds -> map(d -> permutedims(d, (2, 1, 3)), ds)
norm_ψ = ψ_raw |> cs -> map(c -> sum(abs2, c), cs) |> sum |> sqrt
ψ = ψ_raw ./ norm_ψ
# fmt_uv_dbg = w -> w |> real |> w -> reshape(w, (dim_space..., 2, n_mode)) |> w -> permutedims(w, (5, 4, 2, 1, 3))
fmt_uv = w -> w |> real |> w -> reshape(w, (dim_space..., 2, n_mode)) |> w -> permutedims(w, (2, 1, 3, 5, 4)) |> w -> eachslice(w; dims=(4, 5))
u = read(file_uv, "u") |> fmt_uv
v = read(file_uv, "v") |> fmt_uv
ω = read(file_uv, "omega_nv")
norm_uv = map((ud, vd) -> sum(ud .^ 2 .- vd .^ 2), u, v) |>
          ns -> sum(ns; dims=2) |>
                ns -> dropdims(ns; dims=2) |>
                      ns -> sqrt.(ns)
(u, v) = ([u[m, i] ./ norm_uv[m] for m in 1:n_mode, i in 1:2], [v[m, i] ./ norm_uv[m] for m in 1:n_mode, i in 1:2])
mask_relavent_4d = map(c -> c .> maximum(c) / 10, ψ) |> stack
δρ_3d = [(u.-v)[m, i] .* ψ[i] for m in 1:n_mode, i in 1:2]
δφ_3d = [(u.+v)[m, i] ./ ψ[i] for m in 1:n_mode, i in 1:2]
δρ_ti = δρ_3d |> ds -> map(int_z, ds)
δφ_ti = δφ_3d |> ds -> map(cut_z, ds)
δρ_si = δρ_3d |> ds -> map(int_y, ds)
δφ_si = δφ_3d |> ds -> map(cut_y, ds)
max_φ = [[δφ_3d[m, i] .* mask_relavent_4d for i in 1:2] |> stack |> f -> abs.(f) |> maximum for m in 1:n_mode]
close(file_gnd)
close(file_uv)

function set_panel_mode!(gl::GridLayout)
    gl |> clean_gridlayout!
    dict_axs = Dict{String,Vector{Axis}}()
    dict_axs["δρ_si"] = [Axis(gl[1, 1], aspect=DataAspect(), width=500, height=150), Axis(gl[1, 2], aspect=DataAspect(), width=500, height=150)]
    dict_axs["δφ_si"] = [Axis(gl[2, 1], aspect=DataAspect(), width=500, height=150), Axis(gl[2, 2], aspect=DataAspect(), width=500, height=150)]
    dict_axs["δρ_ti"] = [Axis(gl[3, 1], aspect=DataAspect(), width=500, height=150), Axis(gl[3, 2], aspect=DataAspect(), width=500, height=150)]
    dict_axs["δφ_ti"] = [Axis(gl[4, 1], aspect=DataAspect(), width=500, height=150), Axis(gl[4, 2], aspect=DataAspect(), width=500, height=150)]
    dict_axs
end

println("  [$tag] Drawing modes $tag_as")
fig_modes = Figure()

function to_phase_clr_dens(φ, ρ, hue1, hue2; max_ρ=maximum(ρ), max_φ=maximum(φ), thres_alpha=0.1, alpha_base=0.0, l=0.4, l_max=0.8, h=0.24)
    size(φ) == size(ρ) || throw(ArgumentError("φ and ρ must have the same size"))
    ρ_n = (d -> clamp.(d, 0, max_ρ) / max_ρ)(ρ)
    φ_n = (d -> clamp.(d, -max_φ, max_φ) / max_φ)(φ)
    alpha = n -> n > thres_alpha ? 1.0 : (n / thres_alpha * (1 - alpha_base) + alpha_base)
    shader = (f, n) -> Oklch(l_max - (l_max - l) * abs(f), h * abs(f), f > 0 ? hue1 : hue2) |> c -> RGBAf(c, alpha(n))
    return [shader(φ_n[x, y], ρ_n[x, y]) for x in 1:size(ρ, 1), y in 1:size(ρ, 2)]
end

clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360; thres_alpha=0.05, alpha_base=0.05)

for m = 1:n_mode
    print("\r      [$tag_as] Drawing mode $m")
    flush(stdout)
    gl_modes = GridLayout(fig_modes[m, 1])
    axs_modes = set_panel_mode!(gl_modes)
    axs_modes |> clear_axes!

    c_t = maximum(abs, δρ_ti[m, :] |> stack)
    c_s = maximum(abs, δρ_si[m, :] |> stack)

    for i = 1:2
        clr_δφ_si = to_phase_clr_dens(δφ_si[m, i], δρ_si[m, i], 0.57 * 360, 0.96 * 360; max_φ=max_φ[m])
        clr_δφ_ti = to_phase_clr_dens(δφ_ti[m, i], δρ_ti[m, i], 0.57 * 360, 0.96 * 360; max_φ=max_φ[m])
        heatmap!(axs_modes["δρ_si"][i], x_vec, y_vec, δρ_si[m, i]; colormap=clrmap, colorrange=(-c_s, c_s), rasterize=true)
        heatmap!(axs_modes["δρ_ti"][i], x_vec, z_vec, δρ_ti[m, i]; colormap=clrmap, colorrange=(-c_t, c_t), rasterize=true)
        heatmap!(axs_modes["δφ_si"][i], x_vec, y_vec, clr_δφ_si; rasterize=true)
        heatmap!(axs_modes["δφ_ti"][i], x_vec, z_vec, clr_δφ_ti; rasterize=true)
    end
    for ax in values(axs_modes), i = 1:2
        ax[i] |> a -> hidedecorations!(a, ticks=true, ticklabels=true, grid=false)
        ax[i].xgridvisible = true
        ax[i].ygridvisible = true
        ax[i].xminorgridvisible = true
        ax[i].yminorgridvisible = true
        xlims!(ax[i], (-50, 50))
        ylims!(ax[i], (-15, 15))
    end
    fig_modes |> resize_to_layout!
end
fig_modes |> resize_to_layout!
fig_modes |> f -> save(joinpath(path_output, "Modes.$tag_as.png"), f; backend=GLMakie)
println("")


println("  [$tag] Animating modes $tag_as")

fig_gif = Figure()
label_gif = Label(fig_gif[0, 0], "$tag_as"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
ax_ti = Axis(fig_gif[1, 1]; width=500, height=150, aspect=DataAspect())
ax_si = Axis(fig_gif[2, 1]; width=500, height=200, aspect=DataAspect())
ax_prfl = Axis(fig_gif[3, 1]; width=500, height=150)
rowsize!(fig_gif.layout, 0, 20)

for m in 1:15
    print("\r      [$tag_as] Animating mode $m")
    flush(stdout)
    clr_theme = [RGBAf(Oklch(0.52, 0.14, hue_theme_istp[i]), 0.75) for i in ["162" "164"]]
    dens_t = t -> [@. abs2(ψ[i]) + sin(t * 2 * π) * ψ[i] * (u-v)[m, i] for i in 1:2]
    dens_xz_base = map(c -> abs2.(c) |> int_y, ψ)
    dens_xy_base = map(c -> abs2.(c) |> int_z, ψ)
    dens_x_base = map(c -> abs2.(c) |> int_z |> int_y, ψ)
    uv_xz = map(int_y, [@. ψ[i] * (u-v)[m, i] for i in 1:2])
    uv_xy = map(int_z, [@. ψ[i] * (u-v)[m, i] for i in 1:2])
    uv_x = [@. ψ[i] * (u-v)[m, i] for i in 1:2] |> uv -> map(w -> w |> int_z |> int_y, uv)
    mask_relavent = map(c -> c .> maximum(c) / 10, ψ) |> stack
    max_prfl = dens_x_base |> stack |> n -> filter(!isnan, n) |> n -> maximum(n; init=0) |> n -> n * 1.5
    uv_scaler = (abs.(stack((u.-v)[m, :])) ./ stack(ψ)) |>
                ra -> ra[mask_relavent] |>
                      maximum |> max -> 1.0 / max
    dens_xzt = t -> [@. dens_xz_base[i] + uv_scaler * sin(t * 2 * π) * uv_xz[i] for i in 1:2]
    dens_xyt = t -> [@. dens_xy_base[i] + uv_scaler * sin(t * 2 * π) * uv_xy[i] for i in 1:2]
    dens_xt = t -> [@. dens_x_base[i] + uv_scaler * sin(t * 2 * π) * uv_x[i] for i in 1:2]
    # dens_t = t -> [@. sin(t * 2 * π) * ψ[i] * (u-v)[m, i] for i in 1:2]
    function local_draw(t)
        [ax_ti, ax_si, ax_prfl] |> clear_axes!
        label_gif.text = @sprintf("%d a₀ | Mode %d | %.2f Hz | scale = %.03f | t=%.2f", a_s, m, real(ω[m]), uv_scaler, t)
        dens_xz_frame = dens_xzt(t)
        dens_xy_frame = dens_xyt(t)
        dens_x_frame = dens_xt(t)
        clr_misc_xz = to_miscibility_clr(dens_xz_frame[1], dens_xz_frame[2], hue_theme_istp["162"], hue_theme_istp["164"]; to_norm_each=false, max=0.002)
        clr_misc_xy = to_miscibility_clr(dens_xy_frame[1], dens_xy_frame[2], hue_theme_istp["162"], hue_theme_istp["164"]; to_norm_each=false, max=0.002)
        heatmap!(ax_ti, x_vec, y_vec, clr_misc_xy; rasterize=true)
        heatmap!(ax_si, x_vec, z_vec, clr_misc_xz; rasterize=true)
        lines!(ax_prfl, x_vec, dens_x_frame[1]; color=clr_theme[1])
        lines!(ax_prfl, x_vec, dens_x_frame[2]; color=clr_theme[2])
        ax_ti.xticks = LinearTicks(10)
        ax_si.xticks = LinearTicks(10)
        ax_prfl.xticks = LinearTicks(10)
        ax_si.xticklabelsvisible = false
        ax_ti.xticklabelsvisible = false
        xlims!(ax_ti, (-50, 50))
        ylims!(ax_ti, (-15, 15))
        xlims!(ax_si, (-50, 50))
        ylims!(ax_si, (-20, 20))
        xlims!(ax_prfl, (-50, 50))
        ylims!(ax_prfl, (-0.005, max_prfl))
    end
    local_draw(0)
    linkxaxes!(ax_ti, ax_si, ax_prfl)
    fig_gif |> resize_to_layout!
    record(local_draw, fig_gif, joinpath(path_output, "$tag_as.mode-$m.gif"), range(0, 1, 72); framerate=24)
end
println("")
println("  [$tag] Finished $tag_as")
