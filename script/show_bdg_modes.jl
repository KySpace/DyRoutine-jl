using MAT
using Printf
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "visduet.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))
tag = "BdG"
log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))

title = "Anlz.09.Simu-02.[2026.06.12].01"
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"02.[2026.06.12]"
path_output = joinpath(path_root, title)
str_as = 75
name_grn = @sprintf("psi_as12=%s.mat", str_as)
name_uv = @sprintf("spectrum_as12=%s.mat", str_as)
dim_space = (n_x, n_y, n_z) = (256, 256, 64)
n_mode = 15
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
norm_uv = map((ud, vd) -> sum(ud .^ 2 .- vd .^ 2), u, v) |>
          ns -> sum(ns; dims=2) |>
                ns -> dropdims(ns; dims=2) |>
                      ns -> sqrt.(ns)
(u, v) = ([u[m, i] ./ norm_uv[m] for m in 1:n_mode, i in 1:2], [v[m, i] ./ norm_uv[m] for m in 1:n_mode, i in 1:2])
δρ_ti = [(u.-v)[m, i] .* ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(int_z, ds)
δφ_ti = [(u.+v)[m, i] ./ ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(cut_z, ds)
δρ_si = [(u.-v)[m, i] .* ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(int_y, ds)
δφ_si = [(u.+v)[m, i] ./ ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(cut_y, ds)
close(file_gnd)
close(file_uv)

function set_panel_mode!(gl::GridLayout)
    gl |> clean_gridlayout!
    dict_axs = Dict{String,Vector{Axis}}()
    dict_axs["δρ_si"] = [Axis(gl[1, 1], aspect=DataAspect()), Axis(gl[1, 2], aspect=DataAspect())]
    dict_axs["δφ_si"] = [Axis(gl[2, 1], aspect=DataAspect()), Axis(gl[2, 2], aspect=DataAspect())]
    dict_axs["δρ_ti"] = [Axis(gl[3, 1], aspect=DataAspect()), Axis(gl[3, 2], aspect=DataAspect())]
    dict_axs["δφ_ti"] = [Axis(gl[4, 1], aspect=DataAspect()), Axis(gl[4, 2], aspect=DataAspect())]
    dict_axs
end


fig_modes = Figure()
gl_modes = GridLayout(fig_modes[1, 1])
axs_modes = set_panel_mode!(gl_modes)

clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360; thres_alpha=0.05, alpha_base=0.05)
for m = 1:n_mode
    axs_modes |> clear_axes!
    c_t = maximum(abs, δρ_ti[m, :] |> stack)
    c_s = maximum(abs, δρ_si[m, :] |> stack)
    for i = 1:2
        heatmap!(axs_modes["δρ_si"][i], x_vec, y_vec, δρ_si[m, 1]; colormap=clrmap, colorrange=(-c_s, c_s))
        heatmap!(axs_modes["δρ_ti"][i], x_vec, z_vec, δρ_ti[m, 1]; colormap=clrmap, colorrange=(-c_t, c_t))
        heatmap!(axs_modes["δφ_si"][i], x_vec, y_vec, δφ_si[m, 1]; colormap=clrmap, colorrange=(-π, π) .* 0.05)
        heatmap!(axs_modes["δφ_ti"][i], x_vec, z_vec, δφ_ti[m, 1]; colormap=clrmap, colorrange=(-π, π) .* 0.05)
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
    fig_modes |> display
end

fig_gif, ax_gif = set_axis!("Density Evolution")

m = 6

n_t = t -> [@. abs2(ψ[i]) + cos(t * 2 * π) * ψ[i] * (u-v)[m, i] for i in 1:2]
n_y = map(c -> abs2.(c) |> int_y, ψ)
uv_y = map(int_y, [@. ψ[i] * (u-v)[m, i] for i in 1:2])
n_yt = t -> [@. n_y[i] + cos(t * 2 * π) * uv_y[i] for i in 1:2]
# n_t = t -> [@. cos(t * 2 * π) * ψ[i] * (u-v)[m, i] for i in 1:2]
for t = range(0, 6, 120)
    [ax_gif] |> clear_axes!
    ax_gif.title = @sprintf("Mode %d, t=%.2f", m, t)
    n = n_yt(t)
    clr_misc = to_miscibility_clr(n[1], n[2], hue_theme_istp["162"], hue_theme_istp["164"]; to_norm_each=false, max=0.01)
    heatmap!(ax_gif, x_vec, z_vec, clr_misc; rasterize=true)
    xlims!(ax_gif, (-50, 50))
    ylims!(ax_gif, (-15, 15))
    ax_gif.aspect = DataAspect()
    fig_gif |> resize_to_layout!
    fig_gif |> display
    sleep(0.01)
end

# heatmap!(axs, x_vec, y_vec, dropdims(sum(u[m, i, :, :, :] - v[m, i, :, :, :]; dims=3); dims=3))
