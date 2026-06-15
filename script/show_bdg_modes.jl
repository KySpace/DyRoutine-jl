using MAT
using Printf
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
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
int_y = ds -> sum(ds; dims=2) |> ds -> dropdims(ds; dims=2)
int_z = ds -> sum(ds; dims=3) |> ds -> dropdims(ds; dims=3)
file_gnd = matopen(joinpath(path_root, dir_test, name_grn))
file_uv = matopen(joinpath(path_root, dir_test, name_uv))
ψ1 = read(file_gnd, "psi1")
ψ2 = read(file_gnd, "psi2")
ψ = [ψ1, ψ2] |> ds -> map(d -> permutedims(d, (2, 1, 3)), ds)
# fmt_uv_dbg = w -> w |> real |> w -> reshape(w, (dim_space..., 2, n_mode)) |> w -> permutedims(w, (5, 4, 2, 1, 3))
fmt_uv = w -> w |> real |> w -> reshape(w, (dim_space..., 2, n_mode)) |> w -> permutedims(w, (2, 1, 3, 5, 4)) |> w -> eachslice(w; dims=(4, 5))
u = read(file_uv, "u") |> fmt_uv
v = read(file_uv, "v") |> fmt_uv
norm_uv = map((ud, vd) -> sum(@. ud .^ 2 .- vd .^ 2), u, v) |> ns -> sum(ns; dims=2) |> ds -> dropdims(ds; dims=2)
(u, v) = ([u[m, i] ./ norm_uv[m] for m in 1:n_mode, i in 1:2], [v[m, i] ./ norm_uv[m] for m in 1:n_mode, i in 1:2])
δρ_ti = [(u.-v)[m, i] .* ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(int_z, ds)
δφ_ti = [(u.+v)[m, i] ./ ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(int_z, ds)
δρ_si = [(u.-v)[m, i] .* ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(int_y, ds)
δφ_si = [(u.+v)[m, i] ./ ψ[i] for m in 1:n_mode, i in 1:2] |> ds -> map(int_y, ds)
close(file_gnd)
close(file_uv)

fig = Figure()
axs_dens1 = Axis(fig[1, 1]; width=500, height=150)
axs_dens2 = Axis(fig[1, 2]; width=500, height=150)
axs_flow1 = Axis(fig[2, 1]; width=500, height=150)
axs_flow2 = Axis(fig[2, 2]; width=500, height=150)

clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360; thres_alpha=0.05, alpha_base=0.05)
for m = 1:n_mode
    [axs_dens1, axs_dens2, axs_flow1, axs_flow2] |> clear_axes!
    c = maximum(abs, δρ[m, :] |> stack)
    heatmap!(axs_dens1, x_vec, z_vec, δρ_si[m, 1]; colormap=clrmap, colorrange=(-c, c))
    heatmap!(axs_dens2, x_vec, z_vec, δρ_si[m, 2]; colormap=clrmap, colorrange=(-c, c))
    # heatmap!(axs_dens1, x_vec, y_vec, (ψ[1] .* (u.-v)[4, 1]) |> int_z)
    # heatmap!(axs_dens2, x_vec, y_vec, (ψ[2] .* (u.-v)[4, 2]) |> int_z)
    # heatmap!(axs_dens1, x_vec, y_vec, (u.-v)[4, 1] |> int_z)
    # heatmap!(axs_dens2, x_vec, y_vec, (u.-v)[4, 2] |> int_z)
    for ax in [axs_dens1, axs_dens2, axs_flow1, axs_flow2]
        ax.aspect = DataAspect()
        xlims!(ax, (-50, 50))
        ylims!(ax, (-15, 15))
    end
    fig |> resize_to_layout!
    fig |> display
end

# heatmap!(axs, x_vec, y_vec, dropdims(sum(u[m, i, :, :, :] - v[m, i, :, :, :]; dims=3); dims=3))
