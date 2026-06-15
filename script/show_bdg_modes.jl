using MAT
using Printf
include("../src/graphics.jl")
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
file_gnd = matopen(joinpath(path_root, dir_test, name_grn))
file_uv = matopen(joinpath(path_root, dir_test, name_uv))
φ1 = read(file_gnd, "psi1")
φ2 = read(file_gnd, "psi2")
φ = cat(φ1, φ2; dims=1) |> p -> permutedims(p, (2, 3, 1))
u = read(file_uv, "u") |> real |> u -> reshape(u, (dim_space..., 2, n_mode)) |> u -> permutedims(u, (5, 4, 1, 2, 3))
v = read(file_uv, "v") |> real |> u -> reshape(u, (dim_space..., 2, n_mode)) |> u -> permutedims(u, (5, 4, 1, 2, 3))
norm_uv = (@. u^2 - v^2) |> uv -> sum(uv; dims=(2, 3, 4, 5))
(u, v) = [(u[m,:,:,:,:], v[m,:,:,:,:]) ./ norm_uv[m] for m in 1:n_mode]
dim_z = 5
δρ = [(u.-v)[m, i, :, :, :] .* φ[i] for m in 1:n_mode, i in 1:2] |> w -> drodims(sum(w; dims=dim_z); dims=dim_z)
δφ = [(u.+v)[m, i, :, :, :] ./ φ[i] for m in 1:n_mode, i in 1:2] |> w -> drodims(sum(w; dims=dim_z); dims=dim_z)
close(file_gnd)
close(file_uv)

fig, axs = set_axis!("mode")
m = 1
i = 1
heatmap!(axs, x_vec, y_vec, δρ[m, i, :, :])
