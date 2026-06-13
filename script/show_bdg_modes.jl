using MAT
using Printf
tag = "BdG"
log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))

title = "Anlz.09.Simu-02.[2026.06.12].01"
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"02.[2026.06.12]"
path_output = joinpath(path_root, title)
str_as = 75
name_ground = @sprintf("psi_as12=%s.mat", str_as)
name_uv = @sprintf("spectrum_as12=%s.mat", str_as)
x_vec = range(-100, 100, 256);
y_vec = range(-30, 30, 256);
z_vec = range(-40, 40, 64);
dim_space = (256, 256, 64)
# x_vec = range(-100, 100, 128);
# y_vec = range(-30, 30, 128);
# z_vec = range(-40, 40, 64);
#
file_gnd = matopen(joinpath(path_root, dir_test, name_ground))
file_uv = matopen(joinpath(path_root, dir_test, name_uv))
φ1 = read(file_gnd, "psi1")
φ2 = read(file_gnd, "psi2")
φ = cat(φ1, φ2; dims=1) |> p -> permutedims(p, (2, 3, 1))
u1 = read(file_uv, "u") |> real |> u -> reshape(u, (dim_space..., 2, 15)) |> u -> permutedims!(u, (3, 4, 5, 2, 1))
v1 = read(file_uv, "v") |> real |> u -> reshape(u, (dim_space..., 2, 15)) |> u -> permutedims!(u, (3, 4, 5, 2, 1))
close(file_gnd)
close(file_uv)
