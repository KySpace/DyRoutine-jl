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

title = "Anlz.11.BdG-02.[2026.06.16].03"
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"02.[2026.06.12]"
path_show = joinpath(@__DIR__, "show_bdg_modes.jl")
path_output = joinpath(path_root, title)
isdir(path_output) || mkpath(path_output)
n_mode = 15
for a_s_iter in [75, 79, 80, 81, 83, 90]
    global a_s = a_s_iter
    global tag_as = @sprintf("[a_s=%d]", a_s)
    global dim_space = global (n_x, n_y, n_z) = a_s <= 80 ? (256, 256, 64) : (128, 128, 64)
    include(path_show)
end
