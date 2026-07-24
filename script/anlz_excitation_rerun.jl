using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
using ImageFiltering
GLMakie.activate!()

include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "fitmodels.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))

# commit b93372ff41e21a2b410384b55cac1a4687ecf784
# adjusting on the visualization and processing on profile evolutions
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine"
tag = "NTRC"
title_load = "[07.24].100.Extr"
title_anlz = "[07.24].101.PrflProc.[←100]"

path_load = joinpath(path_root, title_load)
path_output = joinpath(path_root, title_anlz)
isdir(path_output) || mkpath(path_output)

path_load_extr = joinpath(path_load, @sprintf("%s_essn_extr.jld2", tag))
path_load_corr = joinpath(path_load, @sprintf("%s_corr.jld2", tag))

## Recomputed correlation settings. Change these here if the rerun should use
# different analysis choices from the saved extraction metadata.

vis_evol_prfl_modl = (height=100, width_to_time=2.5, ylims=(0, 0.6), colorrange=nothing)
vis_evol_prfl_axial = (height=100, width_to_time=2.5, ylims=nothing, colorrange=nothing)
vis_evol_prfl_radial = (height=100, width_to_time=2.5, ylims=nothing, colorrange=nothing)
vis_evol_prfl_core = vis_evol_prfl_axial
selector_t_hold_prfl_modl = t -> true
selector_pos_prfl_modl = k -> 0.2 < k < 0.4
selector_t_hold_prfl_axial = t -> true
selector_pos_prfl_axial = t -> true
selector_t_hold_prfl_radial = t -> true
selector_pos_prfl_radial = x -> true

##
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
copy_and_include = (name_script) -> begin
    path_script = joinpath(@__DIR__, name_script)
    cp(path_script, joinpath(path_output, basename(path_script)); force=true)
    include(path_script)
end
"load_excitation_extr.jl" |> copy_and_include
"load_excitation_corr.jl" |> copy_and_include
# "anlz_excitation_corr.jl" |> copy_and_include
"anlz_excitation_vslz_corr.jl" |> copy_and_include
# "anlz_excitation_vslz_extr.jl" |> copy_and_include
