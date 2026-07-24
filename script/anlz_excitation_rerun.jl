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
# There is going to be a calibrated version of density profiles
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine"
tag = "CFNM"
title_load = "[07.24].100.Extr"
title_anlz = "[07.24].102.PrflProc.[←100]"

path_load = joinpath(path_root, title_load)
path_output = joinpath(path_root, title_anlz)
isdir(path_output) || mkpath(path_output)

path_load_extr = joinpath(path_load, @sprintf("%s_essn_extr.jld2", tag))
path_load_corr = joinpath(path_load, @sprintf("%s_corr.jld2", tag))

## Recomputed correlation settings. Change these here if the rerun should use
# different analysis choices from the saved extraction metadata.

selector_t_pca_dens = t -> 30 .< t .< 80
selector_t_pca_modl = t -> 30 .< t .< 200
n_pca_modes_prfl_modl = 8
freq_query_pca_modl = 1:1:100
selector_t_spectrum = (;
    number=t -> 30 .< t .< 100,
    sp_weight=t -> 30 .< t .< 100,
    sp_height=t -> 30 .< t .< 100,
    sp_width=t -> 30 .< t .< 100,
    sp_wavenum=t -> 30 .< t .< 100,
    nvlp=t -> 30 .< t .< 200,
)
filter_core_pca_sigma = 1.5
filter_core_pca = im -> imfilter(im, Kernel.gaussian(filter_core_pca_sigma))
query_weight_kwargs = NamedTuple()
vis_evol_prfl_modl = (height=200, width_to_time=5, ylims=(0, 0.6), colorrange=nothing)
vis_evol_prfl_axial = (height=200, width_to_time=5, ylims=nothing, colorrange=nothing)
vis_evol_prfl_radial = (height=200, width_to_time=5, ylims=nothing, colorrange=nothing)
vis_evol_prfl_core = vis_evol_prfl_axial
selector_t_hold_prfl_modl = t -> true
selector_pos_prfl_modl = k -> 0.2 < k < 0.4
selector_t_hold_prfl_axial = t -> true
selector_pos_prfl_axial = t -> true
selector_t_hold_prfl_radial = t -> true
selector_pos_prfl_radial = x -> true
quantile_mask_prfl = 0.05
thres_frac_bot_mask_prfl = 0.1

##
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
copy_and_include = (name_script) -> begin
    path_script = joinpath(@__DIR__, name_script)
    cp(path_script, joinpath(path_output, basename(path_script)); force=true)
    include(path_script)
end
"load_excitation_extr.jl" |> copy_and_include
# "load_excitation_corr.jl" |> copy_and_include
"anlz_excitation_corr.jl" |> copy_and_include
"anlz_excitation_vslz_corr.jl" |> copy_and_include
# "anlz_excitation_vslz_extr.jl" |> copy_and_include
