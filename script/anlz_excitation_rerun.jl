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

# commit
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
tag = "CFNM"
title_load = "[06.20].81.Dev.Save"
title_anlz = "[06.20].86.[←81].Dev.PrflEvol.Lite"

path_load = joinpath(path_root, "AnlzRoutine", title_load)
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

path_load_extr = joinpath(path_load, @sprintf("%s_essn_extr.jld2", tag))
path_load_corr = joinpath(path_load, @sprintf("%s_corr.jld2", tag))

# Recomputed correlation settings. Change these here if the rerun should use
# different analysis choices from the saved extraction metadata.
selector_t_spectrum = (;
    number=t -> 0 .< t .< 120,
    sp_weight=t -> 0 .< t .< 80,
    sp_height=t -> 0 .< t .< 80,
    sp_width=t -> 0 .< t .< 80,
    sp_wavenum=t -> 0 .< t .< 80,
    nvlp=t -> 0 .< t .< 80,
)
selector_t_pca = t -> 20 .< t .< 80
filter_core_pca = im -> imfilter(im, Kernel.gaussian(1.5))
query_weight_kwargs = NamedTuple()

cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
copy_and_include = (name_script) -> begin
    path_script = joinpath(@__DIR__, name_script)
    cp(path_script, joinpath(path_output, basename(path_script)); force=true)
    include(path_script)
end
"load_excitation_extr.jl" |> copy_and_include
"anlz_excitation_corr.jl" |> copy_and_include
# "load_excitation_corr.jl" |> copy_and_include
"anlz_excitation_vslz_corr.jl" |> copy_and_include
# "anlz_excitation_vslz_extr.jl" |> copy_and_include
