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

path_runner = @__FILE__
path_load_excitation_extr = joinpath(@__DIR__, "load_excitation_extr.jl")
path_anlz_excitation_corr = joinpath(@__DIR__, "anlz_excitation_corr.jl")
path_load_excitation_corr = joinpath(@__DIR__, "load_excitation_corr.jl")
path_anlz_excitation_vslz = joinpath(@__DIR__, "anlz_excitation_vslz.jl")
path_anlz_excitation_extr = joinpath(@__DIR__, "anlz_excitation_extr.jl")

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
tag_load = get(ENV, "DYROUTINE_EXCITATION_TAG", "CFNM")
mode_rerun = Symbol(get(ENV, "DYROUTINE_EXCITATION_RERUN_MODE", "corr_from_extr"))

if mode_rerun == :corr_from_extr
    title_load = get(ENV, "DYROUTINE_EXCITATION_LOAD_TITLE", "[06.20].81.Dev.Save")
    title_anlz = get(ENV, "DYROUTINE_EXCITATION_OUTPUT_TITLE", "[06.20].82.Dev.Load.Corr")
elseif mode_rerun == :corr_vslz_from_corr
    title_load = get(ENV, "DYROUTINE_EXCITATION_LOAD_TITLE", "[06.20].82.Dev.Load.Corr")
    title_anlz = get(ENV, "DYROUTINE_EXCITATION_OUTPUT_TITLE", "[06.20].83.Dev.Load.Corr.Vslz")
elseif mode_rerun == :extr_vslz_from_extr
    title_load = get(ENV, "DYROUTINE_EXCITATION_LOAD_TITLE", "[06.20].82.Dev.Load.Corr")
    title_anlz = get(ENV, "DYROUTINE_EXCITATION_OUTPUT_TITLE", "[06.20].84.Dev.Load.Extr.Vslz")
else
    throw(ArgumentError("Unknown mode_rerun=$mode_rerun"))
end

path_load = joinpath(path_root, "AnlzRoutine", title_load)
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

for path_script in [
    path_runner,
    path_load_excitation_extr,
    path_anlz_excitation_corr,
    path_load_excitation_corr,
    path_anlz_excitation_vslz,
    path_anlz_excitation_extr,
]
    isfile(path_script) && cp(path_script, joinpath(path_output, basename(path_script)); force=true)
end

tag = tag_load
path_load_extr = joinpath(path_load, @sprintf("%s_essn_extr.jld2", tag_load))
path_load_corr = joinpath(path_load, @sprintf("%s_corr.jld2", tag_load))

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

if mode_rerun == :corr_from_extr
    include(path_load_excitation_extr)
    cp(path_load_extr, joinpath(path_output, basename(path_load_extr)); force=true)
    include(path_anlz_excitation_corr)
elseif mode_rerun == :corr_vslz_from_corr
    include(path_load_excitation_corr)
    plot_corr_figures = true
    plot_extr_figures = false
    include(path_anlz_excitation_vslz)
elseif mode_rerun == :extr_vslz_from_extr
    include(path_load_excitation_extr)
    plot_corr_figures = false
    plot_extr_figures = true
    include(path_anlz_excitation_vslz)
end
