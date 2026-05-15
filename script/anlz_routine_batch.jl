using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))

year_test = 2026
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
runinfos = [
    (date="0325", runid=95, IB=5.311, tag_head="CFNM"),
    (date="0325", runid=82, IB=5.313, tag_head="CFNM"),
    (date="0325", runid=52, IB=5.316, tag_head="CFNM"),
    (date="0325", runid=80, IB=5.318, tag_head="CFNM"),
    (date="0325", runid=96, IB=5.321, tag_head="CFNM"),
    (date="0325", runid=67, IB=5.322, tag_head="CFNM"),
    (date="0325", runid=68, IB=5.328, tag_head="CFNM"),
    (date="0325", runid=50, IB=5.328, tag_head="CFNM"),
    (date="0325", runid=81, IB=5.332, tag_head="CFNM"),
    (date="0325", runid=51, IB=5.333, tag_head="CFNM"),
    (date="0325", runid=79, IB=5.336, tag_head="CFNM"),
    (date="0325", runid=53, IB=5.338, tag_head="CFNM"),
    (date="0322", runid=29, IB=5.314, tag_head="NTRC"),
    (date="0322", runid=28, IB=5.316, tag_head="NTRC"),
    (date="0322", runid=27, IB=5.318, tag_head="NTRC"),
    (date="0322", runid=26, IB=5.322, tag_head="NTRC"),
    (date="0322", runid=25, IB=5.326, tag_head="NTRC"),
    (date="0323", runid=61, IB=5.332, tag_head="NTRC"),
    (date="0323", runid=62, IB=5.336, tag_head="NTRC"),
    (date="0323", runid=63, IB=5.340, tag_head="NTRC"),
    (date="0323", runid=64, IB=5.343, tag_head="NTRC"),
]

title_anlz = "[05.15].49.[0-20ms]"
for runinfo_iter in runinfos
    global runinfo = runinfo_iter
    include(joinpath(@__DIR__, "..", "script", "anlz_routine.jl"))
end
