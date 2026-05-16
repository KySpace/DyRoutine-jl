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
    (date="0513", runids=71:75, IB=5.378, tag_head="ImbaEvol", rep_each=1, bias=0.1:0.05:0.6, t_hold=6:10:116),
    (date="0513", runids=76:80, IB=5.376, tag_head="ImbaEvol", rep_each=1, bias=0.1:0.05:0.6, t_hold=6:10:116),
    (date="0513", runids=81, IB=[5.392, 5.386], tag_head="ImbaEvol", rep_each=1, tbias=0.1:0.1:0.6, t_hold=6:10:116),
    (date="0513", runids=82, IB=[5.378, 5.372], tag_head="ImbaEvol", rep_each=1, tbias=0.1:0.1:0.6, t_hold=6:10:116),
]

title_anlz = "[05.15].01.DevTest"
runinfo = runinfos[1]
str_runids = runinfo.runids |> a -> "$(a)" |> s -> replace(s, ":" => "-")
n_ib = runinfo.IB |> length
name_dims = ["IB" "repeat" "bias" "t_hold" "istp"]
