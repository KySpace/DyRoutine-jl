using HDF5
using CairoMakie
using GLMakie
using JLD2
using Printf
using Statistics

GLMakie.activate!()

include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))
include(joinpath(@__DIR__, "..", "src", "visduet.jl"))

path_runner = @__FILE__
path_anlz_miscibility = joinpath(@__DIR__, "anlz_miscibility.jl")

year_test = 2026
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\SingDrplMisc"
title_anlz = "[05.26].06.DevTest"
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

cp(path_runner, joinpath(path_output, basename(path_runner)); force=true)
cp(path_anlz_miscibility, joinpath(path_output, basename(path_anlz_miscibility)); force=true)

istp = ["162", "164"]
runinfos = [
    (date="0513", runids=71:75, tag_head="ImbaEvol", vars=(IB=5.378, rep=1:5, bias=0.1:0.05:0.6, t_hold=6:5:56, istp)),
    (date="0513", runids=76:80, tag_head="ImbaEvol", vars=(IB=5.376, rep=1:5, bias=0.1:0.05:0.6, t_hold=6:5:56, istp)),
    (date="0513", runids=81:82, tag_head="ImbaEvol", vars=(IB=[5.392, 5.386, 5.378, 5.372], rep=1:1, bias=0.1:0.1:0.6, t_hold=6:10:116, istp)),
]
ids_runinfo = eachindex(runinfos)
sel_vars = NamedTuple()
# sel_vars = (; IB=(; index=i -> i < 20, val=ib -> 5.305 < ib < 5.355))

wh_corner = (10, 10)
smwh_roi = (30, 30)
smwh_core = (20, 20)
smwh_strip = (2, 20)
wh_peak = smwh_roi .* 2 .+ 1
smw_peak, smh_peak = smwh_roi
smw_ft = 5
px_in_um = 6.5 / 22.06
len_avg_peak = 10

step_posi = px_in_um
step_modl = 1 ./ (2 .* smwh_roi .* px_in_um)
x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi
x_modl, y_modl = (x_vec, y_vec) .* step_modl

idx_rep_axis = 2
idx_istp_axis = 5
n_istp_per_condition = length(istp)

proc_sidepeak = false
proc_envelope = true
selector_moment = y -> (y .> 0.10) .& (y .< 0.50)
selector_sidepeak = y -> (y .> 0.1) .& (y .< 0.5)

fit_round_kwargs = (A_hint=(max=25.0, min=0, init=10.0),)
fit_asymm_kwargs = (
    θ_hint=(max=20.0 / 180 * π, min=-10.0 / 180 * π, init=10.0 / 180 * π),
    A_hint=(max=25.0, min=0, init=10.0),
)

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    include(path_anlz_miscibility)
end
