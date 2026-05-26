using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
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

path_runner = @__FILE__
path_anlz_routine = joinpath(@__DIR__, "anlz_routine.jl")

year_test = 2026
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
istp = ["162", "164"]
runinfos = [
    (date="0325", runid=95, tag_head="CFNM", vars=(IB=5.311, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=82, tag_head="CFNM", vars=(IB=5.313, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=52, tag_head="CFNM", vars=(IB=5.316, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=80, tag_head="CFNM", vars=(IB=5.318, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=96, tag_head="CFNM", vars=(IB=5.321, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=67, tag_head="CFNM", vars=(IB=5.322, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=68, tag_head="CFNM", vars=(IB=5.328, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=50, tag_head="CFNM", vars=(IB=5.328, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=81, tag_head="CFNM", vars=(IB=5.332, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=51, tag_head="CFNM", vars=(IB=5.333, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=79, tag_head="CFNM", vars=(IB=5.336, rep=1:3, t_hold=6:2:200, istp)),
    (date="0325", runid=53, tag_head="CFNM", vars=(IB=5.338, rep=1:3, t_hold=6:2:200, istp)),
    (date="0322", runid=29, tag_head="NTRC", vars=(IB=5.314, rep=1:3, t_hold=6:2:200, istp)),
    (date="0322", runid=28, tag_head="NTRC", vars=(IB=5.316, rep=1:3, t_hold=6:2:200, istp)),
    (date="0322", runid=27, tag_head="NTRC", vars=(IB=5.318, rep=1:3, t_hold=6:2:200, istp)),
    (date="0322", runid=26, tag_head="NTRC", vars=(IB=5.322, rep=1:3, t_hold=6:2:200, istp)),
    (date="0322", runid=25, tag_head="NTRC", vars=(IB=5.326, rep=1:3, t_hold=6:2:200, istp)),
    (date="0323", runid=61, tag_head="NTRC", vars=(IB=5.332, rep=1:3, t_hold=6:2:200, istp)),
    (date="0323", runid=62, tag_head="NTRC", vars=(IB=5.336, rep=1:3, t_hold=6:2:200, istp)),
    (date="0323", runid=63, tag_head="NTRC", vars=(IB=5.340, rep=1:3, t_hold=6:2:200, istp)),
    (date="0323", runid=64, tag_head="NTRC", vars=(IB=5.343, rep=1:3, t_hold=6:2:200, istp)),
]
ids_runinfo = 1:5

title_anlz = "[05.26].51.DevTests"
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

cp(path_runner, joinpath(path_output, basename(path_runner)); force=true)
cp(path_anlz_routine, joinpath(path_output, basename(path_anlz_routine)); force=true)

wh_corner = (10, 10)
smwh_roi = (30, 60)
smwh_core = smwh_roi
wh_peak = smwh_roi .* 2 .+ 1
smw_peak, smh_peak = smwh_roi
smw_ft = 5
px_in_um = 6.5 / 22.06
len_avg_peak = 10

step_posi = px_in_um
step_modl = 1 / (2 * smwh_roi[2] * px_in_um)
x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi
x_modl, y_modl = (x_vec, y_vec) .* step_modl

idx_ib_axis = 1
idx_rep_axis = 2
idx_t_hold_axis = 3
idx_istp_axis = 4
idx_ib = 1
n_pca_modes = 8
freq_query = 1:1:100

proc_sidepeak = true
proc_envelope = true
selector_moment = y -> (y .> 0.10) .& (y .< 0.50)
selector_sidepeak = y -> (y .> 0.1) .& (y .< 0.5)
selector_t_sidepeak = t -> 0 .< t .< 20
selector_t_envelope = t -> 0 .< t .< 20
selector_tail_stack = y -> y .> 0.02

fit_stack_kwargs = NamedTuple()
fit_tailess_kwargs = NamedTuple()
fit_asymm_kwargs = NamedTuple()
fit_round_kwargs = NamedTuple()
query_weight_kwargs = NamedTuple()

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    include(path_anlz_routine)
end
