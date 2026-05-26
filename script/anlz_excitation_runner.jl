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
path_anlz_excitation = joinpath(@__DIR__, "anlz_excitation.jl")

year_test = 2026
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
istp = ["162", "164"]
runinfos = [
    (
        tag_head="CFNM",
        date_runid=[
            ("0325", 95),
            ("0325", 82),
            ("0325", 52),
            ("0325", 80),
            ("0325", 96),
            ("0325", 67),
            ("0325", 68),
            ("0325", 50),
            ("0325", 81),
            ("0325", 51),
            ("0325", 79),
            ("0325", 53),
        ],
        vars=(
            runid_IB=[
                (95, 5.311),
                (82, 5.313),
                (52, 5.316),
                (80, 5.318),
                (96, 5.321),
                (67, 5.322),
                (68, 5.328),
                (50, 5.328),
                (81, 5.332),
                (51, 5.333),
                (79, 5.336),
                (53, 5.338),
            ],
            rep=1:3,
            t_hold=6:2:200,
            istp,
        ),
    ),
    (
        tag_head="NTRC",
        date_runid=[
            ("0322", 29),
            ("0322", 28),
            ("0322", 27),
            ("0322", 26),
            ("0322", 25),
            ("0323", 61),
            ("0323", 62),
            ("0323", 63),
            ("0323", 64),
        ],
        vars=(
            runid_IB=[
                (29, 5.314),
                (28, 5.316),
                (27, 5.318),
                (26, 5.322),
                (25, 5.326),
                (61, 5.332),
                (62, 5.336),
                (63, 5.340),
                (64, 5.343),
            ],
            rep=1:3,
            t_hold=6:2:200,
            istp,
        ),
    ),
]
ids_runinfo = eachindex(runinfos)

title_anlz = "[05.26].52.DevTests"
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

cp(path_runner, joinpath(path_output, basename(path_runner)); force=true)
cp(path_anlz_excitation, joinpath(path_output, basename(path_anlz_excitation)); force=true)

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

idx_runid_IB_axis = 1
idx_rep_axis = 2
idx_t_hold_axis = 3
idx_istp_axis = 4
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
    include(path_anlz_excitation)
end
