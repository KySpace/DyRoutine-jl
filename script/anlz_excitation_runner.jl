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
# Grouped version for comparing against the one-folder-per-runinfo list above.
runinfos_grouped = [
    (
        tag_head="CFNM",
        bind_id=:IB,
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
            IB=[5.311, 5.313, 5.316, 5.318, 5.321, 5.322, 5.328, 5.328, 5.332, 5.333, 5.336, 5.338],
            rep=1:3,
            t_hold=6:2:200,
            istp,
        ),
    ),
    (
        tag_head="NTRC",
        bind_id=:IB,
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
            IB=[5.314, 5.316, 5.318, 5.322, 5.326, 5.332, 5.336, 5.340, 5.343],
            rep=1:3,
            t_hold=6:2:200,
            istp,
        ),
    ),
]

runinfos = runinfos_grouped
# runinfos = runinfos_separated

# ids_runinfo = eachindex(runinfos)
ids_runinfo = 1:1
# sel_vars = NamedTuple()
sel_vars = (; t_hold=t-> 20 .<= t .<= 75)

title_anlz = "[05.28].54.TrendOverIB"
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

idx_IB_axis = 1
idx_rep_axis = 2
idx_t_hold_axis = 3
idx_istp_axis = 4
n_pca_modes = 8
freq_query = 1:1:100

proc_sidepeak = true
proc_envelope = true
selector_moment = y -> (y .> 0.10) .& (y .< 0.50)
selector_sidepeak = y -> (y .> 0.1) .& (y .< 0.5)
selector_t_sidepeak = t -> 20 .< t .< 80
selector_t_envelope = t -> 20 .< t .< 80
selector_tail_stack = y -> y .> 0.02

fit_stack_kwargs = NamedTuple()
fit_tailess_kwargs = NamedTuple()
fit_asymm_kwargs = NamedTuple()
fit_round_kwargs = NamedTuple()
query_weight_kwargs = NamedTuple()

##

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    # global tag = gen_run_tag(runinfo)
    global tag = tag_head = runinfo.tag_head
    println("Processing: $tag")
    include(path_anlz_excitation)
end
