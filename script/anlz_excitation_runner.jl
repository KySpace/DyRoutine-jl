using HDF5
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

# With density profiles drawn, extractions as well
# commit e794fd69f788f0ae5ce2a74503900d8542567941
title_anlz = "[07.24].100.Extr"

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
            ("0325", 67),
            ("0325", 96),
            ("0325", 68),
            ("0325", 50),
            ("0325", 81),
            ("0325", 51),
            ("0325", 79),
            ("0325", 53),
        ],
        vars=(
            IB=[5.311, 5.313, 5.316, 5.318, 5.322, 5.325, 5.326, 5.328, 5.332, 5.333, 5.336, 5.338],
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
ids_runinfo = 1:2
sel_vars = NamedTuple()
# sel_vars = (; t_hold=t -> 0 .<= t .<= 100)
# sel_vars = (; IB=b -> 5.316 .<= b .<= 5.318, t_hold=t -> 0 .<= t .<= 20)


path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

wh_corner = (10, 10)
smwh_roi = (50, 100)
smwh_essn = (30, 60)
smwh_core = (30, 60)
wh_peak = smwh_roi .* 2 .+ 1
smw_peak, smh_peak = smwh_roi
px_in_um = 6.5 / 22.06
len_avg_peak = 10

step_posi = px_in_um
step_modl = 1 ./ (2 .* smwh_core .* px_in_um)
x_modl, y_modl = smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_modl
x_posi, y_posi = smwh_roi |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_posi

idx_IB_axis = 1
idx_rep_axis = 2
idx_t_hold_axis = 3
idx_istp_axis = 4
n_pca_modes = 16
n_pca_modes_prfl_modl = 8
freq_query = 1:1:140
freq_query_pca = 1:1:140
freq_query_pca_modl = 1:1:100

proc_sidepeak = true
proc_envelope = true
selector_moment = y -> (y .> 0.10) .& (y .< 0.50)
selector_sidepeak = y -> (y .> 0.1) .& (y .< 0.5)
selector_t_spectrum = (;
    number=t -> 0 .< t .< 300,
    sp_weight=t -> 0 .< t .< 300,
    sp_height=t -> 0 .< t .< 300,
    sp_width=t -> 0 .< t .< 300,
    sp_wavenum=t -> 0 .< t .< 300,
    nvlp=t -> 0 .< t .< 300,
)
selector_t_pca_dens = t -> 20 .< t .< 80
selector_t_pca_modl = t -> 30 .< t .< 170
selector_tail_sidepeak = y -> y .> 0.2
filter_core_pca_sigma = 1.5
filter_core_pca = im -> imfilter(im, Kernel.gaussian(filter_core_pca_sigma))

mask_modl = (;
    fringe=(x, y) -> ((x - 0.37) / 0.16)^2 + ((y - 0.32) / 0.1)^2 < 1,
    center=(x, y) -> begin
        θ = -20 * π / 180
        x_r = cos(θ) * x - sin(θ) * y
        y_r = sin(θ) * x + cos(θ) * y
        ((x_r / 0.4)^2 + (y_r / 0.15)^2 < 1) || ((x / 0.55)^2 + (y / 0.1)^2 < 1)
    end,
    sidepeak=(x, y) -> abs(x) < 0.51,
    main=(x, y) -> abs(x) < 0.51,
)

fit_stack_kwargs = NamedTuple()
fit_tailess_kwargs = NamedTuple()
fit_asymm_kwargs = (;
    preprocess=ds -> imfilter(ds, Kernel.gaussian(3)),
    θ_hint=(
        max=0.0 / 180 * π,
        min=0.0 / 180 * π,
        init=0.0 / 180 * π
    )
)
fit_round_kwargs = NamedTuple()
query_weight_kwargs = NamedTuple()
trend_property_specs = [
    (
        name="number",
        ylabel="density sum",
        ylim=nothing,
        selection_key="t_vec_sel_number",
        overlay_evol_col=1,
        fit_evol=nothing,
        variants=[(name="dens-sum", evol_spct=("all", "sel"), color=:theme, label="sum", extra=false)],
    ),
    (
        name="weight",
        ylabel="side peak \nweight",
        ylim=(-0.02, 0.17),
        selection_key="t_vec_sel_sp_weight",
        overlay_evol_col=1,
        fit_evol=(
            model=:oscillation_decay,
            kwargs=(
                C_hint=(max=0.10, min=0.02, init=0.08),
                A_hint=(max=0.07, min=0.00, init=0.05),
            ),
        ),
        variants=[
            (name="fit-weight", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
            (name="moment-weight", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
        ],
    ),
    (
        name="height",
        ylabel="side peak \nheight",
        ylim=(-0.1, 1.1),
        selection_key="t_vec_sel_sp_height",
        overlay_evol_col=1,
        fit_evol=(
            model=:oscillation_decay,
            kwargs=(
                C_hint=(max=0.7, min=0.10, init=0.5),
                A_hint=(max=0.5, min=0.00, init=0.5),
            ),
        ),
        variants=[
            (name="fit-height", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
            (name="moment-height", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
        ],
    ),
    (
        name="width",
        ylabel="side peak \nwidth (μm⁻¹)",
        ylim=(0.02, 0.205),
        selection_key="t_vec_sel_sp_width",
        overlay_evol_col=1,
        fit_evol=(
            model=:oscillation_decay,
            kwargs=(
                C_hint=(max=0.10, min=0.02, init=0.05),
                A_hint=(max=0.10, min=0.00, init=0.05),
            ),
        ),
        variants=[
            (name="fit-width", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
            (name="moment-width", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
        ],
    ),
    (
        name="wavenum",
        ylabel="side peak \nwavenum (μm⁻¹)",
        ylim=(0.22, 0.38),
        selection_key="t_vec_sel_sp_wavenum",
        overlay_evol_col=1,
        fit_evol=(
            model=:oscillation_decay,
            kwargs=(
                C_hint=(max=0.30, min=0.25, init=0.28),
                A_hint=(max=0.07, min=0.00, init=0.04),
            ),
        ),
        variants=[
            (name="fit-wavenum", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
            (name="moment-wavenum", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
        ],
    ),
    (
        name="nvlp-size-axial",
        ylabel="envelope size\naxial (μm)",
        ylim=(5, 8),
        selection_key="t_vec_sel_nvlp_size",
        overlay_evol_col=1,
        fit_evol=(
            model=:oscillation_decay,
            kwargs=(
                C_hint=(max=8, min=5, init=7),
                A_hint=(max=5, min=0, init=1),
            ),
        ),
        variants=[
            (name="fit-size-y", evol_spct=("all", "sel"), color=:variant_high, label="fit size axial", extra=false),
        ],
    ),
    (
        name="nvlp-size-radial",
        ylabel="envelope size\nradial (μm)",
        ylim=(0.5, 3.5),
        selection_key="t_vec_sel_nvlp_size",
        overlay_evol_col=1,
        fit_evol=(
            model=:oscillation_decay,
            kwargs=(
                C_hint=(max=3.5, min=0.5, init=2.5),
                A_hint=(max=2.0, min=0.0, init=1.0),
            ),
        ),
        variants=[
            (name="fit-size-x", evol_spct=("all", "sel"), color=:variant_low, label="fit size radial", extra=false),
        ],
    ),
    (
        name="nvlp-center",
        ylabel="envelope center \n (μm)",
        ylim=(-10, 10),
        selection_key="t_vec_sel_nvlp_cent",
        overlay_evol_col=2,
        fit_evol=nothing,
        variants=[
            (name="fit-cent-x", evol_spct=("all", "sel"), color=:variant_low, label="fit cent radial", extra=false),
            (name="fit-cent-y", evol_spct=("all", "sel"), color=:variant_high, label="fit cent axial", extra=false),
        ],
    ),
]
trend_panel_per_IB_kwargs = (width_evol=400, width_spct=400, height=200)
trend_panel_per_prop_kwargs = (width_evol=400, width_spct=400, height=120)
trend_all_IB_groups = (:stacked, :all)
trend_spectrum_IB_groups = (:stacked, :all)
trend_spectrum_IB_kwargs = (width=360, height=180)
trend_spectrum_IB_plot_kwargs = (colorrange=(0.3, 1.00),)
plot_corr_figures = true
plot_extr_figures = false
draw_solo_modl_kwargs = NamedTuple()
vis_evol_prfl_modl = (height=400, width_to_time=10, ylims=(0, 0.6), colorrange=nothing)
vis_evol_prfl_axial = (height=400, width_to_time=10, ylims=nothing, colorrange=nothing)
vis_evol_prfl_radial = (height=400, width_to_time=10, ylims=nothing, colorrange=nothing)
vis_evol_prfl_core = vis_evol_prfl_axial
quantile_mask_prfl = 0.05
thres_frac_bot_mask_prfl = 0.1
selector_t_hold_prfl_modl = t -> true
selector_pos_prfl_modl = k -> 0.2 < k < 0.4
selector_t_hold_prfl_axial = t -> true
selector_pos_prfl_axial = x -> true
selector_t_hold_prfl_radial = t -> true
selector_pos_prfl_radial = x -> true

##
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
copy_and_include = (name_script) -> begin
    path_script = joinpath(@__DIR__, name_script)
    cp(path_script, joinpath(path_output, basename(path_script)); force=true)
    include(path_script)
end



for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    # global tag = gen_run_tag(runinfo)
    global tag = tag_head = runinfo.tag_head
    println("Processing: $tag")

    "anlz_excitation_extr.jl" |> copy_and_include
    "anlz_excitation_corr.jl" |> copy_and_include
    "anlz_excitation_vslz_corr.jl" |> copy_and_include
    # "anlz_excitation_vslz_extr.jl" |> copy_and_include
end
