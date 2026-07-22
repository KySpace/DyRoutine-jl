using HDF5
using MAT
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

# axial and radial profiles integrated into saved data
# commit f895ed70eb888f4caa9008041f63d776e1fc772d
title_anlz = "Anlz.16.Simu-03.[2025.07.22]"

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"03.[2026.06.10]"
istp = ["162", "164"]
unit_t = 0.1798
unit_in_um = 0.2613
a22 = [80, 85, 90, 95, 100, 105]

runinfos = [
    (
        tag_head="SIMU-NTRC",
        date="2026.06.10",
        runids=eachindex(a22),
        bind_id=:IB,
        dir=dir_test,
        vars=(
            IB=a22,
            rep=1:1,
            istp,
        ),
    ),
]
ids_runinfo = eachindex(runinfos)

sel_vars = (; t_hold=(; index=i -> isodd(i)))

path_output = joinpath(path_root, title_anlz)
isdir(path_output) || mkpath(path_output)

smwh_roi = (127, 76)
smwh_essn = (127, 76)
smwh_core = (127, 76)
xy_peak_px_fixed = (128, 128)
xy_peak_core_fixed = smwh_roi .+ 1
len_avg_peak = 10

fmt_probe = format_dens_simulation_runinfo(
    first(runinfos);
    path_root,
    smwh_roi,
    xy_peak_px=xy_peak_px_fixed,
    unit_t,
    unit_in_um,
    sel_vars,
)
px_in_um = fmt_probe.px_in_um
unit_x = fmt_probe.unit_x
unit_y = fmt_probe.unit_y

step_posi = px_in_um
step_modl = 1 ./ (2 .* smwh_core .* px_in_um)
x_modl, y_modl = smwh_core |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_modl
x_posi, y_posi = smwh_roi |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_posi

format_dens_runinfo_kwargs = (;
    path_root,
    smwh_roi,
    xy_peak_px=xy_peak_px_fixed,
    unit_t,
    unit_in_um,
    sel_vars,
)
format_dens_runinfo_fn = format_dens_simulation_runinfo

idx_IB_axis = 1
idx_rep_axis = 2
idx_t_hold_axis = 3
idx_istp_axis = 4
n_pca_modes = min(16, max(1, length(fmt_probe.val_vars.rep) * length(fmt_probe.val_vars.t_hold) - 1))
n_pca_modes_prfl_modl = min(8, max(1, length(fmt_probe.val_vars.t_hold) - 1))
freq_query = 1:1:140
freq_query_pca = 1:1:140
freq_query_pca_modl = 1:1:100

proc_sidepeak = true
proc_envelope = true
selector_moment = y -> (y .> 0.10) .& (y .< 0.60)
selector_sidepeak = y -> (y .> 0.1) .& (y .< 0.6)
selector_t_spectrum = (;
    number=t -> 30 .< t .< 100,
    sp_weight=t -> 30 .< t .< 100,
    sp_height=t -> 30 .< t .< 100,
    sp_width=t -> 30 .< t .< 100,
    sp_wavenum=t -> 30 .< t .< 100,
    nvlp=t -> 30 .< t .< 200,
)
selector_t_pca_dens = t -> 30 .< t .< 80
selector_t_pca_modl = t -> 30 .< t .< 100
selector_tail_sidepeak = y -> y .> 0.2
filter_core_pca = im -> imfilter(im, Kernel.gaussian(1.5))

mask_modl = (;
    fringe=(x, y) -> false,
    center=(x, y) -> false,
    sidepeak=(x, y) -> abs(x) < 0.7,
    main=(x, y) -> abs(x) < 0.7,
)

fit_stack_kwargs = NamedTuple()
fit_tailess_kwargs = NamedTuple()
fit_asymm_kwargs = (;
    preprocess=ds -> imfilter(ds, Kernel.gaussian(3)),
    θ_hint=(
        max=0.0 / 180 * π,
        min=0.0 / 180 * π,
        init=0.0 / 180 * π,
    ),
)
fit_round_kwargs = NamedTuple()
query_weight_kwargs = NamedTuple()
plot_prfl_modl_evol_kwargs = (; colorrange=(0, 3.0))
prfl_axial_halfwidth_um = 16.0
prfl_radial_halfwidth_um = 4.0
plot_prfl_axial_evol_kwargs = (; pos_lims=(-prfl_axial_halfwidth_um, prfl_axial_halfwidth_um))
plot_prfl_radial_evol_kwargs = (; pos_lims=(-prfl_radial_halfwidth_um, prfl_radial_halfwidth_um))
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
        ylim=(-0.02, 0.32),
        selection_key="t_vec_sel_sp_weight",
        overlay_evol_col=1,
        fit_evol=nothing,
        variants=[
            (name="fit-weight", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
            (name="moment-weight", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
        ],
    ),
    (
        name="height",
        ylabel="side peak \nheight",
        ylim=(-0.1, 3.6),
        selection_key="t_vec_sel_sp_height",
        overlay_evol_col=1,
        fit_evol=nothing,
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
        fit_evol=nothing,
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
        fit_evol=nothing,
        variants=[
            (name="fit-wavenum", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
            (name="moment-wavenum", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
        ],
    ),
    (
        name="nvlp-size-axial",
        ylabel="envelope size\naxial (μm)",
        ylim=(6, 10),
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
        ylim=(0.5, 1.5),
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
draw_solo_modl_kwargs = (; dens_max=64.0, peak_height_max=3.0)

cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
copy_and_include = (name_script) -> begin
    path_script = joinpath(@__DIR__, name_script)
    cp(path_script, joinpath(path_output, basename(path_script)); force=true)
    include(path_script)
end

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag = runinfo.tag_head
    println("Processing: $tag")
    println("  [$tag] px_in_um=$(px_in_um), selected t_hold=$(extrema(fmt_probe.val_vars.t_hold)) ms, n_t=$(length(fmt_probe.val_vars.t_hold))")

    "anlz_excitation_extr.jl" |> copy_and_include
    "anlz_excitation_corr.jl" |> copy_and_include
    "anlz_excitation_vslz_corr.jl" |> copy_and_include
    # "anlz_excitation_vslz_extr.jl" |> copy_and_include
end
