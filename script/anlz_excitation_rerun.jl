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

# commit f895ed70eb888f4caa9008041f63d776e1fc772d
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
tag = "SIMU-NTRC"
title_load = "Anlz.16.Simu-03.[2025.07.22]"
title_anlz = "Anlz.17.Simu-03.[2025.07.22].[←16]"

path_load = joinpath(path_root, title_load)
path_output = joinpath(path_root, title_anlz)
isdir(path_output) || mkpath(path_output)

path_load_extr = joinpath(path_load, @sprintf("%s_essn_extr.jld2", tag))
path_load_corr = joinpath(path_load, @sprintf("%s_corr.jld2", tag))

## Recomputed correlation settings. Change these here if the rerun should use
# different analysis choices from the saved extraction metadata.

selector_t_pca_dens = t -> 30 .< t .< 200
selector_t_pca_modl = t -> 30 .< t .< 200
n_pca_modes_prfl_modl = 8
freq_query_pca_modl = 1:1:100
selector_t_spectrum = (;
    number=t -> 30 .< t .< 100,
    sp_weight=t -> 30 .< t .< 100,
    sp_height=t -> 30 .< t .< 100,
    sp_width=t -> 30 .< t .< 100,
    sp_wavenum=t -> 30 .< t .< 100,
    nvlp=t -> 30 .< t .< 200,
)
filter_core_pca_sigma = 1.5
filter_core_pca = im -> imfilter(im, Kernel.gaussian(filter_core_pca_sigma))
query_weight_kwargs = NamedTuple()
vis_evol_prfl_modl = (height=400, width_to_time=10, ylims=(0, 0.6), colorrange=nothing)
vis_evol_prfl_core = (height=400, width_to_time=10, ylims=nothing, colorrange=nothing)
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
        ylim=(5, 15),
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

##
cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
copy_and_include = (name_script) -> begin
    path_script = joinpath(@__DIR__, name_script)
    cp(path_script, joinpath(path_output, basename(path_script)); force=true)
    include(path_script)
end
"load_excitation_extr.jl" |> copy_and_include
# "load_excitation_corr.jl" |> copy_and_include
"anlz_excitation_corr.jl" |> copy_and_include
"anlz_excitation_vslz_corr.jl" |> copy_and_include
# "anlz_excitation_vslz_extr.jl" |> copy_and_include
