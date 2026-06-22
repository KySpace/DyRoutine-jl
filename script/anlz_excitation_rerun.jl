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

# commit 180909b483deb890c66bd7ddd72d5223509bd2e5
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
tag = "CFNM"
title_load = "[06.21].92.Extr.Table.Nvlp.SmallBlur.NoRot"
title_anlz = "[06.21].94.Demo.Mask"

path_load = joinpath(path_root, "AnlzRoutine", title_load)
path_output = joinpath(path_root, "AnlzRoutine", title_anlz)
isdir(path_output) || mkpath(path_output)

path_load_extr = joinpath(path_load, @sprintf("%s_essn_extr.jld2", tag))
path_load_corr = joinpath(path_load, @sprintf("%s_corr.jld2", tag))

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
            kwargs = (
                C_hint=(max=0.10, min=0.02, init=0.08),
                A_hint=(max=0.07, min=0.00, init=0.05),
            )
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
            kwargs = (
                C_hint=(max=0.7, min=0.10, init=0.5),
                A_hint=(max=0.5, min=0.00, init=0.5),
            )
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
            kwargs = (
                C_hint=(max=0.10, min=0.02, init=0.05),
                A_hint=(max=0.10, min=0.00, init=0.05),
            )
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
            kwargs = (
                C_hint=(max=0.30, min=0.25, init=0.28),
                A_hint=(max=0.07, min=0.00, init=0.04),
            )
        ),
        variants=[
            (name="fit-wavenum", evol_spct=("all", "sel"), color=:fit, label="fit", extra=false),
            (name="moment-wavenum", evol_spct=("all", "sel"), color=:moment, label="moment", extra=true),
        ],
    ),
    (
        name="nvlp-size",
        ylabel="envelope \nsize (μm)",
        ylim=(1, 8),
        selection_key="t_vec_sel_nvlp_size",
        overlay_evol_col=2,
        fit_evol=(
            model=:oscillation_decay,
            kwargs = (
                C_hint=(max=8, min=1, init=4),
                A_hint=(max=5, min=0, init=1),
            )
        ),
        variants=[
            (name="fit-size-x", evol_spct=("all", "sel"), color=:variant_low, label="fit size radial", extra=false),
            (name="fit-size-y", evol_spct=("all", "sel"), color=:variant_high, label="fit size axial", extra=false),
        ],
    ),
    (
        name="nvlp-center",
        ylabel="envelope center\n (μm)",
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
selector_t_pca_dens = t -> 20 .< t .< 80
selector_t_pca_modl = t -> 00 .< t .< 80
freq_query_pca_modl = 1:1:100
filter_core_pca = im -> imfilter(im, Kernel.gaussian(1.5))
query_weight_kwargs = NamedTuple()
draw_solo_modl_kwargs = NamedTuple()

cp(@__FILE__, joinpath(path_output, basename(@__FILE__)); force=true)
copy_and_include = (name_script) -> begin
    path_script = joinpath(@__DIR__, name_script)
    cp(path_script, joinpath(path_output, basename(path_script)); force=true)
    include(path_script)
end
"load_excitation_extr.jl" |> copy_and_include
# "anlz_excitation_corr.jl" |> copy_and_include
# "load_excitation_corr.jl" |> copy_and_include
# "anlz_excitation_vslz_corr.jl" |> copy_and_include
# "anlz_excitation_vslz_extr.jl" |> copy_and_include
