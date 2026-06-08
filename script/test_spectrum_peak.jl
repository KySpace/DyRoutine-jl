using HDF5
using CairoMakie: Figure, Axis, Colorbar, DataAspect, heatmap!, lines!, scatter!, save, text!, rowgap!, colgap!
using GLMakie
using JLD2
using Printf
using ImageFiltering
include(joinpath(@__DIR__, "..", "src", "helper.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "loadfmt.jl"))
include(joinpath(@__DIR__, "..", "src", "percond.jl"))
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "corr.jl"))
include(joinpath(@__DIR__, "..", "src", "vissolo.jl"))
include(joinpath(@__DIR__, "..", "src", "viscorr.jl"))
include(joinpath(@__DIR__, "..", "src", "vispca.jl"))
path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations"
title_anlz = "[06.08].64.PCATests.[t=0-100ms]"

path_pca = joinpath(path_root, "AnlzRoutine", title_anlz, "PCA modes")
path_pca_data = joinpath(path_root, "AnlzRoutine", title_anlz, "CFNM_pca_modes.jld2")
@load path_pca_data modes_pca_dens2d runinfo val_vars n_pca_modes freq_query_pca px_in_um smwh_core pca_spectra
selector_t_pca = t -> 20 .< t .< 80

tag = "CFNM"
get_bind_date(runinfo, idx_bind) = hasproperty(runinfo, :date_runid) ? first(runinfo.date_runid[idx_bind]) : runinfo.date
get_bind_runid(runinfo, idx_bind) =
    if hasproperty(runinfo, :date_runid)
        last(runinfo.date_runid[idx_bind])
    elseif hasproperty(runinfo, :runids)
        as_vector(runinfo.runids)[idx_bind]
    else
        as_vector(runinfo.runid)[idx_bind]
    end
get_bind_runinfo(runinfo, val_vars, idx_bind) = merge(
    runinfo,
    (;
        date=get_bind_date(runinfo, idx_bind),
        runid=get_bind_runid(runinfo, idx_bind),
        IB=val_vars.IB[idx_bind],
    ),
)
log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))


function get_spectrum_peaks(freq, spct; min_prom=0.2)
    pks = spct |> findmaxima |> peakproms!(; min=min_prom) |> peakwidths!
    height_max = sum(pks.heights)
    pks_record = map(
        (idx, height) -> (freq=freq[idx], value=height, value_reduced=height/height_max), pks.indices, pks.heights
        ) |> p -> sort(p; by=x -> x.value_reduced, rev=true)
    return pks_record
end


fig_live, axs_live = set_axis!("test peaks")
[axs_live] |> clear_axes!
mode = modes_pca_dens2d[2][15]
spectral_weight = calc_spct_rep_evol(eachslice(mode.weight; dims=1), val_vars.t_hold, freq_query_pca; sel_evol=selector_t_pca)
pks = spectral_weight.spct_mean_mask |> findmaxima |> peakproms!(; min=0.2) |> peakwidths!
rcrd_pks = spectral_weight.spct_mean_mask |> s -> record_peaks(freq_query_pca, s; min_prom=0.2)
peaksplot!(axs_live, freq_query_pca, pks)
fig_live |> display


for (c, IB) in enumerate(val_vars.IB)
    tag_IB = gen_run_tag(get_bind_runinfo(runinfo, val_vars, c))
    runinfo_plot = get_bind_runinfo(runinfo, val_vars, c)

    t_stage = log_step("building and saving PCA figure for $tag_IB")
    isdir(path_pca) || mkpath(path_pca)
    fig_pca_mode = Figure()
    for idx_mode in 1:n_pca_modes
        mode = modes_pca_dens2d[c][idx_mode]
        fig_pca_mode.layout |> clean_gridlayout!
        gl_pca_mode = GridLayout()
        fig_pca_mode[1, 1] = gl_pca_mode
        axs_pca_mode = set_panel_pca_duet!(gl_pca_mode)
        gl_pca_mode[0, 1] = Label(fig_pca_mode, "$tag_IB | #$idx_mode"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
        spectral_weight = calc_spct_rep_evol(eachslice(mode.weight; dims=1), val_vars.t_hold, freq_query_pca; sel_evol=selector_t_pca)
        rcrd_pks = spectral_weight.spct_mean_mask |> s -> record_peaks(freq_query_pca, s; min_prom=0.2)
        plot_mode_evol_spct_duet!(axs_pca_mode, mode, spectral_weight, rcrd_pks, val_vars.istp; step_posi=px_in_um, smwh=smwh_core)
        gl_pca_mode |> l -> rowgap!(l, 0)
        resize_to_layout!(fig_pca_mode)
        fig_pca_mode |> f -> save(joinpath(path_pca, @sprintf("%s_%d.png", tag_IB, idx_mode)), f; backend=CairoMakie)
    end
    log_done("saved PCA figure for $tag_IB", t_stage)
end