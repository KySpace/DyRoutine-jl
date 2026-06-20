log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))

path_load_corr = @isdefined(path_load_corr) ? path_load_corr : joinpath(path_load, @sprintf("%s_corr.jld2", tag_load))
t_stage = log_step("loading excitation correlation cache from $path_load_corr")
cache_corr = JLD2.load(path_load_corr)
meta_corr = cache_corr["meta_corr"]

tag = meta_corr.tag
path_output = @isdefined(path_output) ? path_output : meta_corr.path_output
runinfo = meta_corr.runinfo
val_vars = meta_corr.val_vars
name_dims = meta_corr.name_dims
n_dim_vars = meta_corr.n_dim_vars
n_dim_vars_per_IB = meta_corr.n_dim_vars_per_IB
n_IB = meta_corr.n_IB
n_rep = meta_corr.n_rep
n_main = meta_corr.n_main
n_istp = meta_corr.n_istp
n_pca_modes = meta_corr.n_pca_modes
px_in_um = meta_corr.px_in_um
smwh_roi = meta_corr.smwh_roi
smwh_core = meta_corr.smwh_core
y_modl = meta_corr.y_modl
x_posi = meta_corr.x_posi
y_posi = meta_corr.y_posi
freq_query = meta_corr.freq_query
freq_query_pca = meta_corr.freq_query_pca
tag_IBs = meta_corr.tag_IBs
runinfo_plots = meta_corr.runinfo_plots
info_fmt = meta_corr.info_fmt
info_fitting = meta_corr.info_fitting
trend_property_specs = meta_corr.trend_property_specs
trend_panel_per_IB_kwargs = meta_corr.trend_panel_per_IB_kwargs
trend_panel_per_prop_kwargs = meta_corr.trend_panel_per_prop_kwargs
trend_all_IB_groups = meta_corr.trend_all_IB_groups
trend_spectrum_IB_groups = meta_corr.trend_spectrum_IB_groups
trend_spectrum_IB_kwargs = meta_corr.trend_spectrum_IB_kwargs
trend_spectrum_IB_plot_kwargs = meta_corr.trend_spectrum_IB_plot_kwargs

trend_sidepeak_nvlp = cache_corr["trend_sidepeak_nvlp"]
trend_extr_stacked_over_rep = cache_corr["trend_extr_stacked_over_rep"]
trend_stacked_over_rep = cache_corr["trend_stacked_over_rep"]
if haskey(cache_corr, "prfl_evol_stacked")
    prfl_evol = cache_corr["prfl_evol"]
    prfl_evol_stacked = cache_corr["prfl_evol_stacked"]
elseif haskey(cache_corr, "prfl_evol_reps")
    prfl_evol = cache_corr["prfl_evol_reps"]
    prfl_evol_stacked = cache_corr["prfl_evol"]
else
    prfl_evol = nothing
    prfl_evol_stacked = cache_corr["prfl_evol"]
end
modes_pca_dens2d = cache_corr["modes_pca_dens2d"]
pca_spectra = cache_corr["pca_spectra"]
log_done("loaded excitation correlation cache", t_stage)
