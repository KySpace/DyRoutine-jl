log_step(msg) = (println("  [$tag] $msg"); flush(stdout); time())
log_done(msg, t_start) = (println("  [$tag] $msg ($(round(time() - t_start; digits=1)) s)"); flush(stdout))

path_load_extr = @isdefined(path_load_extr) ? path_load_extr : joinpath(path_load, @sprintf("%s_essn_extr.jld2", tag_load))
t_stage = log_step("loading excitation extraction cache from $path_load_extr")
cache_extr = JLD2.load(path_load_extr)
meta_extr = cache_extr["meta_extr"]

tag = meta_extr.tag
runinfo = meta_extr.runinfo
val_vars = meta_extr.val_vars
name_dims = meta_extr.name_dims
n_dim_vars = meta_extr.n_dim_vars
n_dim_vars_per_IB = meta_extr.n_dim_vars_per_IB
n_IB = meta_extr.n_IB
n_rep = meta_extr.n_rep
n_main = meta_extr.n_main
n_istp = meta_extr.n_istp
n_pca_modes = meta_extr.n_pca_modes
px_in_um = meta_extr.px_in_um
smwh_roi = meta_extr.smwh_roi
smwh_core = meta_extr.smwh_core
y_modl = meta_extr.y_modl
x_posi = meta_extr.x_posi
y_posi = meta_extr.y_posi
freq_query = meta_extr.freq_query
freq_query_pca = meta_extr.freq_query_pca
tag_IBs = meta_extr.tag_IBs
runinfo_plots = meta_extr.runinfo_plots
info_fmt = meta_extr.info_fmt
info_fitting = meta_extr.info_fitting
trend_property_specs = meta_extr.trend_property_specs
trend_panel_per_IB_kwargs = meta_extr.trend_panel_per_IB_kwargs
trend_panel_per_prop_kwargs = meta_extr.trend_panel_per_prop_kwargs
trend_all_IB_groups = meta_extr.trend_all_IB_groups
trend_spectrum_IB_groups = meta_extr.trend_spectrum_IB_groups
trend_spectrum_IB_kwargs = meta_extr.trend_spectrum_IB_kwargs
trend_spectrum_IB_plot_kwargs = meta_extr.trend_spectrum_IB_plot_kwargs

essn_2d_fmt = cache_extr["essn_2d_fmt"]
essn_stacked_over_rep = cache_extr["essn_stacked_over_rep"]
fit_prfl_modl_over_rep_1d = cache_extr["fit_prfl_modl_over_rep_1d"]
extr_fmt = cache_extr["extr_fmt"]
extr_stacked_over_rep = cache_extr["extr_stacked_over_rep"]
log_done("loaded excitation extraction cache", t_stage)
