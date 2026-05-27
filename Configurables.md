# Configurables

## Run And Data Selection

Name change notice: the excitation workflow scripts formerly named `script/anlz_routine.jl` and `script/anlz_routine_batch.jl` are now `script/anlz_excitation.jl` and `script/anlz_excitation_runner.jl`.

- `year_test`
- `path_root`
  - excitation runner: `raw"...\Data\Excitations"`
  - miscibility: `raw"...\Data\SingDrplMisc"`
- `title_anlz`
- `runinfos`
  - excitation runner: `tag_head`, `date_runid`, `bind_id`, `vars`
  - `date_runid`: tuples of `(date, runid)` used for loading raw folders; its order should match the first formatted variable axis
  - `bind_id`: variable name whose values are one-to-one with run ids, such as `:IB`
  - miscibility: `date`, `runids`, `tag_head`, `vars`
- variable axes:
  - excitation runner: `IB`, `rep`, `t_hold`, `istp`
  - per-IB plot filenames/log tags combine `date_runid` and `vars.IB` by the shared bound-axis index
  - miscibility: `IB`, `rep`, `bias`, `t_hold`, `istp`
- lite/test data slice:
  - `rng_lite = 1:50`

## Geometry And Calibration

- `wh_corner = (10, 10)`
- `smwh_roi`
  - excitation: `smwh_roi = (30, 60)`
  - miscibility: `smwh_roi = (30, 30)`
- `smwh_core = (20, 20)`
- `smwh_strip = (2, 20)`
- `smw_ft = 5`
- `px_in_um = 6.5 / 22.06`

## Solo Essentials And Extraction

- `calc_solo_essn_2d`
  - default `smwh_strip=smwh`
- `calc_solo_extr`
  - default `proc_sidepeak=false`
  - default `proc_envelope=false`
  - `sel_moment = y -> (y .> 0.10) .& (y .< 0.50)`
  - `sel_sidepeak = (y_modl .> 0.1) .& (y_modl .< 0.5)`
- sidepeak stack fit mask:
  - `fit_prfl_modl_twinpeak_decay_1d(..., (y_modl .> 0.02))`
- peak finding:
  - `find_peak_position_moving(...; len_avg=10)`
  - `find_positive_cluster_center(...; len_avg=10)`

## Profile Fit Hints

- `fit_prfl_modl_twinpeak_decay_1d`
  - `M_hint=(max=Inf, min=2.0, init=3.0)`
  - `σ0_hint=(max=0.30, min=0.02, init=0.1)`
  - `P_hint=(max=2.0, min=0.0, init=0.5)`
  - `σ_hint=(max=0.100, min=0.018, init=0.05)`
  - `p_hint=(max=0.37, min=0.23, init=0.3)`
  - `D_hint=(max=Inf, min=0.0, init=0.5)`
  - `λ_hint=(max=5.0, min=0.5, init=0.8)`
- `fit_prfl_modl_twinpeak_1d`
  - `M_hint=(max=Inf, min=2.0, init=3.0)`
  - `σ0_hint=(max=0.30, min=0.05, init=0.1)`
  - `P_hint=(max=2.0, min=0.0, init=0.5)`
  - `σ_hint=(max=0.200, min=0.018, init=0.05)`
  - `p_hint=(max=0.37, min=0.23, init=0.3)`

## Envelope Fit Hints

- `fit_dens2d_gaussian_elliptic_disk`
  - `θ_hint=(max=20.0/180*π, min=-10.0/180*π, init=10.0/180*π)`
  - `A_hint=(max=25.0, min=0, init=10.0)`
- `fit_dens2d_gaussian_round_disk`
  - amplitude init/upper/lower `[10, 25, 0]`

## PCA And Trend Analysis

- PCA region:
  - modulation sideband row range `smwh_peak[2]+1+8:smwh_peak[2]+1+15`
  - modulation column half-width `smw_ft`
- PCA mode count:
  - `fit_pca_modes(8, m)`
  - script-side `n_mode` values such as `8`
- trend selectors:
  - `selector_t_sidepeak = t -> 0 .< t .< 20`
  - `selector_t_envelope = t -> 0 .< t .< 20`
- frequency query:
  - `freq_query = 1:1:100`
- `query_weight`
  - `scaling = 1000.0`

## Stacking Behavior

- `calc_stacked_essn`
  - averages all numeric fields across `essns`
  - metadata-like fields are taken from `essn_ref`: `smwh`, `smwh_strip`, `smw_modl`, `step_posi`, `step_modl`
  - `offset_cent_core` and `smwh_core` are averaged with `mean_tuple`

<!-- Commented out for now:
- Center finding / duet handling:
  - round Gaussian center fit over `1:wh_peak[1]`, `1:wh_peak[2]`
  - `repeat(..., inner=ntuple(i -> i == 5 ? 2 : 1, length(r.n_dim_vars)))`
- 1D Gaussian center-fit internal heuristics:
  - center-fit width guess `sigma0 = clamp(n / 4, 2.0, float(n))`
  - center-fit amplitude lower bound `amp0 / 100`
  - center-fit sigma lower bound `min(2.0, float(n))`
- 2D Gaussian fit internal scale heuristics:
  - initial size divisor `3`
  - size bounds multipliers `/10` and `*10`
- `calc_solo_essn_2d` internal formulas:
  - modulation step formula `1 / (2 * smwh[2] * px_in_um)`
  - profile normalization `sum(prfl_modl) * step_modl / 2`
- `fit_pca_modes`: `pratio=1.0`
- `query_weight` internals:
  - Fourier phase divisor `1000.0`
  - normalization `weight / sum(weight)`
- `selector_t_sidepeak` and `selector_t_envelope` are externally supplied but good candidates for a trend config object
-->
