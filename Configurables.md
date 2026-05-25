# Configurables
runinfos
title_anlz

wh_corner = (10, 10)
smwh_roi = (30, 30)
smwh_core = (20, 20)
smwh_strip = (2, 20)

smw_ft = 5
px_in_um = 6.5 / 22.06

calc_solo_extr
    sel_moment = y -> (y .> 0.10) .& (y .< 0.50)
    sel_sidepeak = (y_modl .> 0.1) .& (y_modl .< 0.5)

- corr.jl
  - query_weight
    - scaling

- persolo.jl
  - fit_prfl_modl_twinpeak_decay_1d

## Additional candidates found

### script/anlz_routine.jl
- Run/data setup:
  - `year_test = 2026`
  - `path_root = raw"...\Data\Excitations"`
  - default/commented `rep = 1:3`, `t_hold = 6:2:200`, `istp = ["162", "164"]`
  - `runinfos` entries: `date`, `runid`, `IB`, `tag_head`
  - `title_anlz`
- Geometry/calibration:
  - `wh_corner = (10, 10)`
  - `smwh_peak = (30, 60)`
  - `smw_ft = 5`
  - `px_in_um = 6.5 / 22.06`
- Sidepeak stack fit mask:
  - `fit_prfl_modl_twinpeak_decay_1d(..., (y_modl .> 0.02))`
- PCA region and mode count:
  - modulation sideband row range `smwh_peak[2]+1+8:smwh_peak[2]+1+15`
  - modulation column half-width `smw_ft`
  - `fit_pca_modes(8, m)`
- Trend analysis:
  - `selector_t_sidepeak = t -> 0 .< t .< 20`
  - `selector_t_envelope = t -> 0 .< t .< 20`
  - `freq_query = 1:1:100`
- Lite/test data slice:
  - `rng_lite = 1:50`

### script/anlz_routine_batch.jl
- Batch run setup:
  - `year_test = 2026`
  - `path_root = raw"...\Data\Excitations"`
  - `rep = 1:3`
  - `t_hold = 6:2:200`
  - `istp = ["162", "164"]`
  - `runinfos` entries: `date`, `runid`, `IB`, `tag_head`
  - `title_anlz = "[05.20].50.DevTests"`

### script/anlz_miscibility.jl
- Run/data setup:
  - `year_test = 2026`
  - `path_root = raw"...\Data\SingDrplMisc"`
  - `istp = ["162", "164"]`
  - `runinfos` entries: `date`, `runids`, `tag_head`, and `vars`
  - `vars`: `IB`, `rep`, `bias`, `t_hold`, `istp`
  - `title_anlz = "[05.21].05.StackedDuet"`
- Geometry/calibration:
  - `wh_corner = (10, 10)`
  - `smwh_roi = (30, 30)`
  - `smwh_core = (20, 20)`
  - `smwh_strip = (2, 20)`
  - `smw_ft = 5`
  - `px_in_um = 6.5 / 22.06`
- Center finding / duet handling:
  - round Gaussian center fit over `1:wh_peak[1]`, `1:wh_peak[2]`
  - `repeat(..., inner=ntuple(i -> i == 5 ? 2 : 1, length(r.n_dim_vars)))`

### src/persolo.jl
- Peak finding and 1D Gaussian center fit:
  - `find_peak_position_moving(...; len_avg=10)`
  - `find_positive_cluster_center(...; len_avg=10)`
  - center-fit width guess `sigma0 = clamp(n / 4, 2.0, float(n))`
  - center-fit amplitude lower bound `amp0 / 100`
  - center-fit sigma lower bound `min(2.0, float(n))`
- 2D Gaussian envelope fits:
  - `fit_dens2d_gaussian_elliptic_disk`
    - `θ_hint=(max=20.0/180*π, min=-10.0/180*π, init=10.0/180*π)`
    - `A_hint=(max=25.0, min=0, init=10.0)`
    - initial size divisor `3`
    - size bounds multipliers `/10` and `*10`
  - `fit_dens2d_gaussian_round_disk`
    - amplitude init/upper/lower `[10, 25, 0]`
    - initial size divisor `3`
    - size bounds multipliers `/10` and `*10`
- Twin-peak modulation fits:
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
- Essential/extraction processing:
  - `calc_solo_essn_2d`: default `smwh_strip=smwh`
  - modulation step formula `1 / (2 * smwh[2] * px_in_um)`
  - profile normalization `sum(prfl_modl) * step_modl / 2`
  - `calc_solo_extr`: default `proc_sidepeak=false`, `proc_envelope=false`
  - `sel_moment = y -> (y .> 0.10) .& (y .< 0.50)`
  - `sel_sidepeak = (y_modl .> 0.1) .& (y_modl .< 0.5)`

### src/corr.jl
- PCA:
  - `fit_pca_modes`: `pratio=1.0`
  - script-side `n_mode` values such as `8`
- Fourier/query weighting:
  - `query_weight`: `scaling = 1000.0`
  - Fourier phase divisor `1000.0`
  - normalization `weight / sum(weight)`
- Trend extraction:
  - `selector_t_sidepeak` and `selector_t_envelope` are externally supplied but good candidates for a trend config object
  - `freq_query` is externally supplied, currently script-side `1:1:100`

### src/percond.jl
- Stacking behavior:
  - `calc_stacked_essn` averages all numeric fields across `essns`
  - metadata-like fields are taken from `essn_ref`: `smwh`, `smwh_strip`, `smw_modl`, `step_posi`, `step_modl`
  - `offset_cent_core` and `smwh_core` are averaged with `mean_tuple`
