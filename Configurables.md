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

fit_prfl_modl_twinpeak_decay_1d
fit_prfl_modl_twinpeak_1d
    prfl hint:
          M,  σ0,   P,    σ,   p,   D,   λ
        3.0, 0.1, 0.5, 0.05, 0.3, 0.5, 0.8
        
- corr.jl
  - query_weight
    - scaling

- persolo.jl
  - fit_prfl_modl_twinpeak_decay_1d