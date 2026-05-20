# Configurables
runinfos
title_anlz

wh_corner = (10, 10)
smwh_roi = (30, 30)
smwh_strip = (2, 20)
wh_peak = smwh_roi .* 2 .+ 1
smw_peak, smh_peak = smwh_roi
smw_ft = 5
px_in_um = 6.5 / 22.06
step_posi = px_in_um
step_modl = 1 / (2 * smwh_roi[2] * px_in_um)
x_vec, y_vec = smwh_roi |> s -> map(u -> (-u:1:u), s)
x_posi, y_posi = (x_vec, y_vec) .* step_posi
x_modl, y_modl = (x_vec, y_vec) .* step_modl
