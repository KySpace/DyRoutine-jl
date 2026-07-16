using ImageMorphology
using ImageSegmentation
ib, istp = (5, 2)
# ft_eg = abs.(dens_core_ft_cmpx_mean[ib, istp])
clrmap = gen_clrmap_solo(hue_theme_istp[string(val_istp[istp])]; alpha_base=0.2, thres_alpha=0.1)
# 
# fig = Figure()
# axs_2d = Axis(fig[1,1]; aspect=DataAspect(), width=200, height=150)
# axs_1d = Axis(fig[2,1]; width=200, height=150)
# heatmap!(axs_2d, kx_ft, ky_ft, ft_eg'; colorrange=clrrng_ft2d_cmpx_mean, colormap=clrmap)
# lines!(axs_1d, kx_ft, vec(mean(ft_eg; dims=1)))
# ylims!(axs_1d, (0, 10))
# linkxaxes!(axs_1d, axs_2d)
# fig |> display

ft_eg = abs.(ft2d_absl_mean[ib, istp])
ft_eg_x = abs.(ft2d_cmpx_mean[ib, istp])


mask_sidepeak = begin 
    arg_seed_main = argmin([hypot(y, x .- 0) for y in ky_ft, x in kx_ft])
    arg_seed_side = argmin([hypot(y, x .- 0.2) for y in ky_ft, x in kx_ft])
    ft_this = ft2d_absl_mean[ib, istp]
    markers = zeros(Int, size(ft_this)); markers[arg_seed_main] = 1; markers[arg_seed_side] = 2
    seg = watershed(.-ft_this, markers)
    labels_map(seg) .== 2
end

function to_masked_clr(dens, mask, hue; sat_max=0.24, max=16, thres_alpha=0.1, l_max=1.0, l_min=0.0, alpha_base=0.1)
    size(dens) == size(mask) || throw(DimensionMismatch("dens size $(size(dens)) does not match mask size $(size(mask))."))
    dens_norm = clamp.(dens, 0, max) ./ max
    alpha = (n, m) -> m ? (thres_alpha <= 0 ? (n > 0 ? 1.0 : alpha_base) : (n > thres_alpha ? 1.0 : (n / thres_alpha * (1 - alpha_base) + alpha_base))) : 0.0
    shader = (n, m) -> Oklch(l_max - (l_max - l_min) * abs(n), sat_max * abs(n), hue) |> c -> RGBAf(c, alpha(n, m))
    return [shader(dens_norm[x, y], mask[x, y]) for x in 1:size(dens, 1), y in 1:size(dens, 2)]
end

fig = Figure()
axs_2d = Axis(fig[1,1]; aspect=DataAspect(), width=400, height=300)
axs_mask = Axis(fig[2,1]; aspect=DataAspect(), width=400, height=300)
cohr_masked = to_masked_clr(ft_eg_x, mask_sidepeak, hue_theme_istp[string(val_istp[istp])], max=30)
cohr_nonmasked = to_masked_clr(ft_eg_x, .!mask_sidepeak, 0; sat_max=0.0, max=30)

inco_masked = to_masked_clr(ft_eg, mask_sidepeak, hue_theme_istp[string(val_istp[istp])], max=80)
inco_nonmasked = to_masked_clr(ft_eg, .!mask_sidepeak, 0; sat_max=0.0, max=80)
heatmap!(axs_2d, kx_ft, ky_ft, inco_nonmasked')
heatmap!(axs_2d, kx_ft, ky_ft, inco_masked')
heatmap!(axs_mask, kx_ft, ky_ft, cohr_nonmasked')
heatmap!(axs_mask, kx_ft, ky_ft, cohr_masked')
limits!(axs_2d, (0, 0.5), (-0.2, 0.2))
limits!(axs_mask, (0, 0.5), (-0.2, 0.2))

fig |> resize_to_layout!
fig |> display
(ft_eg_x[mask_sidepeak] |> sum) / (ft_eg_x |> sum)
