using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using Colors: Oklch
using LaTeXStrings

function set_axis_pca_4x4!()
    fig = Figure()
    axs_mode = Array{Dict}(undef, 16)
    for r in 1:4, c in 1:4
        gl = GridLayout()
        fig[r, c] = gl
        axs_mode[(r-1)*4+c] = set_panel_pca_duet!(gl)
    end
    return fig, axs_mode
end

function set_axis_pca_dual_4x2!()
    fig = Figure()
    axs_mode = Array{Dict}(undef, (2, 8))
    for r in 1:2, c in 1:4
        gl = GridLayout()
        fig[1, 1][r, c] = gl
        axs_mode[1, (r-1)*4+c] = set_panel_pca_solo!(gl)
        Box(fig[2, 1], color=:black)
        gl = GridLayout()
        fig[3, 1][r, c] = gl
        axs_mode[2, (r-1)*4+c] = set_panel_pca_solo!(gl)
    end
    rowsize!(fig.layout, 2, Fixed(2))
    return fig, axs_mode
end

function set_panel_pca_duet!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_l = Axis(gl[1:2, 1])
    ax_r = Axis(gl[1:2, 2])
    ax_evol = Axis(gl[1, 3])
    ax_freq = Axis(gl[2, 3])
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(100))
    colsize!(gl, 3, Fixed(200))
    rowsize!(gl, 1, Fixed(150))
    rowsize!(gl, 2, Fixed(150))
    return Dict("l" => ax_l, "r" => ax_r, "evol" => ax_evol, "freq" => ax_freq)
end

function set_panel_pca_solo!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_mode = Axis(gl[1:2, 1])
    ax_evol = Axis(gl[1, 2])
    ax_freq = Axis(gl[2, 2])
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(200))
    rowsize!(gl, 1, Fixed(120))
    rowsize!(gl, 2, Fixed(120))
    return Dict("mode" => ax_mode, "evol" => ax_evol, "freq" => ax_freq)
end

function set_panel_pca_duet!(gl::GridLayout)
    gl |> clean_gridlayout!
    gl_dens_modl = GridLayout()
    gl_evol_spct = GridLayout()
    gl[1, 1] = gl_dens_modl
    gl[1, 2] = gl_evol_spct
    ax_mode_l = Axis(gl_dens_modl[1, 1], aspect=DataAspect())
    ax_mode_r = Axis(gl_dens_modl[1, 2], aspect=DataAspect())
    ax_modl_l = Axis(gl_dens_modl[2, 1], aspect=DataAspect())
    ax_modl_r = Axis(gl_dens_modl[2, 2], aspect=DataAspect())
    ax_evol = Axis(gl[1, 2][1, 1])
    ax_spct = Axis(gl[1, 2][2, 1])
    colsize!(gl_dens_modl, 1, Fixed(200))
    colsize!(gl_dens_modl, 2, Fixed(200))
    colsize!(gl, 2, Fixed(400))
    rowsize!(gl_evol_spct, 1, Fixed(240))
    rowsize!(gl_evol_spct, 2, Fixed(240))
    rowsize!(gl_dens_modl, 1, Fixed(400))
    rowsize!(gl_dens_modl, 2, Fixed(120))
    return Dict("mode" => [ax_mode_l, ax_mode_r], "modl" => [ax_modl_l, ax_modl_r], "evol" => ax_evol, "spct" => ax_spct)
end

function gen_clrmap_posneg_nonlin(hue_pos, hue_neg; thres_alpha=0.6, alpha_base=0.2)
    return [
        begin
            alpha = abs(t) > thres_alpha ? 1.0 : (abs(t) / thres_alpha * (1 - alpha_base) + alpha_base)
            Oklch(1 - 0.6 * abs(t), 0.4 * abs2(t), t > 0 ? hue_pos : hue_neg) |> c -> RGBAf(c, alpha)
        end
        for t in range(-1, 1; length=256)
    ]
end

function plot_mode_evol_spct_duet!(axs::Dict{String}, mode::ModeWeight, spectral::NamedTuple, rcrd_pks, val_istp; step_posi::Real=1, smwh=(0, 0))
    step_modl = 1 ./ (2 .* smwh .* step_posi)
    x_vec, y_vec = smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x_vec, y_vec) .* step_posi
    x_modl, y_modl = (x_vec, y_vec) .* step_modl
    n_rep = spectral.n_rep
    val_t = spectral.val_t
    freq_query = spectral.freq_query
    mask_evol = spectral.mask_evol
    evols_weight = spectral.evols
    evol_weight_mean = spectral.evol_mean
    spectra_reps_mask = spectral.spectra_reps_mask
    spct_mean_full = spectral.spct_mean_full
    spct_mean_mask = spectral.spct_mean_mask
    axs |> clear_axes!
    axs["mode"] |> clear_axes!
    length(mode.profile) == 2 || throw(ArgumentError("mode.profile must have 2 components."))
    clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360)
    clr_grid = RGBAf(Oklch(0.84, 0.0, 262), 1)
    c = maximum(abs, mode.profile |> stack)
    step_t = length(val_t) > 1 ? minimum(diff(val_t)) : one(eltype(val_t))
    t_span_lim = any(mask_evol) ? val_t[mask_evol] |> t -> (minimum(t) - step_t / 2, maximum(t) + step_t / 2) : nothing
    for i in 1:2
        clrmap_modl = gen_clrmap_solo(hue_theme_istp[val_istp[i]]; thres_alpha=0.6, alpha_base=0.2)
        ax_dens = axs["mode"][i]
        hm = heatmap!(ax_dens, x_posi, y_posi, mode.profile[i]'; colormap=clrmap, colorrange=(-c, c))
        # translate!(hm, 0, 0, -100)
        ax_dens.aspect = DataAspect()
        ax_dens |> ax -> hidedecorations!(ax, ticks=false, label=true, grid=false, minorgrid=false)
        ax_dens.yticks = -10:5:10
        ax_dens.yminorticks = IntervalsBetween(5)
        ax_dens.yminorgridvisible = true
        ax_dens.ygridcolor = clr_grid
        ax_dens.yminorgridcolor = clr_grid
        ax_modl = axs["modl"][i]
        modl = mode.profile[i] |> d -> d .* gen_win_hann_2d(smwh) |> fft |> fftshift |> c -> abs2.(c)
        hm = heatmap!(ax_modl, x_modl, y_modl, modl'; colormap=clrmap_modl)
        ylims!(ax_modl, (0, 0.6))
        xlims!(ax_modl, (-0.5, 0.5))
        ax_modl |> ax -> hidedecorations!(ax, ticks=false, label=false, grid=false)
        ax_modl.ygridvisible = true
        ax_modl.xgridvisible = true
        ax_modl.yticks = 0:0.1:0.6
        ax_modl.ygridcolor = clr_grid
        if i == 1
            ax_dens.yticklabelsvisible = true
            ax_modl.yticklabelsvisible = true
            ax_dens.ylabelvisible = true
            ax_modl.ylabelvisible = true
            ax_dens.ylabel = "x (μm)"
            ax_modl.ylabel = rich("k", subscript("y"), " (μm⁻¹)")
        end
        ax_modl.xlabel = val_istp[i]
    end
    for pk in rcrd_pks
        scatter!(axs["spct"], pk.freq, pk.value; color=(:darkorchid4, 1))
        str_val_rel = @sprintf("%.2f", pk.value_reduced) |> s -> replace(s, r"^0" => "")
        text!(axs["spct"], pk.freq, pk.value; text=@sprintf("%.0f Hz \n%s", pk.freq, str_val_rel), color=(:darkorchid4, 1), fontsize=14, align=(:left, :bottom))
    end
    isnothing(t_span_lim) || vspan!(axs["evol"], t_span_lim...; color=RGBAf(Oklch(0.4, 0.01, 240), 0.04))
    for rep = 1:n_rep
        clr = Oklch(0.86, 0.053, mod(rep / 6 - 0.1, 1) * 360) |> c -> RGBAf(c, 1)
        scatter!(axs["evol"], val_t, evols_weight[rep]; color=clr)
        lines!(axs["spct"], freq_query, spectra_reps_mask[rep]; color=clr)
    end
    lines!(axs["evol"], val_t, evol_weight_mean; color=(:black, 1))
    lines!(axs["spct"], freq_query, spct_mean_full; color=(:black, 0.5), linestyle=:dash)
    lines!(axs["spct"], freq_query, spct_mean_mask; color=(:black, 1.0))
    for ax in [axs["evol"], axs["spct"]]
        ax.yticklabelspace = 40.0
    end
    axs["evol"].xlabel = "time (ms)"
    axs["spct"].xlabel = "frequency (Hz)"
    ylims!(axs["spct"], (-0.05, 1.15))
end

function plot_mode_evol_freq_solo!(axs::Dict{String,Axis}, mode::ModeWeight, val_t::AbstractVector)
    ndims(mode.profile) == 2 || throw(ArgumentError("mode.profile must be a 2D array. "))
    clrmap = gen_clrmap_posneg(0.60 * 360, 0.96 * 360)
    c = maximum(abs, mode.profile)
    heatmap!(axs["mode"], mode.profile[:, :]; colormap=clrmap, colorrange=(-c, c))
    axs["mode"].aspect = DataAspect()
    axs["mode"] |> hidedecorations!
    for rep = 1:size(mode.weight, 1)
        lines!(axs["evol"], val_t, mode.weight[rep, :]; color=(:black, 0.2))
    end
end
