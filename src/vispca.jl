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
    ax_mode_l = Axis(gl[1:2, 1])
    ax_mode_r = Axis(gl[1:2, 2])
    ax_evol = Axis(gl[1, 3])
    ax_spct = Axis(gl[2, 3])
    colsize!(gl, 1, Fixed(200))
    colsize!(gl, 2, Fixed(200))
    colsize!(gl, 3, Fixed(400))
    rowsize!(gl, 1, Fixed(240))
    rowsize!(gl, 2, Fixed(240))
    return Dict("mode" => [ax_mode_l, ax_mode_r], "evol" => ax_evol, "spct" => ax_spct)
end

function gen_clrmap_posneg_nonlin(hue_pos, hue_neg)
    thres_alpha = 0.6
    alpha_base = 0.2
    return [
        begin
            alpha = abs(t) > thres_alpha ? 1.0 : (abs(t) / thres_alpha * (1 - alpha_base) + alpha_base)
            Oklch(1 - 0.6 * abs(t), 0.4 * abs2(t), t > 0 ? hue_pos : hue_neg) |> c -> RGBAf(c, alpha)
        end
        for t in range(-1, 1; length=256)
    ]
end

function plot_mode_evol_spct_duet!(axs::Dict{String}, mode::ModeWeight, val_t::AbstractVector, freq_query::AbstractVector, sel_evo::Function; step_posi::Real=1, smwh=(0, 0))
    x_vec, y_vec = smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x_vec, y_vec) .* step_posi
    axs |> clear_axes!
    axs["mode"] |> clear_axes!
    length(mode.profile) == 2 || throw(ArgumentError("mode.profile must have 2 components."))
    clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360)
    clr_grid = RGBAf(Oklch(0.84, 0.0, 262), 1)
    c = maximum(abs, mode.profile |> stack)
    mask_evo = map(sel_evo, val_t)
    step_t = val_t |> diff |> minimum
    t_span_lim = val_t[mask_evo] |> t -> (minimum(t) - step_t / 2, maximum(t) + step_t / 2)
    for i in 1:2
        ax = axs["mode"][i]
        hm = heatmap!(ax, x_posi, y_posi, mode.profile[i]'; colormap=clrmap, colorrange=(-c, c))
        # translate!(hm, 0, 0, -100)
        ax.aspect = DataAspect()
        ax |> ax -> hidedecorations!(ax, ticks=false, label=true, grid=false, minorgrid=false)
        ax.yticks = -10:5:10
        ax.yminorticks = IntervalsBetween(5)
        ax.yminorgridvisible = true
        ax.ygridcolor = clr_grid
        ax.yminorgridcolor = clr_grid
    end
    n_rep = size(mode.weight, 1)
    spectra = [
        mode.weight[r, :] |> evo -> query_weight(evo, mask_evo, val_t, freq_query)
        for r in 1:n_rep
    ]
    weight_mean_evol = mean(mode.weight, dims=1) |> vec
    weight_mean_spec_full = weight_mean_evol |> evo -> query_weight(evo, :, val_t, freq_query)
    weight_mean_spec_mask = weight_mean_evol |> evo -> query_weight(evo, mask_evo, val_t, freq_query)
    vspan!(axs["evol"], t_span_lim...; color=RGBAf(Oklch(0.4, 0.01, 240), 0.04))
    for rep = 1:n_rep
        clr = Oklch(0.86, 0.053, mod(rep / 6 - 0.1, 1) * 360) |> c -> RGBAf(c, 1)
        scatter!(axs["evol"], val_t, mode.weight[rep, :]; color=clr)
        lines!(axs["spct"], freq_query, spectra[rep]; color=clr)
    end
    lines!(axs["evol"], val_t, weight_mean_evol; color=(:black, 1))
    lines!(axs["spct"], freq_query, weight_mean_spec_full; color=(:black, 0.5))
    lines!(axs["spct"], freq_query, weight_mean_spec_mask; color=(:black, 1.0))
    for ax in [axs["evol"], axs["spct"]]
        ax.yticklabelspace=40.0
    end
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
