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

function set_panel_pca_solo!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_mode = Axis(gl[1:2, 1])
    ax_evol = Axis(gl[1, 2])
    ax_spct = Axis(gl[2, 2])
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(200))
    rowsize!(gl, 1, Fixed(120))
    rowsize!(gl, 2, Fixed(120))
    return Dict("mode" => ax_mode, "evol" => ax_evol, "spct" => ax_freq)
end

function set_panel_pca_duet!(gl::GridLayout)
    gl |> clean_gridlayout!
    gl_dens_modl = GridLayout()
    gl_evol_spct = GridLayout()
    gl[1, 1] = gl_dens_modl
    gl[1, 2] = gl_evol_spct
    axs_mode = set_pca_mode_axes_2d_duet!(gl_dens_modl)
    ax_evol = Axis(gl[1, 2][1, 1])
    ax_spct = Axis(gl[1, 2][2, 1])
    colsize!(gl, 2, Fixed(400))
    rowsize!(gl_evol_spct, 1, Fixed(240))
    rowsize!(gl_evol_spct, 2, Fixed(240))
    return merge(axs_mode, Dict("evol" => ax_evol, "spct" => ax_spct))
end

function set_panel_pca_duet_params!(
    gl::GridLayout,
    val_params::AbstractVector;
    mode_kind::Symbol=:profile1d,
    width_evol::Real=400,
    width_spct::Real=400,
    height_evol::Real=160,
    height_spct::Real=160,
)
    isempty(val_params) && throw(ArgumentError("val_params must not be empty"))
    gl |> clean_gridlayout!
    gl_mode = GridLayout()
    gl_rows = GridLayout()
    gl[1, 1] = gl_mode
    gl[2, 1] = gl_rows
    axs_mode = mode_kind == :profile1d ? set_pca_mode_axes_profile_duet!(gl_mode) :
               mode_kind == :dens2d ? set_pca_mode_axes_2d_duet!(gl_mode) :
               throw(ArgumentError("unknown mode_kind $mode_kind"))

    axs_evol = Vector{Axis}(undef, length(val_params))
    axs_spct = Vector{Axis}(undef, length(val_params))
    for (idx_param, param) in enumerate(val_params)
        Label(gl_rows[idx_param, 0]; text="$(@sprintf("%.3f", param))", rotation=pi / 2, tellwidth=true, tellheight=false)
        axs_evol[idx_param] = Axis(gl_rows[idx_param, 1]; width=width_evol, height=height_evol)
        axs_spct[idx_param] = Axis(gl_rows[idx_param, 2]; width=width_spct, height=height_spct)
    end
    colsize!(gl_rows, 0, Fixed(24))
    colgap!(gl_rows, 8)
    rowgap!(gl_rows, 4)
    return merge(axs_mode, Dict("evol" => axs_evol, "spct" => axs_spct))
end

function set_pca_mode_axes_2d_duet!(gl::GridLayout)
    ax_mode_l = Axis(gl[1, 1], aspect=DataAspect())
    ax_mode_r = Axis(gl[1, 2], aspect=DataAspect())
    ax_modl_l = Axis(gl[2, 1], aspect=DataAspect())
    ax_modl_r = Axis(gl[2, 2], aspect=DataAspect())
    colsize!(gl, 1, Fixed(200))
    colsize!(gl, 2, Fixed(200))
    rowsize!(gl, 1, Fixed(400))
    rowsize!(gl, 2, Fixed(120))
    return Dict("mode" => [ax_mode_l, ax_mode_r], "modl" => [ax_modl_l, ax_modl_r])
end

function set_pca_mode_axes_profile_duet!(gl::GridLayout)
    ax_mode_l = Axis(gl[1, 1])
    ax_mode_r = Axis(gl[2, 1])
    colsize!(gl, 1, Fixed(400))
    rowsize!(gl, 1, Fixed(72))
    rowsize!(gl, 2, Fixed(72))
    rowgap!(gl, 4)
    return Dict("mode" => [ax_mode_l, ax_mode_r])
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

function draw_pca_mode_2d_duet!(axs::Dict, mode::ModeWeight, val_istp; step_posi::Real=1, smwh=(0, 0))
    step_modl = 1 ./ (2 .* smwh .* step_posi)
    x_vec, y_vec = smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x_vec, y_vec) .* step_posi
    x_modl, y_modl = (x_vec, y_vec) .* step_modl
    axs["mode"] |> clear_axes!
    length(mode.profile) == 2 || throw(ArgumentError("mode.profile must have 2 components."))
    clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360)
    clr_grid = RGBAf(Oklch(0.84, 0.0, 262), 1)
    c = maximum(abs, mode.profile |> stack)
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
    return nothing
end

function draw_pca_mode_profile_duet!(
    axs::Dict,
    mode::ModeWeight,
    val_istp,
    y_modl::AbstractVector;
    x_lims=(0.06, 0.6),
)
    axs["mode"] |> clear_axes!
    ndims(mode.profile) == 2 || throw(ArgumentError("profile PCA mode must be a 2D matrix of istp × wavenum."))
    size(mode.profile, 1) == length(val_istp) ||
        throw(DimensionMismatch("profile mode istp count $(size(mode.profile, 1)) does not match val_istp length $(length(val_istp))"))
    size(mode.profile, 2) == length(y_modl) ||
        throw(DimensionMismatch("profile mode wavenum count $(size(mode.profile, 2)) does not match y_modl length $(length(y_modl))"))
    c = maximum(abs, mode.profile)
    clrmap = gen_clrmap_posneg_nonlin(0.57 * 360, 0.96 * 360)
    for i in axes(mode.profile, 1)
        ax = axs["mode"][i]
        heatmap!(ax, y_modl, [0.0], reshape(mode.profile[i, :], :, 1); colormap=clrmap, colorrange=(-c, c))
        xlims!(ax, x_lims)
        ylims!(ax, (-0.5, 0.5))
        ax.yticks = ([0.0], [string(val_istp[i])])
        ax.yticklabelspace = 32.0
        ax.xticklabelspace = 28.0
        ax.xlabel = i == length(val_istp) ? "wavenum (μm⁻¹)" : ""
        ax.xgridvisible = true
        ax.ygridvisible = false
        ax.leftspinevisible = false
        ax.rightspinevisible = false
    end
    return nothing
end

function plot_pca_evol_spct!(ax_evol::Axis, ax_spct::Axis, spectral::NamedTuple, rcrd_pks)
    n_rep = spectral.n_rep
    val_t = spectral.val_t
    freq_query = spectral.freq_query
    mask_evol = spectral.mask_evol
    evols_weight = spectral.evols
    evol_weight_mean = spectral.evol_mean
    spectra_reps_mask = spectral.spectra_reps_mask
    spct_mean_full = spectral.spct_mean_full
    spct_mean_mask = spectral.spct_mean_mask
    step_t = length(val_t) > 1 ? minimum(diff(val_t)) : one(eltype(val_t))
    t_span_lim = any(mask_evol) ? val_t[mask_evol] |> t -> (minimum(t) - step_t / 2, maximum(t) + step_t / 2) : nothing
    for pk in rcrd_pks
        scatter!(ax_spct, pk.freq, pk.value; color=(:darkorchid4, 1))
        str_val_rel = @sprintf("%.2f", pk.value_reduced) |> s -> replace(s, r"^0" => "")
        text!(ax_spct, pk.freq, pk.value; text=@sprintf("%.0f Hz \n%s", pk.freq, str_val_rel), color=(:darkorchid4, 1), fontsize=14, align=(:left, :bottom))
    end
    isnothing(t_span_lim) || vspan!(ax_evol, t_span_lim...; color=RGBAf(Oklch(0.4, 0.01, 240), 0.04))
    for rep = 1:n_rep
        clr = Oklch(0.86, 0.053, mod(rep / 6 - 0.1, 1) * 360) |> c -> RGBAf(c, 1)
        scatter!(ax_evol, val_t, evols_weight[rep]; color=clr)
        lines!(ax_spct, freq_query, spectra_reps_mask[rep]; color=clr)
    end
    lines!(ax_evol, val_t, evol_weight_mean; color=(:black, 1))
    lines!(ax_spct, freq_query, spct_mean_full; color=(:black, 0.5), linestyle=:dash)
    lines!(ax_spct, freq_query, spct_mean_mask; color=(:black, 1.0))
    for ax in [ax_evol, ax_spct]
        ax.yticklabelspace = 44.0
        ax.xticklabelspace = 30.0
    end
    ax_evol.xlabel = "time (ms)"
    ax_spct.xlabel = "frequency (Hz)"
    ylims!(ax_spct, (-0.05, 1.15))
    return nothing
end

function plot_mode_evol_spct_duet!(axs::Dict, mode::ModeWeight, spectral::NamedTuple, rcrd_pks, val_istp; step_posi::Real=1, smwh=(0, 0))
    axs |> clear_axes!
    draw_pca_mode_2d_duet!(axs, mode, val_istp; step_posi, smwh)
    plot_pca_evol_spct!(axs["evol"], axs["spct"], spectral, rcrd_pks)
    return nothing
end

function plot_mode_evol_spct_duet_params!(
    axs::Dict,
    mode::ModeWeight,
    spectra_params::AbstractVector,
    val_params::AbstractVector,
    val_istp;
    mode_kind::Symbol=:profile1d,
    y_modl=nothing,
    step_posi::Real=1,
    smwh=(0, 0),
)
    length(spectra_params) == length(val_params) ||
        throw(DimensionMismatch("spectra_params length $(length(spectra_params)) does not match val_params length $(length(val_params))"))
    axs |> clear_axes!
    if mode_kind == :profile1d
        isnothing(y_modl) && throw(ArgumentError("y_modl is required for profile1d PCA mode plotting."))
        draw_pca_mode_profile_duet!(axs, mode, val_istp, y_modl)
    elseif mode_kind == :dens2d
        draw_pca_mode_2d_duet!(axs, mode, val_istp; step_posi, smwh)
    else
        throw(ArgumentError("unknown mode_kind $mode_kind"))
    end
    for idx_param in eachindex(val_params)
        spectral, peaks = spectra_params[idx_param]
        plot_pca_evol_spct!(axs["evol"][idx_param], axs["spct"][idx_param], spectral, peaks)
    end
    return nothing
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
