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
    colsize!(gl, 1, Fixed(100))
    colsize!(gl, 2, Fixed(100))
    colsize!(gl, 3, Fixed(200))
    rowsize!(gl, 1, Fixed(120))
    rowsize!(gl, 2, Fixed(120))
    return Dict("mode" => [ax_mode_l, ax_mode_r], "evol" => ax_evol, "spct" => ax_spct)
end

function plot_mode_evol_freq_duet!(axs::Dict{String,Axis}, mode::ModeWeight, val_t::AbstractVector)
    ndims(mode.profile) == 3 && size(mode.profile, 1) == 2 || throw(ArgumentError("mode.profile must be a 3D array with size[1]==2."))
    clrmap = gen_clrmap_posneg(0.60 * 360, 0.96 * 360)
    c = maximum(abs, mode.profile)
    heatmap!(axs["l"], mode.profile[1, :, :]'; colormap=clrmap, colorrange=(-c, c))
    heatmap!(axs["r"], mode.profile[2, :, :]'; colormap=clrmap, colorrange=(-c, c))
    axs["l"].aspect = DataAspect()
    axs["r"].aspect = DataAspect()
    axs["l"] |> hidedecorations!
    axs["r"] |> hidedecorations!
    for rep = 1:size(mode.weight, 1)
        lines!(axs["evol"], val_t, mode.weight[rep, :]; color=(:black, 0.2))
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
