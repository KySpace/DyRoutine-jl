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
