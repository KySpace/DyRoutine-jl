using GLMakie: draw_atomic
using Printf
using Colors: Oklch
using CairoMakie: extract_attributes!
using CairoMakie, GLMakie
using Colors: Oklch
using LaTeXStrings

# intend to layout a grid of axes with 2 main variables
# time varies in the horizontal direction, another variable vert varies in the vertical direction
# possible repetitions are all grouped together, horizontally
# dimensions have to be specified as (n_reps, n_vert, n_time, n_istp)
function set_axes_v_t_rep!(n_dim_vars::Tuple{<:Integer,<:Integer,<:Integer,<:Integer}, panel_setter::Function, runinfo, info_fmt)
    fig = Figure()
    fig[0, 1] = Label(fig, text="$(runinfo.date) run$(runinfo.runids) IB=$(@sprintf("%.3f", runinfo.IB))A $(runinfo.tag_head)"; tellwidth=false, tellheight=true, halign=:left, valign=:top)
    axs = Array{Dict}(undef, n_dim_vars[1:3]...)
    for v in 1:n_dim_vars[2], t in 1:n_dim_vars[3]
        fig[v, t][0, 1:n_dim_vars[1]] = Label(fig, text="bias=$(info_fmt[1,v,t,1].bias), t hold=$(info_fmt[1,v,t,1].t_hold) ms)"; tellwidth=true, tellheight=true, halign=:center, valign=:bottom)
        for r in 1:n_dim_vars[1]
            gl = GridLayout()
            fig[v, t][1, r] = gl
            axs[r, v, t] = panel_setter(gl)
        end
        fig[v, t].layout |> l -> colgap!(l, 0)
    end
    return fig, axs
end

function set_panel_misc_duet_2d!(gl::GridLayout)
    gl |> clean_gridlayout!
    ax_dens_1 = Axis(gl[1, 1])
    ax_dens_2 = Axis(gl[2, 1])
    ax_misc = Axis(gl[3, 1])
    colsize!(gl, 1, Fixed(100))
    colgap!(gl, 0)
    rowgap!(gl, 0)
    return Dict("dens_1" => ax_dens_1, "dens_2" => ax_dens_2, "misc" => ax_misc)
end

function to_miscibility_clr(dens1, dens2, hue1, hue2; max=16, to_norm_each=false)
    size(dens1) == size(dens2) || throw(ArgumentError("dens1 and dens2 must have the same size"))
    dens_norm_1, dens_norm_2 = (dens1, dens2) |> d -> clamp.(d, 0, max)
    shader = (a, b) -> Oklch(1 - (a + b) / 2, abs(a - b) * 0.24, a > b ? hue1 : hue2) |> RGBf
    return [shader(dens_norm_1[x, y], dens_norm_2[x, y]) for x in 1:size(dens1, 1), y in 1:size(dens1, 2)]
end

function draw_misc_duet_2d!(axs::Dict{String,Axis}, essn::AbstractVector{SoloEssentials})
    length(essn) == 2 || throw(ArgumentError("essn duet must have length 2"))
    foreach(empty!, values(axs))
    x, y = essn.smwh |> s -> map(u -> (-u:1:u), s)
    x_posi, y_posi = (x, y) .* essn.step_posi
    misc = to_miscibility_clr(essn[1].dens2d, essn[2].dens2d, hue_theme_istp[1], hue_theme_istp[2]; map=16)
    heatmap!(axs["dens_1"], x_posi, y_posi, essn[1].dens2d'; colorrange=(0, 16.0), colormap=clrmap[1])
    heatmap!(axs["dens_2"], x_posi, y_posi, essn[2].dens2d'; colorrange=(0, 10.0), colormap=clrmap[2])
    heatmap!(axs["misc"], x_posi, y_posi, misc')
    axs["dens_1"].aspect = DataAspect()
    axs["dens_2"].aspect = DataAspect()
    axs["misc"].aspect = DataAspect()
end
