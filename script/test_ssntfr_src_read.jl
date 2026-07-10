using GLMakie
using Printf
using Statistics

GLMakie.activate!()

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

# selections
ib, istp, idx_profile = (5, 1, 1)
x_row = 0.0
x_col = 0.0
ylims_profile = (-1.0, 15.0)

# after running partially from anlz_ssntfr_src.jl
function find_closest_index(x::AbstractVector{<:Real}, x_target::Real)
    isempty(x) && throw(ArgumentError("Cannot select closest index from an empty coordinate vector."))
    return argmin(abs.(x .- x_target))
end

function calc_center_span(n::Integer, smidx::Integer)
    n > 0 || throw(ArgumentError("Profile span length must be positive, got $n."))
    smidx >= 0 || throw(ArgumentError("Mean half-width must be nonnegative, got $smidx."))
    idx_center = cld(n, 2)
    return max(1, idx_center - smidx):min(n, idx_center + smidx)
end

function calc_span_edges(x::AbstractVector{<:Real}, idxs::AbstractUnitRange{<:Integer})
    length(x) >= 2 || throw(ArgumentError("Need at least two coordinate values to infer strip edges."))
    step = median(diff(x))
    return (x[first(idxs)] - step / 2, x[last(idxs)] + step / 2)
end

function disable_rectangle_zoom!(ax::Axis)
    try
        deregister_interaction!(ax, :rectanglezoom)
    catch err
        err isa KeyError || rethrow()
    end
    return ax
end

cycle_index(idx::Integer, n::Integer, step::Integer) = mod1(idx + step, n)

function calc_src_mean_image(dens_src_core::AbstractMatrix, idx_IB::Integer, idx_istp::Integer)
    dens_vec = dens_src_core[idx_IB, idx_istp]
    isempty(dens_vec) && throw(ArgumentError("dens_src_core[$idx_IB, $idx_istp] has no selected crops."))
    return mean(dens_vec)
end

function draw_src_profile_inspector!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    dens_src_core::AbstractMatrix,
    val_istp::AbstractVector;
    ib::Integer,
    istp::Integer,
    idx_profile::Integer,
    x_row::Real,
    x_col::Real,
    smidx_mean_profile::Integer,
    ylims_profile::Tuple{<:Real,<:Real}=(-1.0, 15.0),
)
    ib in axes(dens_src_core, 1) || throw(ArgumentError("ib must be in $(axes(dens_src_core, 1)), got $ib."))
    istp in axes(dens_src_core, 2) || throw(ArgumentError("istp must be in $(axes(dens_src_core, 2)), got $istp."))
    size(dens_src_core, 2) == length(val_istp) || throw(DimensionMismatch(
        "dens_src_core second dimension $(size(dens_src_core, 2)) must match length(val_istp) $(length(val_istp)).",
    ))
    for idx in CartesianIndices(dens_src_core)
        isempty(dens_src_core[idx]) && continue
        size(first(dens_src_core[idx])) == (length(x_dens), length(x_dens)) || throw(DimensionMismatch(
            "dens_src_core[$(Tuple(idx)...)] crop size $(size(first(dens_src_core[idx]))) must match " *
            "(length(x_dens), length(x_dens)) $((length(x_dens), length(x_dens))).",
        ))
    end

    dens_vec = dens_src_core[ib, istp]
    isempty(dens_vec) && throw(ArgumentError("dens_src_core[$ib, $istp] has no selected crops."))
    idx_profile = mod1(idx_profile, length(dens_vec))
    idx_row = find_closest_index(x_dens, x_row)
    idx_col = find_closest_index(x_dens, x_col)
    idxs_center = calc_center_span(length(x_dens), smidx_mean_profile)
    dens2d = dens_vec[idx_profile]
    dens_mean = calc_src_mean_image(dens_src_core, ib, istp)

    gen_theme_clr(idx_istp::Integer, alpha::Real) =
        RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
    gen_theme_clrmap(idx_istp::Integer) =
        gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)

    clr_mean = RGBAf(0.35, 0.35, 0.35, 0.62)
    clr_strip = RGBAf(0.86, 0.86, 0.86, 0.22)
    x_strip_min, x_strip_max = calc_span_edges(x_dens, idxs_center)
    y_strip_min, y_strip_max = x_strip_min, x_strip_max

    obs_idx_IB = Observable(ib)
    obs_idx_istp = Observable(istp)
    obs_idx_profile = Observable(idx_profile)
    obs_idx_row = Observable(idx_row)
    obs_idx_col = Observable(idx_col)
    obs_val_row = Observable(x_dens[idx_row])
    obs_val_col = Observable(x_dens[idx_col])
    obs_dens2d = Observable(dens2d)
    obs_dens2d_hm = lift(ds -> ds', obs_dens2d)
    obs_colorrange = Observable((0.0, maximum(dens2d)))
    obs_clrmap = Observable(gen_theme_clrmap(istp))
    obs_clr_theme = Observable(gen_theme_clr(istp, 0.92))
    obs_clr_theme_faint = Observable(gen_theme_clr(istp, 0.70))
    obs_profile_row = Observable(vec(@view dens2d[idx_row, :]))
    obs_profile_col = Observable(vec(@view dens2d[:, idx_col]))
    obs_profile_row_mean = Observable(vec(mean(@view(dens_mean[idxs_center, :]); dims=1)))
    obs_profile_col_mean = Observable(vec(mean(@view(dens_mean[:, idxs_center]); dims=2)))
    obs_title = lift(obs_idx_IB, obs_idx_istp, obs_idx_profile, obs_val_row, obs_val_col) do idx_IB_live, idx_istp_live, idx_profile_live, val_row_live, val_col_live
        @sprintf(
            "IB idx=%d, istp=%s, profile=%d/%d, x_row=%.3f μm, x_col=%.3f μm",
            idx_IB_live,
            string(val_istp[idx_istp_live]),
            idx_profile_live,
            length(dens_src_core[idx_IB_live, idx_istp_live]),
            val_row_live,
            val_col_live,
        )
    end
    obs_title_row = lift(obs_val_row) do val_row_live
        @sprintf("x_row=%.3f μm", val_row_live)
    end
    obs_title_col = lift(obs_val_col) do val_col_live
        @sprintf("x_col=%.3f μm", val_col_live)
    end

    Label(fig[0, 1:2]; text=obs_title, tellwidth=false, halign=:left)

    ax_hm = Axis(
        fig[1, 1];
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect(),
        xgridvisible=true,
        ygridvisible=true,
    )
    disable_rectangle_zoom!(ax_hm)
    hspan!(ax_hm, y_strip_min, y_strip_max; color=clr_strip)
    vspan!(ax_hm, x_strip_min, x_strip_max; color=clr_strip)
    hm = heatmap!(ax_hm, x_dens, x_dens, obs_dens2d_hm; colormap=obs_clrmap, colorrange=obs_colorrange, rasterize=true)
    hlines!(ax_hm, lift(x -> [x], obs_val_row); color=obs_clr_theme, linewidth=0.9)
    vlines!(ax_hm, lift(x -> [x], obs_val_col); color=obs_clr_theme_faint, linewidth=0.9)

    ax_row = Axis(
        fig[2, 1];
        xlabel="x (μm)",
        ylabel="density",
        title=obs_title_row,
    )
    disable_rectangle_zoom!(ax_row)
    lines!(ax_row, x_dens, obs_profile_row_mean; color=clr_mean, linewidth=2.2)
    lines!(ax_row, x_dens, obs_profile_row; color=obs_clr_theme, linewidth=1.7)
    xlims!(ax_row, extrema(x_dens))
    ylims!(ax_row, ylims_profile)

    ax_col = Axis(
        fig[1, 2];
        xlabel="density",
        ylabel="y (μm)",
        title=obs_title_col,
    )
    disable_rectangle_zoom!(ax_col)
    lines!(ax_col, obs_profile_col_mean, x_dens; color=clr_mean, linewidth=2.2)
    lines!(ax_col, obs_profile_col, x_dens; color=obs_clr_theme, linewidth=1.7)
    xlims!(ax_col, ylims_profile)
    ylims!(ax_col, extrema(x_dens))
    ax_col.xreversed = true

    function update_profiles!()
        dens_vec_live = dens_src_core[obs_idx_IB[], obs_idx_istp[]]
        isempty(dens_vec_live) && return nothing
        obs_idx_profile[] = mod1(obs_idx_profile[], length(dens_vec_live))
        dens2d_live = dens_vec_live[obs_idx_profile[]]
        dens_mean_live = calc_src_mean_image(dens_src_core, obs_idx_IB[], obs_idx_istp[])
        obs_dens2d[] = dens2d_live
        obs_colorrange[] = (0.0, maximum(dens2d_live))
        obs_clrmap[] = gen_theme_clrmap(obs_idx_istp[])
        obs_clr_theme[] = gen_theme_clr(obs_idx_istp[], 0.92)
        obs_clr_theme_faint[] = gen_theme_clr(obs_idx_istp[], 0.70)
        obs_profile_col[] = vec(@view dens2d_live[:, obs_idx_col[]])
        obs_profile_row[] = vec(@view dens2d_live[obs_idx_row[], :])
        obs_profile_row_mean[] = vec(mean(@view(dens_mean_live[idxs_center, :]); dims=1))
        obs_profile_col_mean[] = vec(mean(@view(dens_mean_live[:, idxs_center]); dims=2))
        return nothing
    end

    function update_cut_profiles!(x_click::Real, y_click::Real)
        idx_col_live = find_closest_index(x_dens, x_click)
        idx_row_live = find_closest_index(x_dens, y_click)
        obs_idx_col[] = idx_col_live
        obs_idx_row[] = idx_row_live
        obs_val_col[] = x_dens[idx_col_live]
        obs_val_row[] = x_dens[idx_row_live]
        update_profiles!()
        return nothing
    end

    function update_data_index!(step_IB::Integer, step_istp::Integer, step_profile::Integer)
        obs_idx_IB[] = cycle_index(obs_idx_IB[], size(dens_src_core, 1), step_IB)
        obs_idx_istp[] = cycle_index(obs_idx_istp[], size(dens_src_core, 2), step_istp)
        dens_vec_live = dens_src_core[obs_idx_IB[], obs_idx_istp[]]
        isempty(dens_vec_live) && return nothing
        obs_idx_profile[] = cycle_index(obs_idx_profile[], length(dens_vec_live), step_profile)
        update_profiles!()
        return nothing
    end

    click_handler = on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press && is_mouseinside(ax_hm.scene)
            xy_click = mouseposition(ax_hm)
            update_cut_profiles!(xy_click[1], xy_click[2])
        end
        return Consume(false)
    end

    gl_ctrl = GridLayout(fig[2, 2])
    labels = ("IB", "istp", "profile")
    steps = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    button_handlers = map(enumerate(labels)) do (idx_ctrl, label_ctrl)
        step = steps[idx_ctrl]
        btn_prev = Button(gl_ctrl[1, 3 * idx_ctrl - 2]; label="←", width=34, height=30)
        Label(gl_ctrl[1, 3 * idx_ctrl - 1]; text=label_ctrl, tellwidth=false, tellheight=false, halign=:center)
        btn_next = Button(gl_ctrl[1, 3 * idx_ctrl]; label="→", width=34, height=30)
        (
            on(btn_prev.clicks) do _
                update_data_index!((-step[1]), (-step[2]), (-step[3]))
            end,
            on(btn_next.clicks) do _
                update_data_index!(step...)
            end,
        )
    end

    colsize!(fig.layout, 1, Fixed(360))
    colsize!(fig.layout, 2, Fixed(300))
    rowsize!(fig.layout, 1, Fixed(360))
    rowsize!(fig.layout, 2, Fixed(260))
    resize_to_layout!(fig)
    return (;
        ax_hm,
        ax_row,
        ax_col,
        hm,
        idx_IB=obs_idx_IB,
        idx_istp=obs_idx_istp,
        idx_profile=obs_idx_profile,
        idx_row=obs_idx_row,
        idx_col=obs_idx_col,
        x_row=obs_val_row,
        x_col=obs_val_col,
        click_handler,
        button_handlers,
    )
end

@isdefined(dens_src_core) || error("Run the crop/profile part of script/anlz_ssntfr_src.jl first so dens_src_core is defined.")
@isdefined(x_dens) || error("Run the setup part of script/anlz_ssntfr_src.jl first so x_dens is defined.")
@isdefined(val_istp) || error("Run the setup part of script/anlz_ssntfr_src.jl first so val_istp is defined.")
@isdefined(smwh_src) || error("Run the setup part of script/anlz_ssntfr_src.jl first so smwh_src is defined.")

cfg = get_prfl_modl_1d_config(smwh_src)
fig_live = Figure(fontsize=14)
profile_axes = draw_src_profile_inspector!(
    fig_live,
    x_dens,
    dens_src_core,
    val_istp;
    ib,
    istp,
    idx_profile,
    x_row,
    x_col,
    smidx_mean_profile=cfg.smh_dens_strip,
    ylims_profile,
)
display(fig_live)
