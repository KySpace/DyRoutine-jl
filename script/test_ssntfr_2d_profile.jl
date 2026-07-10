using GLMakie
using HDF5
using Printf
using Statistics
GLMakie.activate!()
include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
include(joinpath(@__DIR__, "..", "src", "persolo.jl"))
include(joinpath(@__DIR__, "..", "src", "modlntfr.jl"))

# selections
ib, istp = (5, 1)
x_row = 0.0
x_col = 0.0
smidx_mean_profile = 20
ylims_profile = (-1.0, 15.0)

# after running partially from anlz_ssntfr_2d.jl
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

function draw_ntfr2d_profile_inspector!(
    fig::Figure,
    x_dens::AbstractVector{<:Real},
    ntfr2d_fmt::AbstractMatrix{<:AbstractMatrix},
    val_istp::AbstractVector;
    ib::Integer,
    istp::Integer,
    x_row::Real,
    x_col::Real,
    smidx_mean_profile::Integer=20,
    ylims_profile::Tuple{<:Real,<:Real}=(-1.0, 15.0),
)
    size(ntfr2d_fmt, 2) == length(val_istp) || throw(DimensionMismatch(
        "ntfr2d_fmt second dimension $(size(ntfr2d_fmt, 2)) must match length(val_istp) $(length(val_istp)).",
    ))
    ib in axes(ntfr2d_fmt, 1) || throw(ArgumentError("ib must be in $(axes(ntfr2d_fmt, 1)), got $ib."))
    istp in axes(ntfr2d_fmt, 2) || throw(ArgumentError("istp must be in $(axes(ntfr2d_fmt, 2)), got $istp."))
    for idx in CartesianIndices(ntfr2d_fmt)
        size(ntfr2d_fmt[idx]) == (length(x_dens), length(x_dens)) || throw(DimensionMismatch(
            "ntfr2d_fmt[$(Tuple(idx)...)] size $(size(ntfr2d_fmt[idx])) must match " *
            "(length(x_dens), length(x_dens)) $((length(x_dens), length(x_dens))).",
        ))
    end

    idx_row = find_closest_index(x_dens, x_row)
    idx_col = find_closest_index(x_dens, x_col)
    idxs_center = calc_center_span(length(x_dens), smidx_mean_profile)
    val_row = x_dens[idx_row]
    val_col = x_dens[idx_col]
    dens2d = ntfr2d_fmt[ib, istp]
    istp_label = val_istp[istp]
    gen_theme_clr(idx_istp::Integer, alpha::Real) =
        RGBAf(Oklch(0.52, 0.14, hue_theme_istp[string(val_istp[idx_istp])]), alpha)
    gen_theme_clrmap(idx_istp::Integer) =
        gen_clrmap_solo(hue_theme_istp[string(val_istp[idx_istp])]; alpha_base=0.2, thres_alpha=0.1)

    clr_mean = RGBAf(0.35, 0.35, 0.35, 0.62)
    clr_strip = RGBAf(0.86, 0.86, 0.86, 0.22)

    profile_row = vec(@view dens2d[idx_row, :])
    profile_col = vec(@view dens2d[:, idx_col])
    profile_row_mean = vec(mean(@view(dens2d[idxs_center, :]); dims=1))
    profile_col_mean = vec(mean(@view(dens2d[:, idxs_center]); dims=2))
    x_strip_min, x_strip_max = calc_span_edges(x_dens, idxs_center)
    y_strip_min, y_strip_max = x_strip_min, x_strip_max
    obs_idx_IB = Observable(ib)
    obs_idx_istp = Observable(istp)
    obs_idx_row = Observable(idx_row)
    obs_idx_col = Observable(idx_col)
    obs_dens2d = Observable(dens2d)
    obs_dens2d_hm = lift(ds -> ds', obs_dens2d)
    obs_colorrange = Observable((0.0, maximum(dens2d)))
    obs_clrmap = Observable(gen_theme_clrmap(istp))
    obs_clr_theme = Observable(gen_theme_clr(istp, 0.92))
    obs_clr_theme_faint = Observable(gen_theme_clr(istp, 0.70))
    obs_val_row = Observable(val_row)
    obs_val_col = Observable(val_col)
    obs_profile_row = Observable(profile_row)
    obs_profile_col = Observable(profile_col)
    obs_profile_row_mean = Observable(profile_row_mean)
    obs_profile_col_mean = Observable(profile_col_mean)
    obs_title = lift(obs_idx_IB, obs_idx_istp, obs_val_row, obs_val_col) do idx_IB_live, idx_istp_live, val_row_live, val_col_live
        @sprintf(
            "IB idx=%d, istp=%s, x_row = %.3f μm, x_col = %.3f μm",
            idx_IB_live,
            string(val_istp[idx_istp_live]),
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

    Label(
        fig[0, 1:2];
        text=obs_title,
        tellwidth=false,
        halign=:left,
    )

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
    ax_col.xreversed = true;

    function update_profiles!()
        dens2d_live = ntfr2d_fmt[obs_idx_IB[], obs_idx_istp[]]
        obs_dens2d[] = dens2d_live
        obs_colorrange[] = (0.0, maximum(dens2d_live))
        obs_clrmap[] = gen_theme_clrmap(obs_idx_istp[])
        obs_clr_theme[] = gen_theme_clr(obs_idx_istp[], 0.92)
        obs_clr_theme_faint[] = gen_theme_clr(obs_idx_istp[], 0.70)
        obs_profile_col[] = vec(@view dens2d_live[:, obs_idx_col[]])
        obs_profile_row[] = vec(@view dens2d_live[obs_idx_row[], :])
        obs_profile_row_mean[] = vec(mean(@view(dens2d_live[idxs_center, :]); dims=1))
        obs_profile_col_mean[] = vec(mean(@view(dens2d_live[:, idxs_center]); dims=2))
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

    cycle_index(idx::Integer, n::Integer, step::Integer) = mod1(idx + step, n)

    function update_data_index!(step_IB::Integer, step_istp::Integer)
        obs_idx_IB[] = cycle_index(obs_idx_IB[], size(ntfr2d_fmt, 1), step_IB)
        obs_idx_istp[] = cycle_index(obs_idx_istp[], size(ntfr2d_fmt, 2), step_istp)
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
    btn_IB_up = Button(gl_ctrl[1, 2]; label="↑", width=44, height=34)
    btn_istp_left = Button(gl_ctrl[2, 1]; label="←", width=44, height=34)
    btn_istp_right = Button(gl_ctrl[2, 3]; label="→", width=44, height=34)
    btn_IB_down = Button(gl_ctrl[3, 2]; label="↓", width=44, height=34)
    Label(gl_ctrl[2, 2]; text="IB\nistp", tellwidth=false, tellheight=false, halign=:center)
    button_handlers = (
        on(btn_IB_up.clicks) do _
            update_data_index!(1, 0)
        end,
        on(btn_IB_down.clicks) do _
            update_data_index!(-1, 0)
        end,
        on(btn_istp_left.clicks) do _
            update_data_index!(0, -1)
        end,
        on(btn_istp_right.clicks) do _
            update_data_index!(0, 1)
        end,
    )
    colsize!(fig.layout, 1, Fixed(360))
    colsize!(fig.layout, 2, Fixed(260))
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
        idx_row=obs_idx_row,
        idx_col=obs_idx_col,
        x_row=obs_val_row,
        x_col=obs_val_col,
        click_handler,
        button_handlers,
    )
end

@isdefined(ntfr2d_mean) || error("Run the data-loading part of script/anlz_ssntfr_2d.jl first so ntfr2d_mean is defined.")
@isdefined(x_dens) || error("Run the data-loading part of script/anlz_ssntfr_2d.jl first so x_dens is defined.")
@isdefined(val_istp) || error("Run the setup part of script/anlz_ssntfr_2d.jl first so val_istp is defined.")

fig_live = Figure(fontsize=14)
profile_axes = draw_ntfr2d_profile_inspector!(
    fig_live,
    x_dens,
    ntfr2d_mean,
    val_istp;
    ib,
    istp,
    x_row,
    x_col,
    smidx_mean_profile,
    ylims_profile,
)
display(fig_live)
