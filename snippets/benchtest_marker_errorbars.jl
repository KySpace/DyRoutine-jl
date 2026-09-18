using GLMakie

include(joinpath(@__DIR__, "marker_errorbars.jl"))

GLMakie.activate!()

x = collect(1.0:4.0)
markers = [:circle, :rect, :utriangle, :diamond]
markersizes = [2, 12, 30, 50]
markercolors = Makie.wong_colors()[[2,2,2,2]]
stroke_colors = [:black, :darkred, :darkgreen, :navy]
tick_labels = ["circle\n2", "rect\n12", "utriangle\n30", "diamond\n50"]

linear_y = zeros(4)
linear_error_low = [0.005, 0.018, 0.025, 0.045]
linear_error_high = [0.008, 0.018, 0.035, 0.060]
linear_limits = [(-0.8, 0.8), (-0.45, 0.45), (-0.10, 0.10)]

log_y = fill(1_000.0, 4)
log_error_low = [50.0, 120.0, 250.0, 400.0]
log_error_high = [80.0, 180.0, 350.0, 600.0]
log_limits = [(100.0, 10_000.0), (300.0, 3_000.0), (500.0, 2_000.0)]

fig = Figure(size=(1080, 1160))
linear_axes = []
log_axes = []
linear_results = []
log_results = []

for column in eachindex(linear_limits)
    ax_linear = Axis(
        fig[1, column];
        title=("wide", "medium", "tight")[column] * " y range",
        ylabel=column == 1 ? "linear y" : "",
        xticks=(x, tick_labels),
    )
    push!(linear_axes, ax_linear)
    xlims!(ax_linear, 0.45, 4.55)
    ylims!(ax_linear, linear_limits[column]...)
    push!(linear_results, marker_errorbars!(
        ax_linear,
        x,
        linear_y,
        linear_error_low,
        linear_error_high;
        marker=markers,
        markersize=markersizes,
        color=markercolors,
        strokecolor=stroke_colors,
        strokewidth=2,
        errorlinewidth=2,
    ))
    hlines!(ax_linear, 0; color=(:gray, 0.25), linewidth=1)

    ax_log = Axis(
        fig[2, column];
        yscale=log10,
        ylabel=column == 1 ? "logarithmic y" : "",
        xticks=(x, tick_labels),
    )
    push!(log_axes, ax_log)
    xlims!(ax_log, 0.45, 4.55)
    ylims!(ax_log, log_limits[column]...)
    push!(log_results, marker_errorbars!(
        ax_log,
        x,
        log_y,
        log_error_low,
        log_error_high;
        marker=markers,
        markersize=markersizes,
        color=markercolors,
        strokecolor=stroke_colors,
        strokewidth=2,
        errorlinewidth=2,
    ))
end

Label(
    fig[0, :],
    "Marker-edge error bars on linear and logarithmic axes";
    fontsize=22,
    font=:bold,
)
Label(
    fig[4, :],
    "Every segment joins the glyph edge to the fixed data endpoint y ± error.";
    fontsize=14,
)

ax_overlap = Axis(
    fig[3, :];
    title="overlap order: each later marker covers earlier error bars",
    ylabel="linear y",
    xticksvisible=false,
    xticklabelsvisible=false,
)
xlims!(ax_overlap, 0.7, 1.3)
ylims!(ax_overlap, -0.13, 0.13)
x_overlap = [1.000, 1.006, 1.012]
y_overlap = zeros(3)
error_overlap = [0.10, 0.075, 0.05]
overlap_result = marker_errorbars!(
    ax_overlap,
    x_overlap,
    y_overlap,
    error_overlap;
    marker=[:diamond, :rect, :circle],
    markersize=[50, 36, 24],
    color=fill(markercolors[1], 3),
    strokecolor=[:navy, :darkred, :black],
    strokewidth=2,
    errorlinewidth=2,
)

resize_to_layout!(fig)

function assert_error_endpoints(result, x, y, error_low, error_high)
    @assert length(result.segments[]) == 4 * length(x)
    for i in eachindex(x)
        j = 4i - 3
        points = result.segments[]
        @assert all(point -> point[1] ≈ x[i], points[j:(j + 3)])
        @assert points[j + 1][2] ≈ y[i] - error_low[i]
        @assert points[j + 3][2] ≈ y[i] + error_high[i]
    end
end

function assert_marker_edges(result, ax, y)
    limits = ax.finallimits[]
    viewport = ax.scene.viewport[]
    yscale = ax.yscale[]
    y_low = limits.origin[2]
    y_high = y_low + widths(limits)[2]
    pixels_per_scaled = widths(viewport)[2] / (yscale(y_high) - yscale(y_low))
    markers = result.scatter.marker[]
    markersizes = result.scatter.markersize[]
    strokewidths = result.scatter.strokewidth[]

    for i in eachindex(y)
        j = 4i - 3
        marker_i = _marker_value(markers, i)
        markersize_i = _marker_size_value(markersizes, i)
        strokewidth_i = _marker_value(strokewidths, i)
        expected_bottom, expected_top = _visible_marker_y_edges(
            marker_i, markersize_i, strokewidth_i
        )
        actual_bottom = (
            yscale(result.segments[][j][2]) - yscale(y[i])
        ) * pixels_per_scaled
        actual_top = (
            yscale(result.segments[][j + 2][2]) - yscale(y[i])
        ) * pixels_per_scaled
        @assert isapprox(actual_bottom, expected_bottom; atol=0.01)
        @assert isapprox(actual_top, expected_top; atol=0.01)
    end
end

for (result, ax) in zip(linear_results, linear_axes)
    assert_error_endpoints(result, x, linear_y, linear_error_low, linear_error_high)
    assert_marker_edges(result, ax, linear_y)
end
for (result, ax) in zip(log_results, log_axes)
    assert_error_endpoints(result, x, log_y, log_error_low, log_error_high)
    assert_marker_edges(result, ax, log_y)
end
assert_error_endpoints(
    overlap_result, x_overlap, y_overlap, error_overlap, error_overlap
)
assert_marker_edges(overlap_result, ax_overlap, y_overlap)
@assert length(overlap_result.marker_overlays) == length(x_overlap) - 1
@assert length(overlap_result.errorbars) == length(x_overlap)

# Confirm that existing linear and log plots recompute their marker edges after zooming.
fig_reactive = Figure(size=(500, 250))
ax_linear = Axis(fig_reactive[1, 1])
ylims!(ax_linear, -1, 1)
linear_result = marker_errorbars!(ax_linear, [1.0], [0.0], [0.02]; markersize=30)
linear_edge_wide = linear_result.segments[][1][2]
ylims!(ax_linear, -0.1, 0.1)
linear_edge_tight = linear_result.segments[][1][2]
@assert abs(linear_edge_wide) > abs(linear_edge_tight)

ax_log = Axis(fig_reactive[1, 2]; yscale=log10)
ylims!(ax_log, 100, 10_000)
log_result = marker_errorbars!(ax_log, [1.0], [1_000.0], [200.0]; markersize=30)
log_edge_wide = log_result.segments[][1][2]
ylims!(ax_log, 500, 2_000)
log_edge_tight = log_result.segments[][1][2]
@assert abs(log10(log_edge_wide / 1_000)) > abs(log10(log_edge_tight / 1_000))
@assert log_result.segments[][2][2] ≈ 800
@assert log_result.segments[][4][2] ≈ 1_200

ax_invalid_log = Axis(Figure()[1, 1]; yscale=log10)
invalid_log_rejected = try
    marker_errorbars!(ax_invalid_log, [1.0], [100.0], [100.0])
    false
catch error
    error isa ArgumentError
end
@assert invalid_log_rejected

path_output = isempty(ARGS) ?
    joinpath(tempdir(), "benchtest_marker_errorbars.png") :
    only(ARGS)
save(path_output, fig)
println("Saved marker-errorbar bench test to: $path_output")

fig
