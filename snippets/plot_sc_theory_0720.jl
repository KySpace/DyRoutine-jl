using CairoMakie
using DelimitedFiles: readdlm
using Colors: Oklch

path_data = raw"C:\Users\ky\OneDrive\Postdoc\Midterm\BinarySupersolid\SC_theory_0720.csv"
path_output_dir = dirname(path_data)

# Isotope colors use the project's 164 and 162 hues; line 3 is green.
hue_theme_istp = Dict("162" => 11.0, "164" => 250.0)
line_colors = (
    Oklch(0.52, 0.14, hue_theme_istp["162"]),
    Oklch(0.52, 0.14, hue_theme_istp["164"]),
    Oklch(0.52, 0.14, 145.0),
)
line_labels = ("162", "164", "162-164")
xlabel = "x"
ylabel = "Value"
x_limits = (9.3, 9.55)
y_limits = (0, 250)
x_major_ticks = 9.3:0.05:9.55
x_minor_ticks = [
    x for x in 9.31:0.01:9.54
    if !any(x_major -> isapprox(x, x_major; atol=1e-10), x_major_ticks)
]

data_text = readdlm(path_data, ',', String)
size(data_text, 2) == 4 || throw(ArgumentError(
    "Expected 4 CSV columns (x plus 3 y values), got $(size(data_text, 2)) in $path_data",
))
data = map(data_text) do token
    token = strip(token)
    token == "ComplexInfinity" ? NaN : parse(Float64, token)
end

fig = Figure()
ax = Axis(
    fig[1, 1];
    width=800,
    height=240,
    xlabel,
    ylabel,
    limits=(x_limits, y_limits),
    xticks=x_major_ticks,
    xminorticks=x_minor_ticks,
    xminorticksvisible=true,
)

for idx_line in 1:3
    lines!(
        ax,
        data[:, 1],
        data[:, idx_line + 1];
        color=line_colors[idx_line],
        linewidth=2,
        label=line_labels[idx_line],
    )
end

axislegend(ax; position=:lt, framevisible=false)
fig |> resize_to_layout!
save(joinpath(path_output_dir, "SC_theory_0720_lines.svg"), fig)
save(joinpath(path_output_dir, "SC_theory_0720_lines.pdf"), fig)
