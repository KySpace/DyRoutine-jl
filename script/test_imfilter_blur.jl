using GLMakie
using ImageFiltering
using Printf

GLMakie.activate!()

px_in_um = 6.5 / 22.06
step_posi = px_in_um
smwh_roi = (50, 100)
x_posi, y_posi = smwh_roi |> s -> map(u -> (-u:1:u), s) |> xy -> xy .* step_posi

# Array rows are axial (y, height); columns are radial (x, width).
# The radial x profile is the modulated narrow direction.
# The radial blob is broad and modulated; the axial blob is narrow.
radial_sigma_um = 1.8
axial_sigma_um = 6.0
wavelength_um = 3.1

gaussian_elliptical(y::Real, x::Real) = exp(-0.5 * ((y / axial_sigma_um)^2 + (x / radial_sigma_um)^2))
axial_modulation(x::Real) = 0.5 * (1 + cos(2π * x / wavelength_um))

dens_dummy = [
    gaussian_elliptical(y, x) * axial_modulation(y)
    for y in y_posi, x in x_posi
]

function filter_gaussian(dens::AbstractMatrix{<:Real}, sigma_px::Real)
    sigma_px >= 0 || throw(ArgumentError("Gaussian sigma must be nonnegative, got $sigma_px."))
    sigma_px == 0 && return Float64.(dens)
    return imfilter(dens, Kernel.gaussian(sigma_px))
end

sigma_px_initial = 0.0
dens_filtered = Observable(filter_gaussian(dens_dummy, sigma_px_initial)')
sigma_px = Observable(sigma_px_initial)

fig = Figure(size=(920, 700))
ax = Axis(
    fig[1, 1];
    xlabel="radial x (μm)",
    ylabel="axial y (μm)",
    title=lift(sigma_px) do σ
        @sprintf("imfilter Gaussian blur: σ = %.2f px (%.3f μm)", σ, σ * px_in_um)
    end,
    aspect=DataAspect(),
)
heatmap!(ax, x_posi, y_posi, dens_filtered; colormap=:viridis, colorrange=(0, maximum(dens_dummy)))

slider = Slider(fig[2, 1]; range=0:0.25:12, startvalue=sigma_px_initial, tellwidth=false)
Label(fig[3, 1]; text="Gaussian kernel σ (pixels)", halign=:center)

on(slider.value) do σ
    sigma_px[] = σ
    dens_filtered[] = filter_gaussian(dens_dummy, σ)'
end

display(fig)
fig
