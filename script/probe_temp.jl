using HDF5
include(joinpath(@__DIR__, "..", "src", "pershot.jl"))

path = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\2026-03\0324\run64\d0324r64.h5"
corner_height = 10
corner_width = 10

name = ["repeat", "t_hold", "istp"]
val = (
    collect(1:3),
    collect(6:2:200),
    [5, 0],
)
variation = length(val[1]) * length(val[2]) * length(val[3])

h5open(path, "r") do f
    global dens = f["/od"] |>
                  read |>
                  x -> permutedims(x, (3, 2, 1)) |>
                  x -> stack(
                      map(d -> subtract_corner_mean(d, corner_height, corner_width), eachslice(x; dims=1));
                      dims=1,
                  )
end

ndims(dens) == 3 || error("Expected /od to have 3 dimensions, got $(ndims(dens)).")
size(dens) == (variation, 401, 201) || error(
    "Expected permuted /od size to be ($variation, 401, 201), got $(size(dens)).",
)

_, height, width = size(dens)

dens_reshaped = reshape(dens, length(val[3]), length(val[2]), length(val[1]), height, width)
dens_by_variation = permutedims(dens_reshaped, (3, 2, 1, 4, 5))

println("name = ", name)
println("val = ", val)
println("variation = ", variation)
println("permuted /dens size = ", size(dens))
println("corner subtraction = ", (corner_height, corner_width))
println("dens_by_variation size = ", size(dens_by_variation))
println("Indexing order: dens_by_variation[repeat, t_hold, istp, height, width]")
