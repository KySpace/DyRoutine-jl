using HDF5

path = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\2026-03\0324\run64\d0324r64.h5"

name = ["repeat", "t_hold", "istp"]
val = (
    collect(1:3),
    collect(6:2:200),
    [5, 0],
)
variation = length(val[1]) * length(val[2]) * length(val[3])

h5open(path, "r") do f
    global dens = read(f["/od"])
end

ndims(dens) == 3 || error("Expected /od to have 3 dimensions, got $(ndims(dens)).")
size(dens) == (201, 401, variation) || error(
    "Expected /od size to be (201, 401, $variation), got $(size(dens)).",
)

width, height, _ = size(dens)

dens_reshaped = reshape(dens, width, height, length(val[3]), length(val[2]), length(val[1]))
dens_by_variation = permutedims(dens_reshaped, (5, 4, 3, 2, 1))

println("name = ", name)
println("val = ", val)
println("variation = ", variation)
println("raw /dens size = ", size(dens))
println("dens_by_variation size = ", size(dens_by_variation))
println("Indexing order: dens_by_variation[repeat, t_hold, istp, height, width]")
