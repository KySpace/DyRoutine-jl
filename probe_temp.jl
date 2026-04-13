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
    global od = read(f["/od"])
end

ndims(od) == 3 || error("Expected /od to have 3 dimensions, got $(ndims(od)).")
size(od) == (201, 401, variation) || error(
    "Expected /od size to be (201, 401, $variation), got $(size(od)).",
)

width, height, _ = size(od)

od_reshaped = reshape(od, width, height, length(val[3]), length(val[2]), length(val[1]))
od_by_variation = permutedims(od_reshaped, (5, 4, 3, 2, 1))

println("name = ", name)
println("val = ", val)
println("variation = ", variation)
println("raw /od size = ", size(od))
println("od_by_variation size = ", size(od_by_variation))
println("Indexing order: od_by_variation[repeat, t_hold, istp, height, width]")
