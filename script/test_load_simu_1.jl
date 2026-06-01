using MAT

path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"01.[2025.06.01]"
step_t = 0.1798
step_space = 0.2613

path_test = joinpath(path_root, dir_test)

pattern_filename_data = Regex(raw"nxy_t=(?<idx_time>\d+).mat")

filenames_data = readdir(path_test) |> fs -> filter(f -> occursin(pattern_filename_data, f), fs)
ids_t = Vector{Int}(undef, length(filenames_data))
dens_raw = Array{Float64}(undef, (length(filenames_data), 2, 256, 256))

for (i, fn) in enumerate(filenames_data)
    file = matopen(joinpath(path_test, fn))
    ids_t[i] = parse(Int, match(pattern_filename_data, fn)["idx_time"])
    dens_raw[i, 1, :, :] = read(file, "n1")
    dens_raw[i, 2, :, :] = read(file, "n2")
    close(file)
end

perm_t = sortperm(ids_t)
ids_t = ids_t[perm_t]
dens_raw = dens_raw[perm_t, :, :, :]
istp = ["162", "164"]

var_vals = (;
    t_hold=ids_t .* step_t,
    istp,
)
imgs = [copy(@view dens_raw[t, i, :, :]) for t in axes(dens_raw, 1), i in axes(dens_raw, 2)]
