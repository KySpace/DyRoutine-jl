path_root = raw"C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\Simulations"
dir_test = raw"01.[2025.06.01]"

path_test = joinpath(path_root, dir_test)

pattern_filename_data = Regex("nxy_t=(<?idx_time>%d+).mat")

filenames_data = readdir(path_test) |> fs -> filter(f -> occursin(pattern_filename_data, f), fs)
for file in readdir(path_test)
end
