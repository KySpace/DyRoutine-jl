using YAML
using MAT
using CairoMakie
using Printf

include(joinpath(@__DIR__, "..", "src", "graphics.jl"))
CairoMakie.activate!()

path_folder = length(ARGS) >= 1 ? ARGS[1] :
    raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\insitu samples\162-164 MOT"
path_plot = length(ARGS) >= 2 ? ARGS[2] : joinpath(path_folder, "MOT in-situ samples.svg")
path_plot_misc = length(ARGS) >= 3 ? ARGS[3] :
    joinpath(path_folder, "MOT in-situ samples miscibility.svg")

config = YAML.load_file(joinpath(path_folder, "config.yaml"))
config isa AbstractVector && length(config) == 1 ||
    throw(ArgumentError("expected one processing entry in config.yaml"))
blocks = get(only(config), "data", nothing)
blocks isa AbstractVector && length(blocks) == 1 ||
    throw(ArgumentError("expected one rectangular data block"))
block = only(blocks)

name_acq = Symbol.(block["vars"])
Set(name_acq) == Set((:f_626_MZM, :rep, :loadcfg, :istp)) && length(name_acq) == 4 ||
    throw(ArgumentError("vars must contain f_626_MZM, rep, loadcfg, and istp once each"))
val_vars = Dict(name => block[string(name)] for name in name_acq)
all(values -> values isa AbstractVector && !isempty(values), values(val_vars)) ||
    throw(ArgumentError("every configured variable must have a nonempty value list"))
val_vars[:loadcfg] == ["SCS"] ||
    throw(ArgumentError("expected the SCS load configuration"))
string.(val_vars[:istp]) == ["162", "164"] ||
    throw(ArgumentError("isotopes must be ordered 162, 164"))

sources = block["source"]
sources isa AbstractVector && length(sources) == 1 ||
    throw(ArgumentError("expected one liferes source file"))
name_allimg = replace(only(sources), r"^liferes(?= )" => "allimg")
name_allimg != only(sources) ||
    throw(ArgumentError("source must be named 'liferes <run>.mat'"))
path_allimg = joinpath(path_folder, name_allimg)
images = vec(matread(path_allimg)["allimg"])
n_dims_acq = length.(getindex.(Ref(val_vars), name_acq))
length(images) == prod(n_dims_acq) ||
    throw(DimensionMismatch("$(length(images)) images, expected $(prod(n_dims_acq)) from config.yaml"))
all(img -> img isa AbstractMatrix{<:Real}, images) ||
    throw(ArgumentError("allimg cells must be numeric matrices"))
size_image = size(first(images))
all(img -> size(img) == size_image && all(isfinite, img), images) ||
    throw(ArgumentError("allimg cells must have equal dimensions and finite values"))

# Acquisition variables run slowest to fastest in config.yaml.
images_acq = reshape(images, Tuple(reverse(n_dims_acq))) |>
    data -> permutedims(data, reverse(1:length(n_dims_acq)))
name = (:f_626_MZM, :rep, :loadcfg, :istp)
images_fmt = permutedims(images_acq, indexin(collect(name), collect(name_acq)))
val_f = val_vars[:f_626_MZM]
val_rep = val_vars[:rep]
max_image = maximum(maximum, images)
max_image > 0 || throw(ArgumentError("all images are nonpositive"))

fig = Figure(size=(1450, 2800), fontsize=17, backgroundcolor=:white)
Label(fig[1, 1], "f_626_MZM")
for (idx_rep, rep) in enumerate(val_rep)
    col_left = 2idx_rep
    Label(fig[1, col_left:col_left+1], "rep = $rep", fontsize=20)
    Label(fig[2, col_left], "162", fontsize=17)
    Label(fig[2, col_left+1], "164", fontsize=17)
end

hue_istp = (21, 259)
alpha_base = 0.0
thres_alpha = 0.01
clrmaps = map(hue -> gen_clrmap_solo(hue; alpha_base, thres_alpha), hue_istp)
for (idx_f, f) in enumerate(val_f)
    row = idx_f + 2
    Label(fig[row, 1], @sprintf("%.4f", f), fontsize=16)
    for idx_rep in eachindex(val_rep), idx_istp in 1:2
        col = 2idx_rep + idx_istp - 1
        ax = Axis(fig[row, col]; aspect=DataAspect(), yreversed=true)
        hidedecorations!(ax)
        hidespines!(ax)
        heatmap!(ax, images_fmt[idx_f, idx_rep, 1, idx_istp]';
            colormap=clrmaps[idx_istp], colorrange=(0, max_image), rasterize=true)
    end
end
colgap!(fig.layout, 3, 18)
colgap!(fig.layout, 5, 18)
colsize!(fig.layout, 1, Fixed(125))
for col in 2:7
    colsize!(fig.layout, col, Fixed(180))
end
rowsize!(fig.layout, 1, Fixed(34))
rowsize!(fig.layout, 2, Fixed(30))
for row in 3:length(val_f)+2
    rowsize!(fig.layout, row, Fixed(180))
end
save(path_plot, fig)
println("Saved $path_plot ($(length(images)) images; shared color range 0–$max_image)")

fig_misc = Figure(size=(850, 2800), fontsize=17, backgroundcolor=:white)
Label(fig_misc[1, 1], "f_626_MZM")
for (idx_rep, rep) in enumerate(val_rep)
    col = idx_rep + 1
    Label(fig_misc[1, col], "rep = $rep", fontsize=20)
    Label(fig_misc[2, col], "162 + 164", fontsize=17)
end
for (idx_f, f) in enumerate(val_f)
    row = idx_f + 2
    Label(fig_misc[row, 1], @sprintf("%.4f", f), fontsize=16)
    for idx_rep in eachindex(val_rep)
        ax = Axis(fig_misc[row, idx_rep+1]; aspect=DataAspect(), yreversed=true)
        hidedecorations!(ax)
        hidespines!(ax)
        clr_misc = to_miscibility_clr(
            images_fmt[idx_f, idx_rep, 1, 1],
            images_fmt[idx_f, idx_rep, 1, 2],
            hue_istp...;
            max=max_image,
            alpha_base,
            thres_alpha,
        )
        heatmap!(ax, clr_misc'; rasterize=true)
    end
end
colsize!(fig_misc.layout, 1, Fixed(125))
for col in 2:4
    colsize!(fig_misc.layout, col, Fixed(180))
end
rowsize!(fig_misc.layout, 1, Fixed(34))
rowsize!(fig_misc.layout, 2, Fixed(30))
for row in 3:length(val_f)+2
    rowsize!(fig_misc.layout, row, Fixed(180))
end
save(path_plot_misc, fig_misc)
println("Saved $path_plot_misc ($(length(val_f) * length(val_rep)) paired images)")
