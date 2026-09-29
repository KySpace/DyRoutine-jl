include(joinpath(@__DIR__, "dualmotcommons.jl"))

path_root = raw"C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data\MOT loading 421"
val_pair = ["162-164", "160-162", "161-162", "161-164", "163-164",
    "162-163", "161-163"]
var_specs_num = DUALMOT_LOADING_VAR_SPECS
key_x_num = :t_load
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8
γ_active = 0.43

runinfos_grouped = [read_num_evol_runinfos(path_root, pair;
    var_specs=var_specs_num,
    validate_vars=validate_dualmot_vars,
) for pair in val_pair]
runinfos = vcat(runinfos_grouped...)
ids_runinfo = eachindex(runinfos)
figure_data_sheets = Dict{String,Vector{Pair{String,Matrix{Any}}}}()
plot_num_evol_base = dualmot_num_evol_plot_spec("MOT";
    key_x=key_x_num,
    scale_x=1,
    xlabel=bias -> iszero(bias) ?
        "Effective loading time (s)" : "Equivalent loading time (s)",
    transform_x=(values, condition, bias, idx_istp) -> begin
        if iszero(bias)
            values .* γ_active
        elseif condition.loadcfg == :DCS
            idx_istp in (1, 2) ||
                throw(ArgumentError("DCS time scaling requires isotope index 1 or 2"))
            values ./ (idx_istp == 1 ? 1 - bias : 1 + bias)
        else
            values
        end
    end,
    xautolimits=(condition, bias) -> iszero(bias) || condition.loadcfg != :DCS,
    yautolimits=(condition, bias) -> iszero(bias) || condition.loadcfg != :DCS,
    ylabel="CMOT number",
    file_head="MOT.loading.421",
    legend_position=:rb,
    loadcfg_plot=(:DDM, :DIS, :DCS),
)
plot_num_evol_target = nothing

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_root, runinfo.folder)
    local count_folder_tags = count(info -> info.folder == runinfo.folder, runinfos)
    local idx_folder_tag = count(info -> info.folder == runinfo.folder,
        @view runinfos[firstindex(runinfos):idx_runinfo])
    local style_tag = count_folder_tags == 1 || idx_folder_tag == 1 ?
        nothing : filename_token(runinfo.tag)
    local filename_run = (scale, bias) -> begin
        style = isnothing(style_tag) ? string(scale) : "$(scale).$(style_tag)"
        "[MOT.loading.421].[$style].[$(balance_tag(bias))-balanced]"
    end
    local folder_run = runinfo.folder
    local limits_linear_run = isdefined(@__MODULE__, :limits_loading_linear) ?
        bias -> get(limits_loading_linear,
            (folder_run, Symbol(balance_tag(bias))), nothing) : nothing
    global plot_num_evol = merge(plot_num_evol_base, (;
        filename=filename_run,
        limits_linear=limits_linear_run))
    println("Processing: $tag_head")
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
    include(joinpath(@__DIR__, "anlz_num_evol_output.jl"))
end
write_pair_figure_workbooks(path_root, figure_data_sheets, "MOT loading 421.xlsx")
