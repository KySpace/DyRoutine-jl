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

function loading_421_curve_stats(name_stat::Tuple, vars::NamedTuple,
    values::AbstractArray, idx_bias::Integer, loadcfg::Symbol, istp::Symbol)
    idx_axis = Dict(key => idx for (idx, key) in enumerate(name_stat))
    indices = Any[Colon() for _ in name_stat]
    indices[idx_axis[:β_MOT]] = idx_bias
    indices[idx_axis[:loadcfg]] = findfirst(==(loadcfg), vars.loadcfg)
    indices[idx_axis[:istp]] = findfirst(==(istp), vars.istp)
    vec(@view values[indices...])
end

function save_loading_ratio_figures(runinfo::NamedTuple, vars::NamedTuple,
    name_stat::Tuple, num_stat::AbstractArray, std_num_stat::AbstractArray,
    n_rep_stat::AbstractArray, path_output::AbstractString)
    val_t_load = vars.t_load
    val_istp = vars.istp
    for bias in vars.β_MOT
        idx_bias = findfirst(==(bias), vars.β_MOT)
        for (numerator, denominator) in ((:DIS, :DCS), (:DDM, :DIS))
            bias_text = @sprintf("%+.2f", bias)
            fig = Figure(size=plot_num_evol_base.size,
                fontsize=plot_num_evol_base.fontsize, figure_padding=1)
            ax = Axis(fig[1, 1];
                xlabel="Loading time (s)",
                ylabel="$numerator / $denominator number",
                title=rich("$(runinfo.folder) $(runinfo.tag) · $numerator/$denominator · ",
                    rich("β", font=:italic), subscript("MOT"), " = ", bias_text),
                width=plot_num_evol_base.frame_size[1],
                height=plot_num_evol_base.frame_size[2],
                yscale=log10,
                dualmot_axis_kwargs(log_y=true)...,
            )
            hlines!(ax, [1.0]; color=RGBAf(0, 0, 0, 0.4),
                linestyle=:dash, linewidth=0.5)
            for istp in val_istp
                style = dualmot_curve_style((; loadcfg=numerator, istp))
                scatter!(ax, Float64[], Float64[]; style.marker_options...,
                    marker=style.marker, label=string(istp))
            end
            ratio_sheet = Matrix{Any}(undef, length(val_t_load) + 1,
                2 + 10 * length(val_istp))
            ratio_sheet[1, 1:2] .= ("t_load (raw)", "t_load plotted (s)")
            for (idx_istp, istp) in enumerate(val_istp)
                col = 3 + 10 * (idx_istp - 1)
                ratio_sheet[1, col:col + 9] .= (
                    "mean $numerator [$istp]", "std $numerator [$istp]",
                    "n_rep $numerator [$istp]", "source $numerator [$istp]",
                    "mean $denominator [$istp]", "std $denominator [$istp]",
                    "n_rep $denominator [$istp]", "source $denominator [$istp]",
                    "ratio $numerator/$denominator [$istp]",
                    "std ratio [$istp]")
            end
            for (idx_t, t_load) in enumerate(val_t_load)
                ratio_sheet[idx_t + 1, 1:2] .= (t_load, t_load)
            end
            for istp in val_istp
                num = loading_421_curve_stats(name_stat, vars, num_stat,
                    idx_bias, numerator, istp)
                std_num = loading_421_curve_stats(name_stat, vars, std_num_stat,
                    idx_bias, numerator, istp)
                n_num = loading_421_curve_stats(name_stat, vars, n_rep_stat,
                    idx_bias, numerator, istp)
                den = loading_421_curve_stats(name_stat, vars, num_stat,
                    idx_bias, denominator, istp)
                std_den = loading_421_curve_stats(name_stat, vars, std_num_stat,
                    idx_bias, denominator, istp)
                n_den = loading_421_curve_stats(name_stat, vars, n_rep_stat,
                    idx_bias, denominator, istp)
                ratio = fill(NaN, length(val_t_load))
                valid = isfinite.(num) .& isfinite.(den) .& (den .> 0) .& (num .> 0)
                ratio[valid] .= num[valid] ./ den[valid]
                ratio_std = fill(NaN, length(val_t_load))
                valid_std = valid .& isfinite.(std_num) .& isfinite.(std_den)
                ratio_std[valid_std] .= sqrt.((std_num[valid_std] ./ den[valid_std]).^2 .+
                    (num[valid_std] .* std_den[valid_std] ./ den[valid_std].^2).^2)

                style = dualmot_curve_style((; loadcfg=numerator, istp))
                mask_line = isfinite.(ratio) .& (ratio .> 0)
                count(mask_line) >= 2 && lines!(ax, val_t_load[mask_line],
                    ratio[mask_line]; style.line_options...)
                mask_error = mask_line .& isfinite.(ratio_std) .&
                    (ratio .- ratio_std .> 0)
                mask_marker_only = mask_line .& .!mask_error
                if any(mask_error)
                    marker_errorbars!(ax, val_t_load[mask_error], ratio[mask_error],
                        ratio_std[mask_error]; style.marker_options...,
                        style.errorbar_options..., marker=style.marker,
                        label=nothing)
                end
                any(mask_marker_only) && scatter!(ax,
                    val_t_load[mask_marker_only], ratio[mask_marker_only];
                    style.marker_options..., marker=style.marker,
                    label=nothing)

                idx_sheet_col = 3 + 10 * (findfirst(==(istp), val_istp) - 1)
                for idx_t in eachindex(val_t_load)
                    condition_num = (; β_MOT=bias, t_load=val_t_load[idx_t],
                        loadcfg=numerator, istp)
                    condition_den = (; β_MOT=bias, t_load=val_t_load[idx_t],
                        loadcfg=denominator, istp)
                    ratio_sheet[idx_t + 1, idx_sheet_col:idx_sheet_col + 9] .= (
                        num[idx_t], std_num[idx_t], n_num[idx_t] == 0 ? NaN : n_num[idx_t],
                        figure_xlsx_source(runinfo, condition_num), den[idx_t],
                        std_den[idx_t], n_den[idx_t] == 0 ? NaN : n_den[idx_t],
                        figure_xlsx_source(runinfo, condition_den), ratio[idx_t],
                        ratio_std[idx_t])
                end
            end
            set_time_minor_ticks!(ax)
            axislegend(ax; position=:rt, DUALMOT_LEGEND_OPTIONS...)
            balance = "$(balance_tag(bias))-balanced"
            variant = "$(balance).$(filename_token(runinfo.tag))"
            style_tag = "ratio.$numerator-$denominator"
            stem = "[MOT.loading.421].[$style_tag].[$variant]"
            for format in plot_num_evol_base.formats
                save_options = format == "png" ? (; px_per_unit=4) : (;)
                save(joinpath(path_output, "$stem.$format"), fig; save_options...)
            end
            sheet_name = figure_xlsx_sheet_name(
                "421 $(filename_token(runinfo.tag)) $(balance_tag(bias)) $numerator-$denominator")
            push!(get!(figure_data_sheets, runinfo.folder,
                Pair{String,Matrix{Any}}[]), sheet_name => ratio_sheet)
        end
    end
    nothing
end

plot_num_evol_target = nothing

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_root, runinfo.folder)
    local filename_run = (scale, bias) ->
        "[MOT.loading.421].[$scale].[$(balance_tag(bias))-balanced.$(filename_token(runinfo.tag))]"
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
    selected_runinfo = select_for_load421_runinfo(
        only(filter(infos -> infos[1].folder == runinfo.folder, runinfos_grouped)),
        runinfo.folder)
    runinfo.tag == selected_runinfo.tag && save_loading_ratio_figures(
        runinfo, vars, name_stat, num_stat, std_num_stat, n_rep_stat, path_output)
end
write_pair_figure_workbooks(path_root, figure_data_sheets, "MOT loading 421.xlsx")
