include(joinpath(@__DIR__, "dualmotcommons.jl"))
include(joinpath(@__DIR__, "anlz_pair_comparison_plot_helpers.jl"))
using CSV

path_dataset = joinpath(path_root, "CMOT lifetime")
path_data_root = path_root
val_pair = ["162-164", "160-162", "161-162", "161-164", "163-164", "162-163", "161-163"]
var_specs_num = DUALMOT_VAR_SPECS
key_x_num = :t_hold
bounds_sigmax_num = (4e-4, Inf)
bounds_sigmay_num = (4e-4, Inf)
num_max_num = 2e8
fit_selector_num = times -> times .<= 2.0
fit_bounds_n0_num = (0.5, 2.0)
fit_bounds_kappa_num = (eps(Float64), Inf)
fit_modes_num = (:full, :kappa)
pair_istp_offset = 0.05
balance_rows = collect(CSV.File(joinpath(path_data_root, "MOT loading balance", "balance.csv");
    header=["Pair", "tbiasmot"], skipto=2, stripwhitespace=true))
balance_bias_by_pair = Dict(join(sort(parse.(Int, split(strip(string(row.Pair)), "-"))), "-") =>
    Float64(row.tbiasmot) for row in balance_rows)
Set(keys(balance_bias_by_pair)) == Set(val_pair) || throw(ArgumentError(
    "MOT loading balance/balance.csv must contain exactly the configured isotope pairs"))
balance_biases = Dict(
    "t-balanced" => Dict(pair => 0.0 for pair in val_pair),
    "n-balanced" => balance_bias_by_pair,
)

runinfos_grouped = [read_num_evol_runinfos(path_dataset, pair;
    var_specs=var_specs_num,
    validate_vars=validate_dualmot_vars,
) for pair in val_pair]
runinfos = vcat(runinfos_grouped...)
ids_runinfo = eachindex(runinfos)
plot_num_evol = dualmot_lifetime_plot_spec("CMOT";
    field_num="atomnum",
    key_x=key_x_num,
    scale_x=1000,
    xlabel="CMOT holding time (s)",
    ylabel="CMOT number",
)
plot_num_evol_target = nothing
fit_records_num_decay = NamedTuple[]
formats_output = ("svg", "png")

for idx_runinfo_iter in ids_runinfo
    global idx_runinfo = idx_runinfo_iter
    global runinfo = runinfos[idx_runinfo]
    global tag_head = "$(runinfo.folder) $(runinfo.tag)"
    global path_output = joinpath(path_dataset, runinfo.folder)
    println("Processing: $tag_head")
    GC.gc()
    include(joinpath(@__DIR__, "anlz_num_evol.jl"))
    for mode in fit_modes_num
        global fit_num_decay_config = (
            mode=mode,
            bounds_n0=fit_bounds_n0_num,
            bounds_kappa=fit_bounds_kappa_num,
            selector=fit_selector_num,
        )
        global fit_model_tag = string(mode)
        global path_output = joinpath(path_dataset, runinfo.folder, "alternatives")
        mkpath(path_output)
        include(joinpath(@__DIR__, "anlz_num_evol_decay.jl"))
    end
end

path_fit_results = save_num_decay_results(
    path_dataset, "CMOT", :alternatives, fit_records_num_decay)

function selected_model_records(records, variant::AbstractString)
    biases = balance_biases[variant]
    [record for record in records if haskey(biases, record.pair) &&
        isapprox(record.panel, biases[record.pair]; atol=1e-10, rtol=0)]
end

function draw_cmot_model_parameter(records, parameter::Symbol, mode::Symbol,
    variant::AbstractString, limits::Tuple)
    available = parameter == :n0 || (parameter == :tau && mode == :full) ||
        (parameter == :kappa && mode in (:full, :kappa))
    available || return nothing
    variant_records = selected_model_records(records, variant)
    model_records = filter(record -> record.mode == mode, variant_records)
    std_parameter = Symbol("std_", string(parameter))
    fig = Figure(size=(340, 175), fontsize=8, figure_padding=1)
    ax = Axis(fig[1, 1];
        xticks=(eachindex(val_pair), val_pair),
        xlabel="Isotope pair",
        ylabel=parameter == :n0 ? "N₀ (atom)" : parameter == :tau ? "τ (s)" : "κ (atom⁻¹ s⁻¹)",
        title="CMOT $parameter · :$mode · $variant",
        yscale=log10,
        dualmot_axis_kwargs(; log_y=true, text_size=8)...,
    )
    ylims!(ax, limits...)
    plotted = Set{Symbol}()
    for (idx_pair, pair) in enumerate(val_pair)
        pair_records = filter(record -> record.pair == pair, model_records)
        seen_conditions = Set{Tuple{Symbol,Symbol}}()
        for record in pair_records
            condition = (record.loadcfg, record.istp)
            condition in seen_conditions && throw(ArgumentError(
                "$pair $variant has duplicate $mode fit records for $condition"))
            push!(seen_conditions, condition)
            istp_list = Symbol.(split(pair, "-"))
            idx_istp = findfirst(==(record.istp), istp_list)
            isnothing(idx_istp) && continue
            value = Float64(getproperty(record, parameter))
            error = Float64(getproperty(record, std_parameter))
            isfinite(value) && value > 0 || continue
            xpos = pair_istp_position(idx_pair, pair, record.istp)
            style = dualmot_curve_style((; loadcfg=record.loadcfg, istp=record.istp))
            label = record.loadcfg in plotted ? nothing : string(record.loadcfg)
            push!(plotted, record.loadcfg)
            if isfinite(error) && error >= 0
                marker_errorbars_log!(ax, [xpos], [value], [error];
                    floor=limits[1],
                    style.marker_options..., style.errorbar_options...,
                    marker=style.marker, label)
            else
                scatter!(ax, [xpos], [value]; style.marker_options...,
                    marker=style.marker, label)
            end
        end
    end
    stem = "[CMOT.models].[$(parameter).$(mode)].[$variant]"
    for format in formats_output
        save_options = format == "png" ? (; px_per_unit=4) : (;)
        save(joinpath(path_comparison, "$stem.$format"), fig; save_options...)
    end
    fig
end

path_comparison = joinpath(path_data_root, "Isotope pair comparison")
mkpath(path_comparison)
model_records_by_balance = Dict(variant => selected_model_records(fit_records_num_decay, variant)
    for variant in keys(balance_biases))
selected_records_all_balances = vcat(values(model_records_by_balance)...)
for parameter in (:n0, :tau, :kappa)
    available_modes = parameter == :tau ? (:full,) : fit_modes_num
    comparison_records = filter(record ->
        record.mode in available_modes && isfinite(getproperty(record, parameter)) &&
            getproperty(record, parameter) > 0, selected_records_all_balances)
    comparison_values = Float64.(getproperty.(comparison_records, parameter))
    comparison_errors = Float64.(getproperty.(comparison_records, Symbol("std_", string(parameter))))
    limits_parameter = log_plot_limits(comparison_values, comparison_errors; cap_factor=5)
    for variant in keys(balance_biases), mode in available_modes
        draw_cmot_model_parameter(fit_records_num_decay, parameter, mode,
            variant, limits_parameter)
    end
end

# Keep the model-comparison panels uncluttered; put the complete marker key in
# its own table cell, using the same isotope/load-configuration styles.
legend_figure = Figure(size=(95, 75), fontsize=8, figure_padding=2)
draw_number_style_key!(legend_figure, (:SCS, :DCS, :DIS, :DDM))
save(joinpath(path_comparison, "[CMOT.models].[legend].[all].png"), legend_figure;
    px_per_unit=4)

windows_powershell = raw"C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe"
onenote_helper = joinpath(@__DIR__, "..", "helpers", "make_multi_dual_mot_onenote.ps1")
run(Cmd([windows_powershell, "-NoProfile", "-ExecutionPolicy", "Bypass",
    "-File", onenote_helper, "-DataRoot", path_data_root,
    "-OutputPath", joinpath(path_data_root, "multi_dual_mot_table.one"),
    "-CmotModelComparison"]))
