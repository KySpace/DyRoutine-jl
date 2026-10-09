"""Run the analysis dependencies and figure scripts needed for the Result folder."""

# Load the comparison statistics and previously saved lifetime fits needed by Result.
withenv("DUALMOT_RESULTS_ONLY" => "true",
    "DUALMOT_SKIP_PAIR_COMPARISON_EXPORTS" => "true") do
    include(joinpath(@__DIR__, "anlz_isotope_pair_comparison_runner.jl"))
    include(joinpath(@__DIR__, "anlz_lifetime_pair_comparison.jl"))
    include(joinpath(@__DIR__, "anlz_odt_cmot_comparison.jl"))
end

# Export all selected Result figures and the separate legend PDFs.
include(joinpath(@__DIR__, "draw_selected_results.jl"))
include(joinpath(@__DIR__, "draw_selected_pair_results.jl"))
if get(ENV, "DUALMOT_SKIP_RESULT_LEGENDS", "false") != "true"
    include(joinpath(@__DIR__, "draw_result_legends.jl"))
end
include(joinpath(@__DIR__, "draw_selected_additional_results.jl"))
if get(ENV, "DUALMOT_SKIP_RESULT_LEGENDS", "false") != "true"
    include(joinpath(@__DIR__, "draw_result_curve_legends.jl"))
end
