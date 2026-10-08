"""Run the analysis dependencies and figure scripts needed for the Result folder."""

# Populate selected pair-comparison statistics and fitted lifetime records.
include(joinpath(@__DIR__, "anlz_isotope_pair_comparison_runner.jl"))
include(joinpath(@__DIR__, "anlz_cmot_lifetime_runner.jl"))
include(joinpath(@__DIR__, "anlz_mot_lifetime_runner.jl"))
include(joinpath(@__DIR__, "anlz_lifetime_pair_comparison.jl"))
include(joinpath(@__DIR__, "anlz_odt_ib_runner.jl"))
include(joinpath(@__DIR__, "anlz_odt_cmot_comparison.jl"))

# Export all selected Result figures and the separate legend PDFs.
include(joinpath(@__DIR__, "draw_selected_results.jl"))
include(joinpath(@__DIR__, "draw_selected_pair_results.jl"))
include(joinpath(@__DIR__, "draw_result_legends.jl"))
