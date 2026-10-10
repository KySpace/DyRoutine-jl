include(joinpath(@__DIR__, "preprocess_multidualmot.jl"))
include(joinpath(@__DIR__, "anlz_mot_loading_421_runner.jl"))
include(joinpath(@__DIR__, "anlz_mot_loading_421_monofreq_runner.jl"))
include(joinpath(@__DIR__, "anlz_mot_loading_626_runner.jl"))
include(joinpath(@__DIR__, "anlz_mot_loading_balance_runner.jl"))
include(joinpath(@__DIR__, "anlz_cmot_lifetime_runner.jl"))
include(joinpath(@__DIR__, "anlz_mot_lifetime_runner.jl"))
include(joinpath(@__DIR__, "anlz_isotope_pair_comparison_runner.jl"))
include(joinpath(@__DIR__, "anlz_lifetime_pair_comparison.jl"))
include(joinpath(@__DIR__, "anlz_odt_ib_runner.jl"))
include(joinpath(@__DIR__, "anlz_odt_cmot_comparison.jl"))
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

isdefined(@__MODULE__, :main_multi_dual_mot_table) ||
    include(joinpath(@__DIR__, "..", "helpers", "make_multi_dual_mot_table.jl"))
main_multi_dual_mot_table([path_root, "--formats=svg,one"])
