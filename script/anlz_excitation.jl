# Compatibility wrapper for the split excitation workflow.
# The runner now includes these stages directly.
include(joinpath(@__DIR__, "anlz_excitation_extr.jl"))
include(joinpath(@__DIR__, "anlz_excitation_corr.jl"))
include(joinpath(@__DIR__, "anlz_excitation_vslz.jl"))
