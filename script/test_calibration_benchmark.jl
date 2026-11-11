using Test
include(joinpath(@__DIR__, "..", "src", "calibration_benchmark.jl"))
using .CalibrationBenchmark

@testset "Continuous calibration disagreement" begin
    # The cubic's maximum absolute difference occurs strictly between anchors.
    frequency(id, t) = id == 1 ? 0.0 : t^3 - t
    slope(id, t) = id == 1 ? 0.0 : 3t^2 - 1
    result = calc_max_curve_disagreement(frequency, slope, [1, 2], [0.0, 0.4, 0.7, 1.0])
    @test result.disagreement_mhz ≈ 2 / (3sqrt(3)) atol=1e-14
    @test result.time_s ≈ 1 / sqrt(3) atol=1e-14
    @test (result.id_upper, result.id_lower) == (1, 2)

    # Piecewise-linear tails and constant curves exercise degenerate derivatives.
    frequency_linear(id, t) = id == 1 ? 0.0 : id == 2 ? t : -2.0
    slope_linear(id, t) = id == 2 ? 1.0 : 0.0
    result_linear = calc_max_curve_disagreement(frequency_linear, slope_linear,
        [1, 2, 3], [-1.0, 0.0, 2.0])
    @test result_linear.disagreement_mhz == 4.0
    @test result_linear.time_s == 2.0
    @test (result_linear.id_upper, result_linear.id_lower) == (2, 3)
    @test_throws ArgumentError calc_max_curve_disagreement(frequency, slope, [1], [0.0, 1.0])
    @test_throws ArgumentError calc_max_curve_disagreement(frequency, slope, [1, 2], [0.0, 0.0])
end
