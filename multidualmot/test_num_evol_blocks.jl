using Test
include(joinpath(@__DIR__, "dualmotcommons.jl"))

@testset "Number-evolution selectors" begin
    mktempdir() do path
        # rep is deliberately between the two experimental axes in acquisition order.
        vars = (rep=:auto, t_load=[30, 10, 20], loadcfg=[:DDM, :DIS])
        path_mat = joinpath(path, "liferes 1009 run01.mat")
        shots = [Dict("atomnum" => Float64(100t + 10r + c),
            "sigmax" => Float64(r), "sigmay" => Float64(c))
            for t in 1:3 for r in 1:10 for c in 1:2]
        # A rejected shot must retain its selected repetition slot.
        shots[5]["sigmax"] = -1.0
        matwrite(path_mat, Dict("liferes" => Dict("cres" => shots)))
        data = (; vars, var_order=(:t_load, :rep, :loadcfg),
            files=[path_mat], selector=NamedTuple())
        kwargs = (; label="selector fixture", bounds_sigmax_num=(0.0, Inf),
            bounds_sigmay_num=(0.0, Inf), num_max_num=1e6)
        original = calc_num_evol_block(data; kwargs...)
        selector = parse_num_evol_selector([
            Dict("rep" => Dict("value" => "r -> r in (3, 5, 8, 10)")),
            Dict("t_load" => Dict("value" => "t -> t >= 20")),
            Dict("loadcfg" => Dict("index" => "i -> i == 1")),
        ]; var_specs=DUALMOT_LOADING_VAR_SPECS, label="fixture")
        selected = calc_num_evol_block(merge(data, (; selector)); kwargs...)
        @test keys(selected.vars) == keys(vars)
        @test selected.vars.rep == 1:4
        @test selected.vars.t_load == [30, 20]
        @test selected.vars.loadcfg == [:DDM]
        @test selected.n_rep == 4
        @test size(selected.num_fmt) == (4, 2, 1)
        for field in (:num_fmt, :sigmax_fmt, :sigmay_fmt)
            @test isequal(getproperty(selected, field),
                getproperty(original, field)[[3, 5, 8, 10], [1, 3], [1]])
        end
        @test ismissing(selected.num_fmt[1, 1, 1])
        @test selected.n_rep_stat[:, 1] == [3, 4]
        @test selected.num_stat[1, 1] == mean([151, 181, 201])
        @test selected.std_num_stat[1, 1] == std([151, 181, 201])
        @test_throws ArgumentError calc_num_evol_block(
            merge(data, (; selector=(rep=(value=r -> false,),))); kwargs...)
        @test_throws ArgumentError calc_num_evol_block(
            merge(data, (; selector=(rep=(value=r -> r,),))); kwargs...)
    end
end

function make_num_evol_test_block(vars::NamedTuple, values::AbstractArray)
    num_fmt = Array{Union{Missing,Float64}}(values)
    (; vars, num_fmt, sigmax_fmt=copy(num_fmt), sigmay_fmt=copy(num_fmt),
        n_rep=length(vars.rep))
end

@testset "Number-evolution block union" begin
    a = make_num_evol_test_block((rep=[3, 5, 8], t_load=[30, 10], bias=[0]),
        reshape(1:6, 3, 2, 1))
    b = make_num_evol_test_block((rep=1:4, t_load=[20, 30], bias=[-0.9]),
        reshape(11:18, 4, 2, 1))
    combined = combine_num_evol_blocks([a, b]; label="union fixture")
    @test combined.vars == (rep=1:4, t_load=[10, 20, 30], bias=[0, -0.9])
    @test combined.n_rep == 4
    @test size(combined.num_fmt) == (4, 3, 2)
    for field in (:num_fmt, :sigmax_fmt, :sigmay_fmt)
        values = getproperty(combined, field)
        @test isequal(values[1:3, [3, 1], [1]], getproperty(a, field))
        @test isequal(values[:, [2, 3], [2]], getproperty(b, field))
        @test all(ismissing, values[4, :, 1])
        @test all(ismissing, values[:, 2, 1])
        @test all(ismissing, values[:, 1, 2])
    end
    reordered = [make_num_evol_test_block(
        (t_load=block.vars.t_load, rep=block.vars.rep, bias=block.vars.bias),
        permutedims(block.num_fmt, (2, 1, 3))) for block in (a, b)]
    combined_reordered = combine_num_evol_blocks(reordered; label="rep axis 2")
    @test keys(combined_reordered.vars) == (:t_load, :rep, :bias)
    @test isequal(combined_reordered.num_fmt, permutedims(combined.num_fmt, (2, 1, 3)))
    @test_throws ArgumentError combine_num_evol_blocks(NamedTuple[]; label="empty")
    @test_throws ArgumentError combine_num_evol_blocks([a, first(reordered)]; label="axes")
    @test_throws DimensionMismatch combine_num_evol_blocks(
        [merge(a, (; num_fmt=zeros(2, 2, 1)))]; label="dimensions")

    # Only t_load=30 overlaps: preserve both blocks, including a rejected shot.
    c = make_num_evol_test_block((rep=[8, 10], t_load=[20, 30], bias=[0]),
        reshape(21:24, 2, 2, 1))
    a.num_fmt[2, 1, 1] = missing
    overlap = combine_num_evol_blocks([a, b, c]; label="overlap")
    @test overlap.n_rep == 5
    @test overlap.vars.rep == 1:5
    @test isequal(overlap.num_fmt[:, 3, 1], [1, missing, 3, 23, 24])
    @test isequal(overlap.num_fmt[:, 1, 1], [4, 5, 6, missing, missing])
    @test isequal(overlap.num_fmt[:, 2, 1], [21, 22, missing, missing, missing])
    @test overlap.sigmax_fmt[:, 3, 1] == [1, 2, 3, 23, 24]
    @test all(ismissing, overlap.num_fmt[5, :, 2])
end
