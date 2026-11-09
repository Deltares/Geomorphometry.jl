using Geomorphometry
using Geomorphometry: PriorityFlood, Zhou, resolveflats, flowdirection, maxreceivers
using Test
using Random

"""Random DEM with broad features, so it contains real depressions and flats."""
function rough_dem(n, m; seed = 42, levels = 6, quantize = nothing)
    rng = Random.MersenneTwister(seed)
    dem = zeros(Float32, n, m)
    for k in 1:levels
        f = 2.0^k
        a = randn(rng, Float32, ceil(Int, n / f) + 2, ceil(Int, m / f) + 2)
        for i in 1:n, j in 1:m
            dem[i, j] += a[ceil(Int, i / f), ceil(Int, j / f)] / k
        end
    end
    isnothing(quantize) || (dem = round.(dem ./ quantize) .* quantize)
    dem
end

"""Directions alone, dropping the weights."""
directionsof(dem; method, cellsize = (1.0, 1.0)) =
    first(flowdirection(dem; method, cellsize))

const unitcell = (1.0, 1.0)

@testset "filldepressions" begin
    @testset "Zhou reaches the same surface as PriorityFlood" begin
        # Different queue strategies, same definition of a filled surface. Quantizing
        # makes large flats, where the two are most likely to diverge.
        for dem in (rough_dem(64, 64), rough_dem(80, 80; quantize = 0.25f0))
            @test filldepressions(dem; method = Zhou(), flats = false) ==
                  filldepressions(dem; method = PriorityFlood(), flats = false)
        end
    end

    @testset "the filled surface has no depressions left" begin
        filled = filldepressions(rough_dem(50, 50); method = Zhou(), flats = false)
        @test all(filled .>= rough_dem(50, 50))
        interior = CartesianIndices((2:49, 2:49))
        offsets = filter(!iszero, CartesianIndices((-1:1, -1:1)))
        @test all(any(filled[I + o] <= filled[I] for o in offsets) for I in interior)
    end

    @testset "flats are resolved by default" begin
        dem = rough_dem(40, 40; quantize = 0.5f0)
        @test filldepressions(dem) != filldepressions(dem; flats = false)
        # Resolving only nudges by float steps, so the surface is unchanged to Float32.
        @test filldepressions(dem) ≈ filldepressions(dem; flats = false)
    end
end

@testset "resolveflats" begin
    @testset "a tilted plane has no flats" begin
        dem = Float32[100.0f0 - 2i - 3j for i in 1:10, j in 1:10]
        @test resolveflats(dem) == dem
    end

    @testset "a flat drains towards its outlet" begin
        # A 2x2 flat at 5, draining at (3,3) towards the 4 at (4,4).
        dem = Float32[9 9 9 9; 9 5 5 9; 9 5 5 9; 9 9 9 4]
        z = resolveflats(dem)
        flat = dem .== 5
        @test z[.!flat] == dem[.!flat]        # only flat cells move
        @test minimum(z[flat]) == z[3, 3]     # the outlet cell is the lowest of the flat
    end

    @testset "every flat cell drains afterwards" begin
        z = filldepressions(rough_dem(60, 60; quantize = 0.5f0))
        interior = CartesianIndices((2:59, 2:59))
        offsets = filter(!iszero, CartesianIndices((-1:1, -1:1)))
        @test all(any(z[I + o] < z[I] for o in offsets) for I in interior)
    end
end

@testset "flowdirection" begin
    @testset "points downslope on a tilted plane" begin
        dem = Float32[100.0f0 - 2i - 3j for i in 1:8, j in 1:8]
        dirs = directionsof(dem; method = D8())
        @test all(
            CartesianIndex(dirs[I]) == CartesianIndex(1, 1) for
            I in CartesianIndices((2:7, 2:7))
        )
    end

    @testset "takes the steepest neighbour, not the lowest" begin
        # E drops 1.0 over 1.0; NE drops 1.2 over sqrt(2), so it is lower but shallower.
        dem = Float32[10 10 10; 10 10 9; 8.8 10 10]
        @test CartesianIndex(directionsof(dem; method = D8())[2, 2]) ==
              CartesianIndex(0, 1)
    end

    @testset "a cell with no lower neighbour is a pit" begin
        @test all(Geomorphometry.ispit, directionsof(fill(5.0f0, 5, 5); method = D8()))
    end

    @testset "receivers stay within what the method allows" begin
        dem = rough_dem(40, 40)
        for method in (D8(), DInf(), FD8())
            dirs = directionsof(dem; method)
            @test all(d -> Geomorphometry.ndirections(d) <= maxreceivers(method), dirs)
        end
        # A peak sends flow every way FD8 can, and DInf splits in two somewhere.
        peak = Float32[1 1 1; 1 2 1; 1 1 1]
        @test Geomorphometry.ndirections(directionsof(peak; method = FD8())[2, 2]) == 8
        @test any(d -> Geomorphometry.ndirections(d) == 2, directionsof(dem; method = DInf()))
    end

    @testset "all methods agree on which cells are pits" begin
        dem = rough_dem(40, 40; quantize = 0.5f0)
        pits(method) = Geomorphometry.ispit.(directionsof(dem; method))
        @test pits(D8()) == pits(DInf()) == pits(FD8())
    end

    @testset "weights sum to one wherever a cell drains" begin
        dem = rough_dem(40, 40)
        for method in (D8(), DInf(), FD8())
            dirs, weights = flowdirection(dem; method, cellsize = unitcell)
            for cell in CartesianIndices(dem)
                n = Geomorphometry.ndirections(dirs[cell])
                iszero(n) && continue
                @test sum(weights[cell, slot] for slot in 1:n) ≈ 1
            end
        end
    end

    @testset "weights are finite on a truncated neighbourhood" begin
        # At the border a cell is missing neighbours, leaving an angular gap of 180
        # degrees. Taken literally the contour length runs to tand(90) = Inf, and that
        # neighbour would take the cell's whole flow.
        _, weights = flowdirection(rough_dem(20, 20); method = FD8(), cellsize = unitcell)
        @test all(isfinite, weights)
    end

    @testset "a narrower weight type still routes the flow" begin
        dirs, weights = flowdirection(rough_dem(20, 20); method = FD8(), weighttype = Float16)
        @test eltype(weights) == Float16
        @test all(isfinite, weights)
    end
end
