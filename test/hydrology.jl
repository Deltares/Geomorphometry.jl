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

@testset "direction and weight grids follow the input type" begin
    # A direction eltype has no representation for a source missingval, and the weights
    # need a direction dimension that `similar` would otherwise drop. Both go through
    # hooks the wrapper extensions override.
    @test Geomorphometry._directiongrid(zeros(Float32, 3, 3), Geomorphometry.LDD) isa
          Matrix{<:Geomorphometry.FlowDirection}
    dirs, weights = flowdirection(rough_dem(8, 8); method = FD8(), cellsize = unitcell)
    @test size(weights) == (8, 8, maxreceivers(FD8()))
    @test all(>=(0), weights)
end

@testset "filldepressions leaves alone what it should" begin
    @testset "an integer DEM keeps its type and is not nudged" begin
        # Resolving flats moves cells by single float steps, which an integer cannot
        # represent, so it is off by default rather than throwing.
        dem = fill(5, 5, 5)
        @test filldepressions(dem) == dem
        @test eltype(filldepressions(dem)) == Int
    end

    @testset "cells the caller excluded are untouched" begin
        dem = fill(5.0, 5, 5)
        @test filldepressions(dem, trues(5, 5)) == dem
        # Excluding just the interior leaves those cells alone while the rim resolves.
        mask = falses(5, 5)
        mask[2:4, 2:4] .= true
        z = filldepressions(dem, mask)
        @test z[2:4, 2:4] == dem[2:4, 2:4]
    end

    @testset "a flat cannot be raised over terrain that was higher" begin
        # One float step of relief: the centre is above its neighbours by the smallest
        # amount representable, so raising the flat around it must not overtake it.
        dem = fill(1.0, 5, 5)
        dem[3, 3] = nextfloat(1.0)
        z = filldepressions(dem)
        centre = CartesianIndex(3, 3)
        offsets = filter(!iszero, CartesianIndices((-1:1, -1:1)))
        @test any(z[centre + o] < z[centre] for o in offsets)
    end
end

@testset "FD8 weights agree between the stencil and accumulation" begin
    # Both read their raster contour lengths from `_rastercontour`, so the share a cell
    # gives each neighbour must not depend on which route derived it.
    peak = Float64[1 1 1; 1 2 1; 1 1 1]
    centre = CartesianIndex(2, 2)
    method = FD8(; p = 1.0)

    directions, weights = flowdirection(peak; method, cellsize = unitcell)
    acc, = flowaccumulation(copy(peak); method, cellsize = unitcell)

    offsets = Geomorphometry.decompose(CartesianIndex{2}, directions[centre])
    for (slot, offset) in enumerate(offsets)
        iszero(offset) && continue
        # The centre starts with one cell of area, so a neighbour's accumulation above its
        # own area is exactly the share the centre sent it.
        @test weights[centre, slot] ≈ acc[centre + offset] - 1 rtol = 1e-5
    end
end

@testset "carved directions are steepest descent" begin
    @testset "on a drained surface the carve matches the stencil" begin
        # The flood pops in elevation order, so every strictly lower neighbour of a cell
        # is already popped when that cell is. Choosing the steepest of those is D8, so
        # on a surface that drains the two routes agree away from the rim, where they
        # differ by convention rather than by rule.
        z = filldepressions(rough_dem(80, 80))
        _, carved = flowaccumulation(copy(z); method = D8(), cellsize = unitcell)
        stencil = directionsof(z; method = D8())
        @test all(carved[I] == stencil[I] for I in CartesianIndices((2:79, 2:79)))
    end

    @testset "prefers the steeper orthogonal over the lower diagonal" begin
        # NE is lower, E is steeper; the flood reaches (3,2) from NE either way.
        dem = Float32[
            20 20 20 20 20
            20 20 20 20 20
            20 10 9.0 20 20
            20 8.8 20 20 20
            20 20 20 20 20
        ]
        z = filldepressions(dem)
        _, carved = flowaccumulation(copy(z); method = D8(), cellsize = unitcell)
        @test carved[3, 2] == directionsof(z; method = D8())[3, 2]
    end
end

@testset "DInf facets" begin
    @testset "the steepest facet need not touch the steepest edge" begin
        # A facet's slope is that of the plane through its three cells, which can exceed
        # both of its own edges. So searching only the facets beside the steepest single
        # neighbour is not Tarboton's method, and here it would pick the wrong quadrant.
        dem = zeros(Float64, 3, 3)
        centre = CartesianIndex(2, 2)
        dem[CartesianIndex(3, 2)] = -1.0      # E,  edge slope 1.0
        dem[CartesianIndex(3, 3)] = -1.4      # NE, edge slope 0.990
        dem[CartesianIndex(1, 2)] = -1.05     # W,  edge slope 1.05, the steepest edge
        near, far, _, _ = Geomorphometry._dinfreceivers(dem, centre, nothing, unitcell)
        @test Set((near - centre, far - centre)) ==
              Set((CartesianIndex(1, 0), CartesianIndex(1, 1)))
    end

    @testset "flow is split between two neighbours sharing a corner" begin
        dem = rough_dem(40, 40)
        for cell in CartesianIndices(dem)
            receivers = Geomorphometry._dinfreceivers(dem, cell, nothing, unitcell)
            isnothing(receivers) && continue
            near, far, afraction, bfraction = receivers
            @test afraction >= 0 && bfraction >= 0 && afraction + bfraction ≈ 1
            @test maximum(abs, Tuple(far - near)) == 1      # adjacent to each other
            @test maximum(abs, Tuple(near - cell)) == 1     # and both to the cell
        end
    end

    @testset "a planar slope sends all flow one way" begin
        dem = Float64[10.0 - i for i in 1:5, j in 1:5]
        near, _, afraction, bfraction =
            Geomorphometry._dinfreceivers(dem, CartesianIndex(3, 3), nothing, unitcell)
        @test near - CartesianIndex(3, 3) == CartesianIndex(1, 0)
        @test afraction ≈ 1 && bfraction ≈ 0
    end
end

@testset "DInf falls back to bracketing an aspect off a raster" begin
    # Rasters take the facet route above. Grids with no facet geometry derive an aspect
    # per cell and split between the two neighbours bracketing it, which is the only
    # path a DGGS uses, so it is exercised directly here.
    using Geomorphometry: _localaspect, _flowneighbors

    @testset "a plane sends everything one way" begin
        # Drops along dim1, which is east, so the aspect is due east and the two
        # bracketing neighbours collapse onto one.
        plane = Float64[10.0 - i for i in 1:5, j in 1:5]
        cell = CartesianIndex(3, 3)
        aspect = _localaspect(plane, cell, unitcell)
        @test aspect ≈ 90
        near, far, afraction, bfraction = _flowneighbors(plane, cell, aspect, unitcell)
        @test near == far == cell + CartesianIndex(1, 0)
        @test afraction ≈ 1 && bfraction ≈ 0
    end

    @testset "an aspect between two bearings is split between them" begin
        tilt = Float64[10.0 - 2i - j for i in 1:5, j in 1:5]
        cell = CartesianIndex(3, 3)
        aspect = _localaspect(tilt, cell, unitcell)
        near, far, afraction, bfraction = _flowneighbors(tilt, cell, aspect, unitcell)
        @test 45 < aspect < 90                      # between NE and E
        @test Set((near - cell, far - cell)) ==
              Set((CartesianIndex(1, 1), CartesianIndex(1, 0)))
        @test afraction + bfraction ≈ 1
        @test afraction > 0 && bfraction > 0
    end

    @testset "a flat has no aspect, so the cell has no receivers" begin
        # NaN here is what `_dinfreceivers` guards: bracketing it would pick the cell
        # itself as both of its receivers.
        @test isnan(_localaspect(fill(5.0, 5, 5), CartesianIndex(3, 3), unitcell))
    end
end

@testset "FD8 contour lengths" begin
    # A peak whose eight neighbours all lie one below it, so every way down is equally
    # steep and only the cell geometry decides the split: slope carries the distance
    # travelled, the contour length the width crossed.
    peak = fill(0.0f0, 5, 5)
    peak[2:4, 2:4] .= 9
    peak[3, 3] = 10

    function peakshares(cellsize)
        acc, = flowaccumulation(peak; method = FD8(), cellsize)
        area = abs(cellsize[1] * cellsize[2])
        share(I) = (acc[I] - area) / area
        (
            row = share(CartesianIndex(2, 3)) + share(CartesianIndex(4, 3)),
            col = share(CartesianIndex(3, 2)) + share(CartesianIndex(3, 4)),
            diagonal = sum(share(I) for I in CartesianIndex.(((2, 2), (2, 4), (4, 2), (4, 4)))),
        )
    end

    @testset "square cells reproduce Quinn's L1 = 0.5, L2 = 0.354" begin
        s = peakshares(unitcell)
        @test s.row + s.col + s.diagonal ≈ 1
        @test s.row ≈ s.col
        # Each cardinal outweighs each diagonal by 2^p before normalising.
        @test s.row / 2 ≈ (s.diagonal / 4) * 2^FD8().p rtol = 1e-5
    end

    @testset "a long thin cell sends more water out of its long side" begin
        # Three times as tall as wide: a step along the columns is the short one and
        # leaves through the tall face, so it dominates twice over. Invisible on squares.
        tall, wide = peakshares((3.0, 1.0)), peakshares((1.0, 3.0))
        @test tall.col > tall.row
        @test tall.col ≈ wide.row && tall.row ≈ wide.col   # mirror images
    end

    @testset "the split depends on the cell's shape, not its size" begin
        # Slope loses exactly what the contour length gains, so scaling changes nothing.
        @test all(isapprox.(values(peakshares((30.0, 10.0))), values(peakshares((3.0, 1.0))), rtol = 1e-5))
    end
end

@testset "height_above_nearest_drainage over a given stream mask" begin
    z = filldepressions(rough_dem(30, 30))
    mask = falses(30, 30)
    mask[1, 1] = true
    hand = height_above_nearest_drainage(z, mask)
    @test size(hand) == size(z)
    @test all(h -> !isfinite(h) || h >= 0, hand)
    @test hand[1, 1] == 0
end

@testset "the carve drains into -Inf" begin
    # The manual sets sea to -Inf before filling, so -Inf has to be a usable receiver: it
    # is lower than anything. NaN and +Inf are not, and must not be routed into.
    dem = fill(0.0, 3, 3)
    dem[1, :] .= -Inf
    dem[3, :] .= -Inf
    dem[:, 1] .= -Inf
    dem[:, 3] .= -Inf
    dem[2, 2] = 5.0
    _, carved = flowaccumulation(copy(dem); method = D8(), cellsize = unitcell)
    @test !Geomorphometry.ispit(carved[2, 2])

    posinf = copy(dem)
    posinf[1, 2] = Inf
    _, uphill = flowaccumulation(copy(posinf); method = D8(), cellsize = unitcell)
    @test CartesianIndex(uphill[2, 2]) != CartesianIndex(-1, 0)
end

@testset "flowaccumulation over a given direction field" begin
    @testset "every cell counts once" begin
        z = filldepressions(rough_dem(40, 40))
        for method in (D8(), DInf(), FD8())
            dirs, weights = flowdirection(z; method, cellsize = unitcell)
            acc = flowaccumulation(z, dirs, weights; cellsize = unitcell)
            # Each cell contributes its own area, so nothing is created or lost.
            @test minimum(acc) >= 1
            @test maximum(acc) <= length(z)
        end
    end

    @testset "a single chain accumulates linearly" begin
        # A staircase down one column: each cell passes its own area plus all above it.
        dem = Float32[10.0f0 - i for i in 1:6, j in 1:1]
        dirs, weights = flowdirection(dem; method = D8(), cellsize = unitcell)
        acc = flowaccumulation(dem, dirs, weights; cellsize = unitcell)
        @test vec(acc) == Float32[1, 2, 3, 4, 5, 6]
    end

    @testset "it follows the directions it is given, not the terrain" begin
        # A flat DEM has no directions of its own, so these are purely the ones given:
        # send everything east along the row.
        dem = zeros(Float32, 1, 4)
        dirs, _ = flowdirection(dem; method = D8(), cellsize = unitcell)
        dirs[1, 1:3] .= FlowDirection{LDD}(CartesianIndex(0, 1))
        acc = flowaccumulation(dem, dirs; cellsize = unitcell)
        @test vec(acc) == Float32[1, 2, 3, 4]

        # Cutting a cell off makes it a pit: it keeps its own flow and passes none on.
        dirs[1, 2] = FlowDirection{LDD}(CartesianIndex(0, 0))
        @test vec(flowaccumulation(dem, dirs; cellsize = unitcell)) == Float32[1, 2, 1, 2]
    end

    @testset "the methods route the same water off the grid" begin
        # Divergent methods spread flow differently, but the total leaving must match.
        z = filldepressions(rough_dem(50, 50))
        R = CartesianIndices(z)
        function leaving(method)
            dirs, weights = flowdirection(z; method, cellsize = unitcell)
            acc = flowaccumulation(z, dirs, weights; cellsize = unitcell)
            sum(
                acc[cell] for cell in R if
                Geomorphometry.ispit(dirs[cell]) || any(
                    !(cell + Geomorphometry._orient(o, unitcell) in R) for
                    o in Geomorphometry.decompose(CartesianIndex{2}, dirs[cell])
                )
            )
        end
        @test leaving(D8()) ≈ length(z) rtol = 1e-4
        @test leaving(DInf()) ≈ leaving(D8()) rtol = 1e-3
        @test leaving(FD8()) ≈ leaving(D8()) rtol = 1e-3
    end
end

@testset "accumulation follows the facet, not the carved D8 direction" begin
    @testset "a facet that misses the carved direction still routes" begin
        # The window from "the steepest facet need not touch the steepest edge", with a
        # rim to drain to. The carve takes W, the steepest single edge; the steepest facet
        # is E/NE. Both of those receivers are earlier in the pop order than the centre,
        # so both are legal and the facet decides.
        dem = fill(0.5, 5, 5)
        dem[1, :] .= -10
        dem[5, :] .= -10
        dem[:, 1] .= -10
        dem[:, 5] .= -10
        dem[3, 3] = 0.0
        dem[4, 3] = -1.0      # E,  steepest facet with NE
        dem[4, 4] = -1.4      # NE
        dem[2, 3] = -1.05     # W,  steepest single edge, so the carved direction
        acc, = flowaccumulation(copy(dem); method = DInf(), cellsize = unitcell)
        @test acc[4, 3] ≈ 1.5155 atol = 1e-4
        @test acc[4, 4] ≈ 1.4845 atol = 1e-4
        @test acc[2, 3] ≈ 1.0                    # nothing went to the carved direction
    end

    @testset "a spill cell whose facet looks back into its depression" begin
        # Inside a depression the facet points at the floor, which the flood popped after
        # the spill cell, so neither receiver can take the flow. The carved escape route
        # is the only edge guaranteed earlier in the pop order, and taking it is what
        # keeps the single pass sound rather than throwing.
        dem = rough_dem(120, 120)
        acc, = flowaccumulation(copy(dem); method = DInf(), cellsize = unitcell)
        @test all(isfinite, acc)
        @test minimum(acc) >= 1
        @test maximum(acc) <= length(dem)
    end

    @testset "the carve wins where a facet cannot reach nodata" begin
        # No facet can name a -Inf cell, since the vertices must be finite, but the carve
        # can: a coastal cell drains straight to sea rather than along the shore.
        dem = fill(1.0, 4, 4)
        dem[:, 1] .= -Inf
        dem[2, 2] = 0.5
        acc, dirs = flowaccumulation(copy(dem); method = DInf(), cellsize = unitcell)
        @test all(isfinite, acc)
        @test !Geomorphometry.ispit(dirs[2, 2])
        # Straight west into the sea, not north or south along the shore.
        @test CartesianIndex(dirs[2, 2]) == CartesianIndex(0, -1)
    end
end

@testset "watersheds" begin
    @testset "one basin when everything drains to one point" begin
        # A cone draining to its single lowest cell on the rim.
        dem = Float32[hypot(i - 1, j - 1) for i in 1:9, j in 1:9]
        for method in (D8(), DInf(), FD8())
            dirs, weights = flowdirection(dem; method, cellsize = unitcell)
            @test length(unique(watersheds(dem, dirs, weights; cellsize = unitcell))) == 1
        end
    end

    @testset "a label is exactly the set of cells sharing an outlet" begin
        z = filldepressions(rough_dem(40, 40))
        dirs, weights = flowdirection(z; method = D8(), cellsize = unitcell)
        labels = watersheds(z, dirs, weights; cellsize = unitcell)
        @test all(>(0), labels)
        # Every cell carries the label of the cell it drains to.
        for cell in CartesianIndices(z)
            offsets = Geomorphometry.decompose(CartesianIndex{2}, dirs[cell])
            for o in offsets
                iszero(o) && continue
                target = cell + Geomorphometry._orient(o, unitcell)
                target in CartesianIndices(z) || continue
                @test labels[target] == labels[cell]
            end
        end
    end

    @testset "a divergent cell joins the basin taking most of its flow" begin
        z = filldepressions(rough_dem(30, 30))
        dirs, weights = flowdirection(z; method = FD8(), cellsize = unitcell)
        labels = watersheds(z, dirs, weights; cellsize = unitcell)
        for cell in CartesianIndices(z)
            main = Geomorphometry._mainreceiver(
                dirs, weights, cell, unitcell, CartesianIndices(z)
            )
            isnothing(main) && continue
            @test labels[main] == labels[cell]
        end
    end

    @testset "flooding labels an unfilled DEM, following directions does not" begin
        # A hollow with no outlet. Flooding routes through it, so it joins the basin that
        # fills it; the directions stop dead inside it, making it a basin of its own.
        dem = rough_dem(40, 40)
        flooded = watersheds(dem)
        dirs, weights = flowdirection(dem; method = D8(), cellsize = unitcell)
        followed = watersheds(dem, dirs, weights; cellsize = unitcell)
        @test all(>(0), flooded)
        @test length(unique(flooded)) < length(unique(followed))
    end

    @testset "it does not touch the elevation model" begin
        dem = rough_dem(20, 20)
        before = copy(dem)
        watersheds(dem)
        @test dem == before
    end
end

@testset "a minority branch cannot suppress the labels" begin
    # A -> B, with B sending 10% back to A and 90% on to the pit C. The chosen paths
    # A -> B -> C are acyclic, but the 10% back edge is a cycle in the full receiver
    # graph, which once left every cell unlabelled.
    dem = zeros(Float32, 1, 3)
    east = Geomorphometry._directioncode(dem, CartesianIndex(0, 1), unitcell)
    west = Geomorphometry._directioncode(dem, CartesianIndex(0, -1), unitcell)

    directions = fill(FlowDirection{Geomorphometry.D8D}(zero(UInt8)), 1, 3)
    weights = zeros(Float32, 1, 3, 8)
    directions[1, 1] = FlowDirection{Geomorphometry.D8D}(east)
    weights[1, 1, 1] = 1
    directions[1, 2] = FlowDirection{Geomorphometry.D8D}(west | east)
    for (slot, bit) in enumerate(sort([west, east]))
        weights[1, 2, slot] = bit == west ? 0.1 : 0.9
    end

    labels = watersheds(dem, directions, weights; cellsize = unitcell)
    @test all(>(0), labels)
    @test length(unique(labels)) == 1          # all three drain to the same pit
end
