using Geomorphometry
using Test
using Aqua
using ExplicitImports

function Geomorphometry.decompose(
    ::Type{Int},
    direction::FlowDirection{D8D},
    center,
)
    map(decompose(direction)) do branch
        Int(branch) == 0 && return 0
        Int(branch) == 1 && return 1
        Int(branch) == 16 && return -1
        error("Unsupported test direction")
    end
end

@testset "Geomorphometry" begin
    @testset "Aqua" begin
        Aqua.test_all(Geomorphometry)
    end
    @testset "ExplicitImports" begin
        # `LocalFilters.Kernel` is a deliberate dependency on non-public API:
        # it carries the per-neighbour coefficients (1..9) that the terrain
        # update functions switch on, which no public LocalFilters neighbourhood
        # exposes. The public-ness checks use `Base.ispublic` only on Julia
        # 1.11+, so gate them by version to avoid false positives on 1.10.
        test_explicit_imports(
            Geomorphometry;
            all_explicit_imports_are_public = VERSION >= v"1.11",
            all_qualified_accesses_are_public = VERSION >= v"1.11" ?
                                                (; ignore = (:Kernel,)) : false,
        )
    end
    # Must run after ExplicitImports
    include("horizon.jl")
    @testset "pmf" begin
        # Write your own tests here.
        A = rand(25, 25)
        A[2, 2] = NaN
        B, flags = progressive_morphological_filter(A)
        @test (A .<= B) == (flags .== 0.0)

        Bc, flagsc = progressive_morphological_filter(A; circular = true)
        @test (A .<= Bc) == (flagsc .== 0.0)
        A = reshape(A, (size(A)..., 1))

        B3, flags3 = progressive_morphological_filter(A)
        @test length(size(B3)) == 2
        @test isnan.(B3) == isnan.(B)
        @test B3[end] == B[end]
        @test isnan.(flags3) == isnan.(flags)
    end
    @testset "smf" begin
        # Write your own tests here.
        A = rand(25, 25)
        A[2, 2] = NaN
        B = simple_morphological_filter(A)
    end
    @testset "pssm" begin
        B = pssm(rand(25, 25))
    end
    @testset "skb" begin
        B = skewness_balancing(rand(25, 25))
        B = skewness_balancing(rand(25, 25); mean = 0.25)
        B = Geomorphometry.skbr(rand(25, 25); iterations = 5)
    end
    @testset "pitremoval" begin
        B = pitremoval(rand(25, 25))
        B = pitremoval(rand(25, 25); limit = 0.1)
    end
    @testset "spread" begin
        points = [0.0 0 0 0 2; 0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0]
        initial = [8.0 8 8 8 4; 8 8 8 8 8; 8 8 8 8 8; 0 0 8 8 8; 0 0 8 8 8]
        friction = [1.0 200 1 1 1; 200 1 1 4 4; 1 1 4 4 4; 1 1 3 200 200; 1 Inf 3 200 4]
        @test spread(points, initial, friction; method = Tomlin()) ==
              spread(points, initial, friction; method = Eastman())
    end
    @testset "terrain" begin
        A = rand(25, 25)
        terrain_ruggedness_index(A)
        topographic_position_index(A)
        bathymetric_position_index(A)
        roughness_index_elevation(A)
        roughness(A)
        slope(A; method = Horn())
        slope(A; cellsize = (5, 5), method = ZevenbergenThorne())
        slope(A; cellsize = (10, 10), method = MaximumDownwardGradient())
        slope(A; cellsize = (10, 10), method = SymmetricGradient())
        aspect(A)
        aspect(A; method = MaximumDownwardGradient())
        laplacian(A; gis = true)
        hillshade(A)
    end
    @testset "neighbors clip to the domain" begin
        A = zeros(3, 3)
        @test Set(Geomorphometry.neighbors(A, CartesianIndex(1, 1))) ==
              Set([CartesianIndex(2, 1), CartesianIndex(1, 2), CartesianIndex(2, 2)])
    end
    @testset "hydrology" begin
        # Large depression with flow towards cell (1,1)
        dem = Float32[
            10 10 10 10 10
            10 8 7 8 10
            10 7 5 7 10
            10 8 7 8 10
            10 10 10 10 10
        ]

        rdem = height_above_nearest_drainage(dem)
        @test all(iszero, rdem)

        @testset "HAND path clipping" begin
            path_dem = Float32[4, 6, 5, 7, 0]
            output = zeros(Float32, length(path_dem))
            directions = [1, 1, 1, 1, 0]
            stream_mask = falses(length(path_dem))
            stream_mask[end] = true
            Geomorphometry._hand!(
                output,
                reverse(eachindex(path_dem)),
                directions,
                eachindex(path_dem),
                path_dem,
                stream_mask,
            )
            @test output == path_dem

            below_stream_dem = Float32[-2, 1, 0]
            fill!(output, 0)
            resize!(output, length(below_stream_dem))
            directions = [1, 1, 0]
            stream_mask = Bool[false, false, true]
            Geomorphometry._hand!(
                output,
                reverse(eachindex(below_stream_dem)),
                directions,
                eachindex(below_stream_dem),
                below_stream_dem,
                stream_mask,
            )
            @test output == Float32[0, 1, 0]
        end

        @testset "HAND flow directions" begin
            branch_dem = zeros(Float32, 3, 3)
            branch_dem[2, 2] = 10
            branch_dem[3, 2] = 6
            branch_dem[2, 3] = 8
            streams = falses(3, 3)
            streams[3, 2] = true
            streams[2, 3] = true
            order = [filter(!=(5), collect(eachindex(branch_dem))); 5]
            cells = CartesianIndices(branch_dem)
            output = similar(branch_dem)

            d8 = fill(FlowDirection{LDD, UInt8}(5), size(branch_dem))
            d8[2, 2] = FlowDirection{LDD, UInt8}(6)
            Geomorphometry._hand!(
                output,
                order,
                d8,
                cells,
                branch_dem,
                streams,
                (1, 1),
            )
            @test output[2, 2] == 4

            dinf = fill(FlowDirection{D8D, UInt8}(0), size(branch_dem))
            dinf[2, 2] = FlowDirection{D8D, UInt8}(1 | 64)
            Geomorphometry._hand!(
                output,
                order,
                dinf,
                cells,
                branch_dem,
                streams,
                (1, 1),
            )
            @test output[2, 2] == 2

            north_down = fill(FlowDirection{D8D, UInt8}(0), size(branch_dem))
            north_down[2, 2] = FlowDirection{D8D, UInt8}(64)
            north_stream = falses(3, 3)
            north_stream[2, 1] = true
            branch_dem[2, 1] = 7
            Geomorphometry._hand!(
                output,
                order,
                north_down,
                cells,
                branch_dem,
                north_stream,
                (1, -1),
            )
            @test output[2, 2] == 3
        end

        @testset "HAND stream burning and depression filling" begin
            depression_dem = reshape(Float32[0, 7, 5, 6, 4], 5, 1)
            cells = CartesianIndices(depression_dem)
            order = collect(eachindex(depression_dem))
            directions =
                fill(FlowDirection{D8D, UInt8}(16), size(depression_dem))
            directions[1] = FlowDirection{D8D, UInt8}(0)
            priority_directions =
                fill(CartesianIndex(-1, 0), size(depression_dem))
            priority_directions[1] = CartesianIndex(0, 0)

            streams = falses(size(depression_dem))
            streams[end] = true
            Geomorphometry._burn_streams!(
                streams,
                order,
                priority_directions,
                cells,
            )
            @test all(streams)

            filled = similar(depression_dem)
            Geomorphometry._fill_flow_depressions!(
                filled,
                order,
                directions,
                cells,
                depression_dem,
                (1, 1),
            )
            @test vec(filled) == Float32[0, 7, 7, 7, 7]

            streams .= false
            streams[1] = true
            output = similar(depression_dem)
            Geomorphometry._hand!(
                output,
                order,
                directions,
                cells,
                filled,
                streams,
                (1, 1),
            )
            @test vec(output) == Float32[0, 7, 7, 7, 7]

            generic_dem = Float32[4, 6, 5, 7, 0]
            generic_cells = eachindex(generic_dem)
            generic_order = reverse(eachindex(generic_dem))
            generic_directions =
                fill(FlowDirection{D8D, UInt8}(1), length(generic_dem))
            generic_directions[end] = FlowDirection{D8D, UInt8}(0)
            generic_priority = ones(Int, length(generic_dem))
            generic_priority[end] = 0

            generic_streams = falses(length(generic_dem))
            generic_streams[1] = true
            Geomorphometry._burn_streams!(
                generic_streams,
                generic_order,
                generic_priority,
                generic_cells,
            )
            @test all(generic_streams)

            generic_filled = similar(generic_dem)
            Geomorphometry._fill_flow_depressions!(
                generic_filled,
                generic_order,
                generic_directions,
                generic_cells,
                generic_dem,
                nothing,
            )
            @test generic_filled == Float32[7, 7, 7, 7, 0]
        end

        acc, dir = flowaccumulation(dem; method = D8())
        @test maximum(acc) == 10
        @test acc[2, 2] == 9
        @test acc[5, 5] == 1
        acc, dir = flowaccumulation(dem; method = DInf())
        @test maximum(acc) == 10
        @test acc[2, 2] == 9
        @test acc[5, 5] == 1
        acc, dir = flowaccumulation(dem; method = FD8(2))
        @test maximum(acc) == 10
        @test acc[2, 2] == 9
        @test acc[5, 5] == 1

        @testset "FD8 conserves flow" begin
            conservation_dem =
                Float32[mod(i + j + 6i * j, 101) for i in 1:8, j in 1:8]
            acc, _ =
                flowaccumulation(conservation_dem; method = FD8(), cellsize = (1, 1))
            @test maximum(acc) <= length(conservation_dem)
        end

        fdem = filldepressions(dem)
        @test all(==(10), fdem)

        @testset "depression_depth" begin
            bd = depression_depth(dem)
            # Center should be deepest (10 - 5 = 5)
            @test bd[3, 3] ≈ 5.0
            # Edges should be zero (not in depression)
            @test bd[1, 1] ≈ 0.0
            @test all(bd .>= 0)
        end

        @testset "depression_volume" begin
            vol = depression_volume(dem; cellsize = (1.0, 1.0))
            @test vol > 0
            # Volume should equal sum of depths for unit cells
            @test vol ≈ sum(depression_depth(dem))
        end

        @testset "drainage_potential" begin
            dp = drainage_potential(dem; cellsize = (1.0, 1.0))
            # Center (flat, high accumulation) should have low drainage
            @test dp[3, 3] ≈ 0.0
            # Edges should have higher drainage potential
            @test dp[1, 3] > dp[3, 3]
        end

        @testset "percentile_elevation" begin
            pct = percentile_elevation(dem; radius = 1)
            # Center (lowest) should have low percentile
            @test pct[3, 3] ≈ 0.0
            # Corners (highest) should have high percentile
            @test pct[1, 1] > 0.5
        end

        @testset "TWI and SPI" begin
            twi = topographic_wetness_index(dem; cellsize = (1.0, 1.0))
            spi = stream_power_index(dem; cellsize = (1.0, 1.0))
            @test size(twi) == size(dem)
            @test size(spi) == size(dem)
        end

        @testset "FlowDirection types and conventions" begin
            # LDD pit
            d = FlowDirection{LDD}(Int8(5))
            @test Geomorphometry.ispit(d)
            @test Geomorphometry.convention(d) == LDD

            # D8D single direction
            d = FlowDirection{D8D}(Int16(1))
            @test !Geomorphometry.ispit(d)
            @test Geomorphometry.issingle(d)
            @test Geomorphometry.ndirections(d) == 1
            @test Geomorphometry.convention(d) == D8D

            # D8D combined directions
            d = FlowDirection{D8D}(Int16(1 | 2 | 4))
            @test !Geomorphometry.issingle(d)
            @test Geomorphometry.ndirections(d) == 3
            @test Geomorphometry.length(Geomorphometry.decompose(d)) == 3

            # Convention traits
            @test Geomorphometry.ismulti(D8D)
            @test !Geomorphometry.ismulti(LDD)
        end

        @testset "CartesianIndex round-trip" begin
            # For each LDD direction, converting to CI and back should be identity
            for code in UInt8(1):UInt8(9)
                ci = only(decompose(
                    CartesianIndex{2},
                    FlowDirection{LDD}(code),
                ))
                code2 = FlowDirection{LDD}(ci).value
                @test code == code2
            end
            @test collect(decompose(
                CartesianIndex{2},
                FlowDirection{D8D}(1 | 64),
            )) == [CartesianIndex(1, 0), CartesianIndex(0, 1)]
            @test collect(decompose(
                CartesianIndex{2},
                FlowDirection{D8D}(1 | 64),
                CartesianIndex(2, 2),
            )) == [CartesianIndex(1, 0), CartesianIndex(0, 1)]
        end

        @testset "D8 tilted plane" begin
            # Simple tilted plane — no depressions, unambiguous flow
            dem = Float32[100.0f0 - 2.0f0 * i - 3.0f0 * j for i in 1:10, j in 1:10]
            acc, ldd = flowaccumulation(dem; method = D8())

            # The lowest corner (10,10) should have the highest accumulation
            inner = CartesianIndices((2:9, 2:9))
            @test all(ldd[ci] != 5 for ci in inner)  # no pits inside
        end

        @testset "FlowDirection show" begin
            # Pits render as the centre dot in both conventions
            @test repr(FlowDirection{LDD}(Int8(5))) == "·"
            @test repr(FlowDirection{D8D}(Int16(0))) == "·"

            # A single LDD direction renders as a single arrow character
            @test length(repr(FlowDirection{LDD}(Int8(9)))) == 1

            # D8D opposite pairs render as double-headed arrows
            @test repr(FlowDirection{D8D}(Int16(1 | 16))) == "↔"  # E+W
            @test repr(FlowDirection{D8D}(Int16(4 | 64))) == "↕"  # S+N
            @test repr(FlowDirection{D8D}(Int16(2 | 32))) == "⤢"  # SE+NW
            @test repr(FlowDirection{D8D}(Int16(8 | 128))) == "⤡"  # SW+NE

            # Non-opposite pairs render as a circular-average double-stroke arrow
            @test length(repr(FlowDirection{D8D}(Int16(1 | 2)))) == 1

            # Three or more directions render as the asterisk glyph
            @test repr(FlowDirection{D8D}(Int16(1 | 2 | 4))) == "✳"
        end

        @testset "FlowDirection equality and conversion" begin
            d = FlowDirection{LDD}(Int8(6))
            @test d == 6
            @test 6 == d
            @test d == big(6)
            @test big(6) == d
            # Directions in different conventions are never equal
            @test FlowDirection{LDD}(Int8(6)) != FlowDirection{D8D}(Int16(6))
            @test Int(d) == 6

            # Converting a combined D8D direction to a CartesianIndex is unsupported
            @test_throws "cannot be converted to a single CartesianIndex" CartesianIndex(
                FlowDirection{D8D}(Int16(3)),
            )
        end
    end

    @testset "curvature" begin
        # A tilted (linear) plane has zero curvature everywhere in its interior
        tilt = Float64[2.0 * i + 3.0 * j for i in 1:8, j in 1:8]
        @test size(profile_curvature(tilt)) == size(tilt)
        @test profile_curvature(tilt)[4, 4] ≈ 0.0 atol = 1e-9
        @test plan_curvature(tilt)[4, 4] ≈ 0.0 atol = 1e-9
        @test tangential_curvature(tilt)[4, 4] ≈ 0.0 atol = 1e-9

        # A curved (parabolic) surface has non-zero, finite curvature
        bowl = Float64[(i - 5.0)^2 + (j - 5.0)^2 for i in 1:9, j in 1:9]
        @test isfinite(profile_curvature(bowl)[3, 5])
        @test !iszero(profile_curvature(bowl)[3, 5])
        @test isfinite(plan_curvature(bowl)[3, 5])
        @test !iszero(plan_curvature(bowl)[3, 5])
        @test isfinite(tangential_curvature(bowl)[3, 5])
        @test !iszero(tangential_curvature(bowl)[3, 5])

        # The cellsize keyword scales the result
        @test profile_curvature(bowl; cellsize = (2.0, 2.0))[3, 5] !=
              profile_curvature(bowl)[3, 5]
    end

    @testset "aspect ZevenbergenThorne" begin
        tilt = Float64[2.0 * i + 3.0 * j for i in 1:8, j in 1:8]
        a = aspect(tilt; method = ZevenbergenThorne())
        @test size(a) == size(tilt)
        @test all(0 .<= a[2:7, 2:7] .< 360)
        # On a constant-gradient plane ZevenbergenThorne and Horn agree
        @test aspect(tilt; method = ZevenbergenThorne())[4, 4] ≈
              aspect(tilt; method = Horn())[4, 4] ≈ aspect(tilt; method = SymmetricGradient())[4, 4]
              
    end

    @testset "directional slope and laplacian" begin
        tilt = Float64[2.0 * i + 3.0 * j for i in 1:8, j in 1:8]

        # Directional slope differs from the omnidirectional slope
        @test slope(tilt; direction = 45.0)[4, 4] != slope(tilt)[4, 4]
        @test isfinite(slope(tilt; direction = 45.0)[4, 4])
        @test isfinite(
            slope(tilt; method = ZevenbergenThorne(), direction = 45.0)[4, 4],
        )

        # Maximum Downward Gradient does not support a direction
        @test_throws "Direction is not supported for MDG." slope(
            tilt;
            method = MDG(),
            direction = 10.0,
        )

        # Directional laplacian returns a finite field of matching size
        l = laplacian(tilt; direction = 45.0)
        @test size(l) == size(tilt)
        @test all(isfinite, l)
    end

    @testset "rugosity" begin
        # A flat surface has rugosity 1 (surface area equals planimetric area)
        flat = fill(5.0, 9, 9)
        @test rugosity(flat)[5, 5] ≈ 1.0

        # A rough surface has rugosity greater than 1
        rough = Float64[isodd(i + j) ? 0.0 : 10.0 for i in 1:9, j in 1:9]
        @test rugosity(rough)[5, 5] > 1.0
    end

    @testset "entropy" begin
        # A constant neighbourhood has zero Shannon entropy
        flat = fill(5.0, 9, 9)
        @test entropy(flat)[5, 5] ≈ 0.0 atol = 1e-7

        # A varied neighbourhood has positive entropy
        rough = Float64[isodd(i + j) ? 0.0 : 10.0 for i in 1:9, j in 1:9]
        @test entropy(rough)[5, 5] > 0.0
        # The step keyword controls binning; nothing disables it
        @test size(entropy(rough; step = nothing)) == size(rough)
    end

    @testset "terrain_ruggedness_index options" begin
        A = rand(25, 25)
        # Normalising divides by the neighbour count, yielding smaller values
        @test terrain_ruggedness_index(A; normalize = true)[13, 13] <
              terrain_ruggedness_index(A)[13, 13]
        # The non-squared (absolute difference) variant is supported
        @test all(isfinite, terrain_ruggedness_index(A; squared = false, normalize = true))
        # Disabling both squaring and normalisation warns the user
        @test_logs (:warn,) match_mode = :any terrain_ruggedness_index(
            A;
            squared = false,
            normalize = false,
        )
    end

    @testset "multihillshade" begin
        tilt = Float64[2.0 * i + 3.0 * j for i in 1:8, j in 1:8]
        mh = multihillshade(tilt)
        @test size(mh) == size(tilt)
        @test eltype(mh) == Union{Missing, UInt8}
        @test all(0 .<= skipmissing(mh) .<= 255)
    end

    @testset "skewness_balancing non-finite" begin
        B = rand(20, 20)
        B[5, 5] = Inf
        B[6, 6] = -Inf
        B[7, 7] = NaN
        mask = skewness_balancing(B)
        @test size(mask) == size(B)
        @test mask isa BitMatrix
    end
end
