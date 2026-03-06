using Geomorphometry
using Test

include("horizon.jl")

@testset "Geomorphometry" begin
    @testset "pmf" begin
        # Write your own tests here.
        A = rand(25, 25)
        A[2, 2] = NaN
        B, flags = pmf(A)
        @test (A .<= B) == (flags .== 0.0)

        Bc, flagsc = pmf(A; circular = true)
        @test (A .<= Bc) == (flagsc .== 0.0)
        A = reshape(A, (size(A)..., 1))

        B3, flags3 = pmf(A)
        @test length(size(B3)) == 2
        @test isnan.(B3) == isnan.(B)
        @test B3[end] == B[end]
        @test isnan.(flags3) == isnan.(flags)
    end
    @testset "smf" begin
        # Write your own tests here.
        A = rand(25, 25)
        A[2, 2] = NaN
        B = smf(A)
    end
    @testset "pssm" begin
        B = pssm(rand(25, 25))
    end
    @testset "skb" begin
        B = skb(rand(25, 25))
        B = skb(rand(25, 25); mean = 0.25)
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
        TRI(A)
        TPI(A)
        BPI(A)
        RIE(A)
        roughness(A)
        slope(A; method = Horn())
        slope(A; cellsize = (5, 5), method = ZevenbergenThorne())
        slope(A; cellsize = (10, 10), method = MDG())
        aspect(A)
        aspect(A; method = MDG())
        curvature(A)
        hillshade(A)
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
        @test rdem[1, 1] == 0  # all descending
        @test rdem[4, 4] == 3  # highest asecnding point
        @test maximum(rdem) == 3

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
            twi = TWI(dem; cellsize = (1.0, 1.0))
            spi = SPI(dem; cellsize = (1.0, 1.0))
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
                ci = CartesianIndex(FlowDirection{LDD}(code))
                code2 = FlowDirection{LDD}(ci).value
                @test code == code2
            end
        end

        @testset "D8 tilted plane" begin
            # Simple tilted plane — no depressions, unambiguous flow
            dem = Float32[100.0f0 - 2.0f0 * i - 3.0f0 * j for i in 1:10, j in 1:10]
            acc, ldd = flowaccumulation(dem; method = D8())

            # The lowest corner (10,10) should have the highest accumulation
            inner = CartesianIndices((2:9, 2:9))
            @test all(ldd.data[ci] != 5 for ci in inner)  # no pits inside
        end
    end
end
