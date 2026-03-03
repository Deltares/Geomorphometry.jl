using Geomorphometry
using Test

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
        B = skbr(rand(25, 25); iterations = 5)
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
        # DEM with a central depression
        dem = Float32[
            10 10 10 10 10
            10  8  7  8 10
            10  7  5  7 10
            10  8  7  8 10
            10 10 10 10 10
        ]

        @testset "depression_depth" begin
            bd = depression_depth(dem)
            # Center should be deepest (10 - 5 = 5)
            @test bd[3, 3] ≈ 5.0
            # Edges should be zero (not in depression)
            @test bd[1, 1] ≈ 0.0
            @test all(bd .>= 0)
        end

        @testset "depression_volume" begin
            vol = depression_volume(dem; cellsize=(1.0, 1.0))
            @test vol > 0
            # Volume should equal sum of depths for unit cells
            @test vol ≈ sum(depression_depth(dem))
        end

        @testset "drainage_potential" begin
            dp = drainage_potential(dem; cellsize=(1.0, 1.0))
            # Center (flat, high accumulation) should have low drainage
            @test dp[3, 3] ≈ 0.0
            # Edges should have higher drainage potential
            @test dp[1, 3] > dp[3, 3]
        end

        @testset "percentile_elevation" begin
            pct = percentile_elevation(dem; radius=1)
            # Center (lowest) should have low percentile
            @test pct[3, 3] ≈ 0.0
            # Corners (highest) should have high percentile
            @test pct[1, 1] > 0.5
        end

        @testset "flowaccumulation" begin
            acc, dir = flowaccumulation(dem; cellsize=(1.0, 1.0))
            # Center should have highest accumulation
            @test_broken acc[3, 3] == maximum(acc)
        end

        @testset "TWI and SPI" begin
            twi = TWI(dem; cellsize=(1.0, 1.0))
            spi = SPI(dem; cellsize=(1.0, 1.0))
            @test size(twi) == size(dem)
            @test size(spi) == size(dem)
        end
    end
end
