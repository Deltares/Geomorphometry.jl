using GeoArrayOps
using Test

@testset "GeoArrayOps" begin
    @testset "pmf" begin
        # Write your own tests here.
        A = rand(25, 25)
        A[2, 2] = NaN
        B, flags = pmf(A)
        @test (A .<= B) == (flags .== 0.0)

        Bc, flagsc = pmf(A; circular=true)
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
    @testset "psf" begin
        B = psf(rand(25, 25))
    end
    @testset "skb" begin
        B = skb(rand(25, 25))
        B = skb(rand(25, 25); mean=0.25)
        B = skbr(rand(25, 25); iterations=5)
    end
    @testset "pitremoval" begin
        B = pitremoval(rand(25, 25))
        B = pitremoval(rand(25, 25); limit=0.1)
    end
    @testset "spread" begin
        points = [0.0 0 0 0 2; 0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0]
        initial = [8.0 8 8 8 4; 8 8 8 8 8; 8 8 8 8 8; 0 0 8 8 8; 0 0 8 8 8]
        friction = [1.0 200 1 1 1; 200 1 1 4 4; 1 1 4 4 4; 1 1 3 200 200; 1 Inf 3 200 4]
        @test spread(points, initial, friction) == spread2(points, initial, friction)
    end
    @testset "terrain" begin
        A = rand(25, 25)
        TRI(A)
        TPI(A)
        roughness(A)
        slope(A; method=Horn())
        slope(A; cellsize=5, method=ZevenbergenThorne())
        slope(A; cellsize=10, method=MDG())
        aspect(A)
        aspect(A, method=MDG())
        curvature(A)
        hillshade(A)
    end
end
