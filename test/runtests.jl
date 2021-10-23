using GeoArrayOps
using Test

@testset "pmf" begin
    # Write your own tests here.
    A = rand(25, 25)
    A[2,2] = NaN
    B, flags = pmf(A)
    @test (A .<= B) == (flags .== 0.)
end
@testset "smf" begin
    # Write your own tests here.
    A = rand(25, 25)
    A[2,2] = NaN
    B = smf(A)
end
@testset "pssm" begin
    B = pssm(rand(25, 25))
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
end
