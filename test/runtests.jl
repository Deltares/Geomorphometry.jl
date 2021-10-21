using GeoRasterFiltering
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
    B, flags = smf(A)
    @test (A .<= B) == (flags .== 0.)
end
