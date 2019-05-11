using GeoRasterFiltering
using Test

@testset "GeoRasterFiltering.jl" begin
    # Write your own tests here.
    A = rand(25, 25)
    A[2,2] = NaN
    B, flags = pmf(A)
    @test (A .<= B) == (flags .== 0.)
end
