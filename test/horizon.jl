using Geomorphometry
using Test
using Rasters
using GeoArrays
using Rasters.DimensionalData
using Rasters.DimensionalData.Lookups

@testset "horizon" begin
    @testset "4 directions" begin
        flat = ones(10, 10) * 100.0
        h = horizon_angle(flat; directions=4, cellsize=(1.0, 1.0))

        @test size(h) == (10, 10, 4)
        # Interior cells on flat terrain should have 0° horizons
        # Order: N, E, S, W
        @test h[5, 5, 1] ≈ 0.0 atol=1e-5  # N
        @test h[5, 5, 2] ≈ 0.0 atol=1e-5  # E
        @test h[5, 5, 3] ≈ 0.0 atol=1e-5  # S
        @test h[5, 5, 4] ≈ 0.0 atol=1e-5  # W
    end

    @testset "8 directions" begin
        flat = ones(10, 10) * 100.0
        h = horizon_angle(flat; directions=8, cellsize=(1.0, 1.0))

        @test size(h) == (10, 10, 8)
        # Interior cells on flat terrain should have 0° horizons
        # Order: N, NE, E, SE, S, SW, W, NW
        @test h[5, 5, 1] ≈ 0.0 atol=1e-5  # N
        @test h[5, 5, 2] ≈ 0.0 atol=1e-5  # NE
        @test h[5, 5, 4] ≈ 0.0 atol=1e-5  # SE
    end

    @testset "16 directions (needing rotation)" begin
        dem = [Float64(row + sin(col/3)) for row in 1:32, col in 1:32]

        h16 = horizon_angle(dem; directions=16, cellsize=(1.0, 1.0))
        @test size(h16) == (32, 32, 16)

        # All values should be in valid range
        @test all(-90 .<= h16 .<= 90)
    end

    @testset "invalid directions" begin
        flat = ones(5, 5)
        @test_throws ArgumentError horizon_angle(flat; directions=6)
    end

    @testset "sloping terrain" begin
        # Terrain that increases with row index (slopes south)
        slope_dem = [Float64(row) for row in 1:10, col in 1:10]
        h = horizon_angle(slope_dem; directions=4, cellsize=(1.0, 1.0))

        # Order: N, E, S, W
        # Looking south (uphill), horizon should be positive
        @test h[5, 5, 3] > 0  # S
        # Looking north (downhill), horizon should be negative or zero
        @test h[5, 5, 1] <= 0  # N
        # E-W should be flat (0°)
        @test h[5, 5, 2] ≈ 0.0 atol=1e-5  # E
        @test h[5, 5, 4] ≈ 0.0 atol=1e-5  # W
    end

    @testset "sky_view_factor" begin
        # Flat terrain should have SVF = 1.0
        flat = ones(10, 10) * 100.0
        svf = sky_view_factor(flat; cellsize=(1.0, 1.0))
        @test svf[5, 5] ≈ 1.0 atol=1e-5

        # Sloping terrain should have SVF < 1.0
        slope_dem = [Float64(row) for row in 1:10, col in 1:10]
        svf_slope = sky_view_factor(slope_dem; cellsize=(1.0, 1.0))
        @test svf_slope[5, 5] < 1.0
        @test svf_slope[5, 5] > 0.0
    end

    @testset "Raster wrapping is preserved" begin
        A = Float64[row + sin(col/3) for row in 1:16, col in 1:16]
        x = Rasters.X(Sampled(1.0:1.0:16.0; sampling=Intervals(Start()),
                              order=ForwardOrdered(), span=Regular(1.0)))
        y = Rasters.Y(Sampled(1.0:1.0:16.0; sampling=Intervals(Start()),
                              order=ForwardOrdered(), span=Regular(1.0)))
        r = Raster(A, (x, y); crs=Rasters.EPSG(32633))

        h8 = horizon_angle(r; directions=8)
        @test h8 isa Raster
        @test size(h8) == (16, 16, 8)
        @test Rasters.dims(h8)[3] isa Rasters.Band

        h16 = horizon_angle(r; directions=16)
        @test h16 isa Raster
        @test size(h16) == (16, 16, 16)

        svf = sky_view_factor(r; directions=8)
        @test svf isa Raster
        @test size(svf) == (16, 16)
    end

    @testset "GeoArray wrapping is preserved" begin
        A = Float64[row + sin(col/3) for row in 1:16, col in 1:16]
        ga = GeoArray(A)
        GeoArrays.bbox!(ga, (min_x=0.0, min_y=0.0, max_x=16.0, max_y=16.0))

        h8 = horizon_angle(ga; directions=8)
        @test h8 isa GeoArray
        @test size(h8) == (16, 16, 8)

        h16 = horizon_angle(ga; directions=16)
        @test h16 isa GeoArray
        @test size(h16) == (16, 16, 16)

        svf = sky_view_factor(ga; directions=8)
        @test svf isa GeoArray
        @test size(svf) == (16, 16)
    end
end
