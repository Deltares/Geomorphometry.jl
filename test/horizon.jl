using Geomorphometry
using Test

@testset "horizon" begin
    @testset "GridSweep types" begin
        @test GridSweep().directions == 8
        @test GridSweep(directions=4).directions == 4
        @test GridSweep(directions=16).directions == 16
    end

    @testset "GridSweep 4 directions" begin
        flat = ones(10, 10) * 100.0
        h = horizon_angle(flat; method=GridSweep(directions=4), cellsize=(1.0, 1.0))

        @test size(h) == (10, 10, 4)
        # Interior cells on flat terrain should have 0° horizons
        # Order: N, E, S, W
        @test h[5, 5, 1] ≈ 0.0 atol=1e-5  # N
        @test h[5, 5, 2] ≈ 0.0 atol=1e-5  # E
        @test h[5, 5, 3] ≈ 0.0 atol=1e-5  # S
        @test h[5, 5, 4] ≈ 0.0 atol=1e-5  # W
    end

    @testset "GridSweep 8 directions" begin
        flat = ones(10, 10) * 100.0
        h = horizon_angle(flat; method=GridSweep(directions=8), cellsize=(1.0, 1.0))

        @test size(h) == (10, 10, 8)
        # Interior cells on flat terrain should have 0° horizons
        # Order: N, NE, E, SE, S, SW, W, NW
        @test h[5, 5, 1] ≈ 0.0 atol=1e-5  # N
        @test h[5, 5, 2] ≈ 0.0 atol=1e-5  # NE
        @test h[5, 5, 4] ≈ 0.0 atol=1e-5  # SE
    end

    @testset "GridSweep invalid directions" begin
        flat = ones(5, 5)
        @test_throws ArgumentError horizon_angle(flat; method=GridSweep(directions=6))
    end

    @testset "sloping terrain" begin
        # Terrain that increases with row index (slopes south)
        slope_dem = [Float64(row) for row in 1:10, col in 1:10]
        h = horizon_angle(slope_dem; method=GridSweep(directions=4), cellsize=(1.0, 1.0))

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

    @testset "GridSweep 16 directions" begin
        dem = [Float64(row + sin(col/3)) for row in 1:32, col in 1:32]

        h16 = horizon_angle(dem; method=GridSweep(directions=16), cellsize=(1.0, 1.0))
        @test size(h16) == (32, 32, 16)

        # All values should be in valid range
        @test all(-90 .<= h16 .<= 90)
    end

    @testset "maxdist parameter" begin
        # Create terrain with a distant peak
        dem = ones(20, 20) * 100.0
        dem[1, 10] = 200.0  # Peak at north edge

        # With unlimited distance, should see the peak
        h_inf = horizon_angle(dem; method=GridSweep(directions=4), cellsize=(1.0, 1.0), maxdist=Inf)

        # With limited distance, may not see the peak from far cells
        h_short = horizon_angle(dem; method=GridSweep(directions=4), cellsize=(1.0, 1.0), maxdist=5.0)

        # Cell far from peak should have different N horizon with limited distance
        # N is direction 1
        @test h_inf[15, 10, 1] >= h_short[15, 10, 1]
    end
end
