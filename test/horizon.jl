using Geomorphometry
using Test

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
        # Multiples of 8 that aren't powers of 2 (e.g. 24, 40) must also reject.
        @test_throws ArgumentError horizon_angle(flat; directions=24)
        @test_throws ArgumentError horizon_angle(flat; directions=40)
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

    @testset "distant peak beyond flat foreground" begin
        # 10 m peak at the north edge, 9 cells of flat ground in front.
        dem = fill(100.0, 10, 5)
        dem[1, :] .= 110.0
        h = horizon_angle(dem; directions = 4, cellsize = (1.0, 1.0))
        # N, E, S, W. From (10, 3) the peak is 9 cells north.
        @test h[10, 3, 1] ≈ atand(10 / 9) atol = 1e-3
        @test h[10, 3, 2] ≈ 0.0 atol = 1e-5
        @test h[10, 3, 4] ≈ 0.0 atol = 1e-5
    end

    @testset "angle to a spike decays with distance" begin
        # 10 m spike at row 2; observers south of it must see it shrink.
        dem = fill(100.0, 10, 5)
        dem[2, :] .= 110.0
        h = horizon_angle(dem; directions = 4, cellsize = (1.0, 1.0))
        @test h[5, 3, 1] ≈ atand(10 / 3) atol = 1e-3
        @test h[10, 3, 1] ≈ atand(10 / 8) atol = 1e-3
        @test h[10, 3, 1] < h[5, 3, 1]
    end

    @testset "missing cells use missing_elevation (default 0)" begin
        # 200 m peak at row 1, NaN row at row 2, flat 100 m elsewhere.
        dem = fill(100.0, 10, 5)
        dem[1, :] .= 200.0
        dem[2, :] .= NaN
        h = horizon_angle(dem; directions = 4, cellsize = (1.0, 1.0))
        # NaN row treated as 0 m (below observer) so peak 9 cells north dominates.
        @test h[10, 3, 1] ≈ atand(100 / 9) atol = 1e-3
        @test isnan(h[2, 3, 1])
    end

    @testset "missing_elevation = Inf fully occludes" begin
        dem = fill(100.0, 10, 5)
        dem[1, :] .= 200.0
        dem[2, :] .= NaN
        h = horizon_angle(dem; directions = 4, cellsize = (1.0, 1.0), missing_elevation = Inf)
        @test h[10, 3, 1] ≈ 90.0 atol = 1e-3
    end

    @testset "missing_elevation custom value" begin
        # NaN row 8 cells north of observer; substitute 150 m → 50 m rise.
        dem = fill(100.0, 10, 5)
        dem[2, :] .= NaN
        h = horizon_angle(dem; directions = 4, cellsize = (1.0, 1.0), missing_elevation = 150.0)
        @test h[10, 3, 1] ≈ atand(50 / 8) atol = 1e-3
    end

    @testset "horizon scales with cellsize" begin
        # 5 m bump at the south edge, observer 4 cells north.
        dem = fill(100.0, 5, 3)
        dem[5, :] .= 105.0
        h1 = horizon_angle(dem; directions = 4, cellsize = (1.0, 1.0))
        h2 = horizon_angle(dem; directions = 4, cellsize = (2.0, 2.0))
        @test h1[1, 2, 3] ≈ atand(5 / 4) atol = 1e-3
        @test h2[1, 2, 3] ≈ atand(5 / 8) atol = 1e-3
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

end
