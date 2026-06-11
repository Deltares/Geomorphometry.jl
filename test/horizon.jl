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
    @testset "viewshed" begin
        # Flat terrain: every cell is visible from any observer (dense Bool).
        flat = ones(11, 11) * 100.0
        v = viewshed(flat, (6, 6); cellsize = (1.0, 1.0))
        @test size(v) == (11, 11)
        @test eltype(v) == Bool
        @test all(v)                  # everything visible on a flat plane
        @test v[6, 6]                 # observer

        # A wall hides everything directly behind it.
        dem = fill(100.0, 10, 3)
        dem[5, :] .= 110.0
        v = viewshed(dem, (10, 2); cellsize = (1.0, 1.0))
        @test all(v[r, 2] for r in 5:10)        # observer up to the wall
        @test !any(v[r, 2] for r in 1:4)        # hidden behind the wall

        # Off-axis line of sight: a spike between observer and target hides it.
        dem2 = fill(100.0, 7, 7)
        dem2[4, 5] = 130.0                       # spike on the way to (1, 7)
        v2 = viewshed(dem2, (7, 3); cellsize = (1.0, 1.0))
        @test !v2[1, 7]                          # blocked by the off-axis spike
        @test v2[7, 3]                           # observer

        # A peak observer sees more than one in a pit.
        cone = fill(100.0, 21, 21)
        for i in 1:21, j in 1:21
            cone[i, j] = 100.0 + max(0.0, 10.0 - hypot(i - 11, j - 11))
        end
        @test count(viewshed(cone, (11, 11); cellsize = (1.0, 1.0))) >
              count(viewshed(cone, (1, 1); cellsize = (1.0, 1.0)))

        # Accepts a CartesianIndex observer.
        vc = viewshed(flat, CartesianIndex(6, 6); cellsize = (1.0, 1.0))
        @test all(vc)

        # observer_height lets the observer see over a low obstacle.
        demh = fill(100.0, 1, 9)
        demh[1, 5] = 101.0                       # 1 m bump between observer and target
        @test !viewshed(demh, (1, 1); cellsize = (1.0, 1.0))[1, 9]
        @test viewshed(demh, (1, 1); cellsize = (1.0, 1.0), observer_height = 3.0)[1, 9]

        # NaN handling: a NaN target is not visible; a NaN observer hides all.
        demn = fill(100.0, 5, 5)
        demn[1, 1] = NaN
        @test !viewshed(demn, (5, 5); cellsize = (1.0, 1.0))[1, 1]
        @test !any(viewshed(demn, (1, 1); cellsize = (1.0, 1.0)))

        # Invalid observer.
        @test_throws ArgumentError viewshed(flat, (0, 6))
        @test_throws ArgumentError viewshed(flat, (6, 99))
    end

    @testset "total_viewshed" begin
        # Flat terrain: everything visible everywhere.
        flat = ones(12, 12) * 100.0
        tv = total_viewshed(flat; directions = 8, cellsize = (1.0, 1.0))
        @test size(tv) == (12, 12)
        @test all(tv .≈ 1.0f0)

        # Uniform slope: no occlusion anywhere. A small observer height makes the
        # grazing line of sight robust (every cell visible).
        slope = [Float64(r) for r in 1:12, c in 1:12]
        tvs = total_viewshed(slope; directions = 8, cellsize = (1.0, 1.0), observer_height = 1.6)
        @test all(tvs .≈ 1.0f0)

        # A peak sees everything; cells in its shadow see less. Observer height of
        # ~1.6 m lets the peak see over the grazing cone surface.
        dem = fill(100.0, 21, 21)
        for i in 1:21, j in 1:21
            d = hypot(i - 11, j - 11)
            dem[i, j] = 100.0 + max(0.0, 10.0 - d)
        end
        tvp = total_viewshed(dem; directions = 8, cellsize = (1.0, 1.0), observer_height = 1.6)
        @test tvp[11, 11] ≈ 1.0f0          # peak sees all
        @test tvp[1, 1] < tvp[11, 11]      # corner partly blocked by the peak
        @test all(0.0 .<= tvp .<= 1.0)

        # Without an observer height the grazing cone surface is not fully visible.
        tvp0 = total_viewshed(dem; directions = 8, cellsize = (1.0, 1.0))
        @test tvp0[11, 11] < 1.0f0

        # Rotation path (>8 directions) runs and stays in range.
        tv16 = total_viewshed(dem; directions = 16, cellsize = (1.0, 1.0), observer_height = 1.6)
        @test all(0.0 .<= tv16 .<= 1.0)

        @test_throws ArgumentError total_viewshed(flat; directions = 6)
    end

end
