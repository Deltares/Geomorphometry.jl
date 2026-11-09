# Usage

```@meta
CurrentModule = Geomorphometry
```

```@setup plots
using Geomorphometry, CairoMakie, GeoArrays, Statistics # hide
import Eikonal # hide
set_theme!(theme_minimal(); transparency = true) # hide
CairoMakie.activate!(type = "png")
A = GeoArrays.read("saba.tif") # hide
B = GeoArrays.read("saba_dsm.tif") # hide
mask = ismissing.(A) # hide
# dem = coalesce(A, NaN32).A # hide
dtm = coalesce(A, NaN32) # hide
dsm = coalesce(B, NaN32) # hide
# ndem = GeoArrays.flipud!(deepcopy(dtm))

fdtm = deepcopy(dtm)
out = zeros(size(fdtm))
queued = falses(size(fdtm))
fig, ax, hm = heatmap(out, colormap=:turbo, colorrange=(0, 600))

openn = Geomorphometry.PriorityQueue{CartesianIndex{2}, eltype(fdtm)}()
pit = Geomorphometry.DataStructures.Queue{CartesianIndex{2}}()

R = CartesianIndices(fdtm)
I_first, I_last = first(R), last(R)
fps=60
for cell in Geomorphometry.outlets(fdtm)
    Geomorphometry.enqueue!(openn, cell, fdtm[cell])
    queued[cell] = true  # queued
end
record(fig, expanduser("test.mp4")) do io
    i = 0
    while (!isempty(openn) || !isempty(pit))
        cell = !isempty(pit) ? Geomorphometry.DataStructures.dequeue!(pit) : Geomorphometry.dequeue!(openn)
        out[cell] = fdtm[cell]
        # sleep(1/fps)
        # hm[1] = out # update data
        if iszero(i % 5_000)
            hm[1] = out # update data
            recordframe!(io)
        end
        for ncell in max(I_first, cell - Geomorphometry.Δ):min(I_last, cell + Geomorphometry.Δ)
            (queued[ncell] || ncell == cell) && continue
            queued[ncell] = true
            if fdtm[ncell] <= fdtm[cell]
                fdtm[ncell] = fdtm[cell]
                Geomorphometry.DataStructures.enqueue!(pit, ncell)
            else
                Geomorphometry.enqueue!(openn, ncell, fdtm[ncell])
            end
        end
        i += 1
    end
end
```

In Geomorphometry.jl we provide a set of tools to analyze and visualize the shape of the Earth. The package is designed to be fast, flexible, and easy to use.
Moreover, we have implemented several algorithms for a common operation so that you can choose the one that best fits your needs.

In these pages we will use the elevation model of [Saba](https://en.wikipedia.org/wiki/Saba_(island)) to showcase the different categories of operations that are available in Geomorphometry.jl. However, any DEM should work, most algorithms accept an `AbstractMatrix{<:Real}`. Note that you thus need to remove missing values from your DEM, which is ideally done by interpolation, but a quick `coalesce.(dem, NaN32)` will work too (but will cause NaN gaps in your output). Below are the steps to reproduce the example plots.

::: tabs

== Downloading
```julia
using Downloads

dtm_fn = abspath("saba.tif")
isfile(dtm_fn) || Downloads.download(
    "https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/saba.tif",
    dtm_fn,
)
dsm_fn = abspath("saba_dsm.tif")
isfile(dsm_fn) || Downloads.download(
    "https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/saba_dsm.tif",
    dsm_fn,
)
```

== Reading (using Rasters.jl)
```julia
using Rasters
import ArchGDAL

dtm = Raster(dtm_fn)
dsm = Raster(dsm_fn)
mask = ismissing.(dtm)  # used for hiding nodata areas in plots
dtm = replace_missing(dtm, NaN32)
dsm = replace_missing(dsm, NaN32)
```

== Reading (using GeoArrays.jl)
```julia
using GeoArrays

dtm = GeoArrays.read(dtm_fn)
dsm = GeoArrays.read(dsm_fn)
mask = ismissing.(dtm)  # used for hiding nodata areas in plots
dtm = coalesce(dtm, NaN32)
dsm = coalesce(dsm, NaN32)
```

== Plotting
```julia
using CairoMakie  # or GLMakie, Plots

heatmap(dtm)
```
:::

::: details

You can download the data yourself from https://www.beeldmateriaal.nl/dataroom-caribisch-nederland, specifically the [DTM](https://basisdata.nl/hwh-ortho/Caribisch_NL/LiDAR/caribisch/Saba/02a_DTM_50cm/M_SABA2023.TIF) and [DSM](https://basisdata.nl/hwh-ortho/Caribisch_NL/LiDAR/caribisch/Saba/03a_DSM_50cm/R_SABA2023.TIF). I've resampled them to a coarser resolution (from 0.5 m to 5 m), interpolated some nodata areas, and attached them to a Github release https://github.com/Deltares/Geomorphometry.jl/releases/tag/v0.6.0.

:::

## Visualization
Visualization is done using the [`hillshade`](@ref), [`multihillshade`](@ref), and [`pssm`](@ref) functions. The first two shade the terrain by using a single or multiple light source(s) respectively, while `pssm` is a slope map exaggerated for human perception.


:::tabs

== Hillshade
```@example plots
heatmap(hillshade(dtm))
```

== Multihillshade
```@example plots
heatmap(multihillshade(dtm, azimuth=0:60:270, zenith=60))
```

== PSSM
```@example plots
heatmap(pssm(dtm), colormap=Reverse(:viridis))
```
:::

We can also overlay the `pssm` visualization on top of a colored `aspect` map.

```@example plots
f = heatmap(aspect(dtm); colormap=:curl)
heatmap!(pssm(dtm); colormap=Reverse(:greys), alpha=0.5)
f
```

## Derivatives
Common derivatives are implemented in Geomorphometry.jl. These include [`slope`](@ref), [`aspect`](@ref), and `curvature`. The latter is ill-defined, here we provide [`plan_curvature`](@ref) (also called *projected contour curvature*), [`profile_curvature`](@ref) (also called *normal slope line curvature*), and [`tangential_curvature`](@ref) (also called *normal contour curvature*). Note that functions here allow for a custom radius, as demonstrated for `profile_curvature`.


### Slope

:::tabs

== Horn
```@example plots
heatmap(slope(dtm; method=Horn()); colormap=:matter, colorrange=(0, 60))
```
== ZevenbergenThorne
```@example plots
heatmap(slope(dtm; method=ZevenbergenThorne()); colormap=:matter, colorrange=(0, 60))
```
== MaximumDownwardGradient
```@example plots
heatmap(slope(dtm; method=Geomorphometry.MaximumDownwardGradient()); colormap=:matter, colorrange=(0, 60))
```

:::

We can also determine the slope of the terrain along a specific direction.

:::tabs

== Horn (0°)
```@example plots
heatmap(slope(dtm; method=Horn(), direction=0); colormap=:matter, colorrange=(-45, 45))
```
== ZevenbergenThorne (0°)
```@example plots
heatmap(slope(dtm; method=ZevenbergenThorne(), direction=0); colormap=:matter, colorrange=(-45, 45))
```
== Horn (90°)
```@example plots
heatmap(slope(dtm; method=Horn(), direction=90); colormap=:matter, colorrange=(-45, 45))
```
== ZevenbergenThorne (90°)
```@example plots
heatmap(slope(dtm; method=ZevenbergenThorne(), direction=90); colormap=:matter, colorrange=(-45, 45))
```

:::

### Aspect

:::tabs

== Horn
```@example plots
heatmap(aspect(dtm; method=Horn()); colormap=:romaO)
```
== ZevenbergenThorne
```@example plots
heatmap(aspect(dtm; method=ZevenbergenThorne()); colormap=:romaO)
```
== MaximumDownwardGradient
```@example plots
heatmap(aspect(dtm; method=Geomorphometry.MaximumDownwardGradient()); colormap=:romaO)
```

:::

### Curvature

:::tabs

== Profile curvature
```@example plots
heatmap(profile_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

== Plan curvature
```@example plots
heatmap(plan_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

== Tangential curvature
```@example plots
heatmap(tangential_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

== Laplacian
```@example plots
heatmap(laplacian(dtm); colorrange=(-1,1), colormap=:tarn)
```

:::

We can also determine the curvature of the terrain along a specific direction.

:::tabs

== Laplacian in 90° direction
```@example plots
heatmap(laplacian(dtm, radius=1, direction=90); colorrange=(-1,1), colormap=:tarn)
```

== Laplacian in 0° direction
```@example plots
heatmap(laplacian(dtm, radius=1, direction=0); colorrange=(-1,1), colormap=:tarn)
```

:::

## Relative position
There are several terrain descriptors that can be used to analyze the relative position of a point with respect to its neighbors. These include [`topographic_position_index`](@ref), [`terrain_ruggedness_index`](@ref), [`roughness_index_elevation`](@ref), [`bathymetric_position_index`](@ref), [`rugosity`](@ref), [`roughness`](@ref) and [`percentile_elevation`](@ref). Here we use `bathymetric_position_index`, but with a custom sized `Window` (from Stencils.jl). All these functions can be used with a custom window size.

:::tabs

== Bathymetric Position Index (BPI)
```@example plots
heatmap(bathymetric_position_index(dtm, Geomorphometry.Annulus(5, 3)); colormap=:delta, colorrange=(-10,10))
```
== Topographic Position Index (TPI)
```@example plots
heatmap(topographic_position_index(dtm); colormap=:delta, colorrange=(-10,10))
```
== Terrain Ruggedness Index (TRI)
```@example plots
heatmap(terrain_ruggedness_index(dtm); colormap=:speed)
```
== Roughness Index Elevation (RIE)
```@example plots
heatmap(roughness_index_elevation(dtm); colormap=:speed)
```
== rugosity
```@example plots
heatmap(rugosity(dtm); colormap=:speed)
```
== roughness
```@example plots
heatmap(roughness(dtm); colormap=:speed)
```
== Percentile Elevation
```@example plots
heatmap(percentile_elevation(dtm; radius=5); colormap=:delta, colorrange=(0,1))
```

:::

## Horizon
Related to the relative position is the category related to the view (to a horizon) from any cell in a DEM. We have [`horizon_angle`](@ref) to compute the maximum horizon angle for a given number of directions, and [`sky_view_factor`](@ref) to compute fraction of the sky (hemisphere) visible from a cell. The same running-maximum test powers [`viewshed`](@ref), the line-of-sight visibility map of a single observer cell, and [`total_viewshed`](@ref), a normalized visibility index for every cell.

:::tabs

== Horizon Angle (to the west)
```@example plots
hor = horizon_angle(dtm)  # this produces a width,height,directions sized raster
heatmap(hor[:,:,1]; colormap=:curl, colorrange=(-90, 90))  # first direction (east) out of 16
```
== Sky View Factor (SVF) in 16 directions
```@example plots
skf = sky_view_factor(dtm)
heatmap(skf; colormap=:ice)
```
== Sky View Factor (SVF) in 64 directions
```@example plots
skf = sky_view_factor(dtm; directions=64)
heatmap(skf; colormap=:ice)
```
== Viewshed
```@example plots
point = CartesianIndex(400,400)
vs = viewshed(dtm, point; observer_height=500)
f = heatmap(vs; colormap=:ice, colorrange=(-0.5,1.5))
scatter!(GeoArrays.coords(dtm, point); color=:orange)
f
```
== Total Viewshed
```@example plots
tvs = total_viewshed(dtm)
heatmap(tvs; colormap=:ice)
```

:::


## Landforms
The same lines of sight also classify the terrain itself. [`geomorphon`](@ref) reduces the view in eight directions to a ternary pattern of rising, flat and falling directions, and maps it onto one of ten [`Landform`](@ref) classes. The `radius` keyword sets the scale of the classification, while `flatness` sets how much relief a direction needs before it is no longer considered flat.

:::tabs

== Geomorphons
```@example plots
forms = geomorphon(dtm)
heatmap(forms; colormap=:tab10, colorrange=(0.5, 10.5), lowclip=:transparent)
```
== Geomorphons at a larger scale
```@example plots
forms = geomorphon(dtm; radius=25, skip=2)
heatmap(forms; colormap=:tab10, colorrange=(0.5, 10.5), lowclip=:transparent)
```
== Ridges only
```@example plots
forms = geomorphon(dtm)
heatmap(forms .== Ridge; colormap=:ice, colorrange=(-0.5, 1.5))
```

:::


## Hydrology
Hydrological operations are used to analyze the flow of water on the terrain. We provide [`filldepressions`](@ref) to fill depressions, [`depression_depth`](@ref) to calculate the depth of each depression (difference between filled dem and dem) and [`depression_volume`](@ref) that sums all depression depths. The major depression in the example is the caldera of the dormant volcano ([Mount Scenery](https://en.wikipedia.org/wiki/Mount_Scenery)).


:::tabs

== Fill depressions
```@example plots
cdtm = copy(dtm)
cdtm[mask] .= -Inf  # Sea should be lower than terrain
fdtm = filldepressions(cdtm)
heatmap(fdtm)
```
== Depression depth
```@example plots
depth = depression_depth(dtm; filled=fdtm)
depth[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(depth; colormap=:tempo)
```
== Depression volume
```@example plots
depression_volume(dtm; filled=fdtm)
```

:::

### Lowest neighbour is not steepest descent

[`flowaccumulation`](@ref) carves: one Priority-Flood pass gives each cell a direction as
it goes, routing flow through depressions without changing any elevation, after
[barnesPriorityFloodOptimalDepressionFilling2014](@citet) and the least-cost search of
[wangEfficientMethodIdentifying2006](@citet) and [metzEfficientExtractionDrainage2011](@citet).

!!! info "Why the direction is chosen when a cell is popped"
    A flood reaches a cell from whichever neighbour pops first, and it pops in elevation
    order, so the rule it implies is *lowest neighbour*. [`D8`](@ref) is *steepest
    descent*. On a square grid these differ, because a diagonal is `√2` away: it must drop
    41.4% further to win on gradient, but wins on elevation as soon as it drops further at
    all. With a centre at 10.0 and unit cells, an eastern neighbour at 9.0 has gradient
    1.000 while a north-eastern one at 8.8 has only 0.849 — lower, but shallower.

    So a cell settles its direction when it is *popped*, from the neighbours already
    popped, which are exactly the ones below it. On a drained surface the carved and
    derived directions then agree for 99.97% of interior cells of this tile; what remains
    is cells beside nodata, or with nowhere downhill to go.

### Filling and flow directions

[`filldepressions`](@ref) takes a `method`: the default [`PriorityFlood`](@ref) of
[barnesPriorityFloodOptimalDepressionFilling2014](@citet), or [`Zhou`](@ref) after
[zhouEfficientVariantPriorityFlood2016](@citet), which reaches the same surface 1.1-1.4x
faster. Filling leaves flats, which have no flow direction, so by default it also gives
them a gradient with [`resolveflats`](@ref) after
[barnesEfficientAssignmentDrainage2014](@citet); `flats=false` gives the filled surface
alone. [`flowdirection`](@ref) then reads directions off that surface, with the share each
receiver takes beside them.

```@example plots
z = filldepressions(cdtm)
directions, weights = flowdirection(z; method=D8())
parent(directions[400:800, 300:700])  # arrows, as the drainage direction below
```

[`D8`](@ref) has one receiver that takes everything, so its weights are a `Fill` costing
nothing at any grid size. [`DInf`](@ref) and [`FD8`](@ref) store a share per receiver, two
and eight, which `weighttype=Float16` halves.

[`flowaccumulation`](@ref) then accumulates over them, after
[barnesParallelNonDivergentFlow2017](@citet): each cell is released to its receivers once
all of its own contributors have been, so the grid is traversed in dependency order without
sorting it. There is no `method` to pass, as the directions already settle it.

```@example plots
acc = flowaccumulation(z, directions, weights)
nothing # hide
```

Filling and then deriving is the faster of the two routes on this tile, by 20 to 40%, and
the fill is most of it — deriving directions and accumulating over them are cheap by
comparison, so the choice between carving and filling is a hydrological one rather than a
performance one.

| | D8 | D∞ | FD8 |
| --- | --- | --- | --- |
| **carve**, one Priority-Flood | 2.4 s | 3.4 s | 3.0 s |
| **fill**, three passes | 2.0 s | 2.4 s | 2.5 s |
| — `filldepressions` (fill + flats) | 1.6 s | 1.6 s | 1.6 s |
| — `flowdirection` | 0.15 s | 0.37 s | 0.38 s |
| — `flowaccumulation` over those directions | 0.06 s | 0.11 s | 0.20 s |

Filling the depressions in a DEM is not necessary to calculate the flow accumulation. Here we use [`flowaccumulation`](@ref) to do so, and it automatically carves out depressions. Note that the local drainage direction is also returned. By default the FD8 algorithm is used, but you can also use the D∞ or D8 algorithm by setting the `method` keyword argument to `DInf()` or `D8()`.

::: warning

Most hydrology related methods iterate through all DEM cells, thus producing output for all DEM cells, including cells containing nodata values (e.g. `NaN32`, `Inf`).
Here we manually hide these nodata areas for plotting.

:::


:::tabs

== Flow accumulation with FD8
```@example plots
acc, ldd = flowaccumulation(dtm; method=FD8(2))
acc[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(log10.(acc); colormap=:rain)
```
== Flow accumulation with D∞
```@example plots
acc, ldd = flowaccumulation(dtm; method=DInf())
acc[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(log10.(acc); colormap=:rain)
```
== Flow accumulation with D8
```@example plots
acc, ldd = flowaccumulation(dtm; method=D8())
acc[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(log10.(acc); colormap=:rain)
```
== Drainage direction

The ldd is a Matrix of FlowDirections that visualizes the integer encodings.
```@example plots
acc, ldd = flowaccumulation(dtm; method=D8())
parent(ldd[400:800, 300:700])  # overland flow directions
```
== Underlying method

We use the Priority Flood method by [barnesPriorityFloodOptimalDepressionFilling2014](@citet).

```@raw html
<video autoplay loop muted playsinline controls src="./test.mp4" />
```

:::

Combined with the previous analysis, we can calculate [`topographic_wetness_index`](@ref), [`stream_power_index`](@ref) and 
[`drainage_potential`](@ref) indices. These indices are used to analyze the terrain's ability to accumulate water, the terrain's ability to store water respectively, and the ability to drain water respectively.

:::tabs

== Topographic Wetness Index (TWI)
```@example plots
twi = topographic_wetness_index(dtm; method=FD8())
heatmap(twi; colormap=:tempo)
```
== Stream Power Index (SPI)
```@example plots
twi = stream_power_index(dtm; method=FD8())
heatmap(twi; colormap=:tempo)
```
== Drainage Potential
```@example plots
dp = drainage_potential(dtm; method=FD8())
heatmap(dp; colormap=:tempo)
```

:::

We can also calculate the Height Above Nearest Drainage (HAND) using the [`height_above_nearest_drainage`](@ref) function. This function requires a flow accumulation map to determine the stream network. Here we use a threshold based on the flow accumulation to define streams.

:::tabs

== HAND threshold 1e3
```@example plots
hand = height_above_nearest_drainage(dtm; threshold=1e3)
hand[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(hand; colormap=:turbo)
```
== HAND threshold 1e5
```@example plots
hand = height_above_nearest_drainage(dtm; threshold=1e5)
hand[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(hand; colormap=:turbo)
```
== HAND threshold 1e5 with FD8 method
```@example plots
hand = height_above_nearest_drainage(dtm; method=FD8(), threshold=1e5)
hand[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(hand; colormap=:turbo)
```

:::


### Watersheds

[`watersheds`](@ref) labels the sets of cells draining to the same outlet, and there are
two ways to it, matching the two ways of making a DEM drain.

Given a direction field it labels over that, so it answers for whatever flow metric
produced it. It can only send a cell where the directions point, so it needs a DEM that
already drains. Given only a DEM it floods inwards from the rim instead, exactly as
[`filldepressions`](@ref) does, so a hollow with no outlet joins the basin that fills it
rather than becoming one of its own. Nothing needs preparing and the DEM is not modified.

The difference shows on a raw DEM: flooding finds 5400 watersheds on this tile, one per
cell on the rim where water leaves, while following directions on the same raw DEM finds
300822, because every depression is somewhere flow arrives and stops.

A cell that sends flow several ways joins the one taking the largest share, since
watersheds partition the grid. So the labelling depends on the flow metric even though the
outlets do not: on the filled surface all three methods find the same 2135 watersheds, but
[`DInf`](@ref) puts 95.5% of cells in the same one as [`D8`](@ref), and [`FD8`](@ref) only
72.1%.

:::tabs

== Flooded, raw DEM
```@example plots
heatmap(Float32.(watersheds(dtm)); colormap=:glasbey_bw_minc_20_n256)
```
== D8, filled
```@example plots
heatmap(Float32.(watersheds(z, directions, weights)); colormap=:glasbey_bw_minc_20_n256)
```
== FD8, filled
```@example plots
fd8directions, fd8weights = flowdirection(z; method=FD8())
heatmap(Float32.(watersheds(z, fd8directions, fd8weights)); colormap=:glasbey_bw_minc_20_n256)
```

:::

The flood is always [`PriorityFlood`](@ref), unlike [`filldepressions`](@ref). Taking cells
in elevation order is what makes the labelling meaningful: a cell is reached by its lowest
neighbour and takes that label. [`Zhou`](@ref) gives up that ordering to go faster, and a
cell then takes its label from whichever neighbour the region growing arrived from,
draining to it or not, so it is not offered here.

## Terrain filters
While filters seem unrelated to the previously discussed methods, these are used to filter (rasterized) pointclouds i.e. surface models (DSM) to arrive at a terrain model (DTM). We provide [`progressive_morphological_filter`](@ref), [`simple_morphological_filter`](@ref), and [`skewness_balancing`](@ref) filters. Here we use the `progressive_morphological_filter` filter to remove elevations that do not pass the slope filter (normally vegetation and buildings).


:::tabs

== Progressive Morphological Filter (PMF)
```@example plots
B, flags = progressive_morphological_filter(dsm, ωₘ = 15.0, slope = 0.2, dhₘ = 5.0, dh₀ = 0.3)
A = copy(dsm)
A[A .> B] .= NaN32
heatmap(A)
```
== Simple Morphological Filter (SMF)
```@example plots
A = simple_morphological_filter(dsm; ω=15.0, slope=0.2)
heatmap(A)
```
== Skewness Balancing (SKB)
```@example plots
smask = skewness_balancing(dsm)
A = copy(dsm)
A[.!smask] .= NaN32
heatmap(A)
```
:::


## Cost distance
Finally there are so-called cost distance functions. Instead of calculating a distance in the Euclidean sense, these functions calculate the cost of moving from one point to another, given a spatial varying cost. We provide a single [`spread`](@ref) function with multiple `method`s. The function takes starting points, an initial cost, and a cost (friction) map.

Here we use the top of the volcano as the starting point, and the squared slope as the cost map. We compare the `Tomlin` method with the `Eikonal` and `Eastman` methods. `Eastman` requires a high number of iterations to converge (even at 25 it's not enough), and is thus slower.

:::tabs

== Tomlin
```@example plots
tspread = spread(dtm .> 850, 0, slope(dtm).^2)
tspread[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(tspread, colormap=:dense)
```

== Eikonal
```@example plots
fspread = spread(dtm .> 850, 0, slope(dtm).^2, method=FastSweeping())
fspread[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(fspread, colormap=:dense)
```

== Eastman
```@example plots
espread = spread(dtm .> 850, 0, slope(dtm).^2, method=Eastman(iterations=25))
espread[mask] .= NaN32  # ignore nodata areas for plotting
heatmap(espread, colormap=:dense)
```

:::
