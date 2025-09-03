
# Usage {#Usage}



In Geomorphometry.jl we provide a set of tools to analyze and visualize the shape of the Earth. The package is designed to be fast, flexible, and easy to use. Moreover, we have implemented several algorithms for a common operation so that you can choose the one that best fits your needs.

In these pages we will use the elevation model of [Saba](https://en.wikipedia.org/wiki/Saba_(island)) to showcase the different categories of operations that are available in Geomorphometry.jl.

## Visualization {#Visualization}

Visualization is done using the [`hillshade`](/reference#Geomorphometry.hillshade-Tuple{AbstractMatrix{<:Real}}), [`multihillshade`](/reference#Geomorphometry.multihillshade-Tuple{AbstractMatrix{<:Real}}), and [`pssm`](/reference#Geomorphometry.pssm-Tuple{AbstractMatrix{<:Real}}) functions. The first two shade the terrain by using a single or multiple light source(s) respectively, while `pssm` is a slope map exaggerated for human perception.

:::tabs

== Hillshade

```julia
heatmap(hillshade(dtm))
```

![](ndgzksw.png){width=600px height=450px}

== Multihillshade

```julia
heatmap(multihillshade(dtm, azimuth=0:60:270, zenith=60))
```

![](xpatrdk.png){width=600px height=450px}

== PSSM

```julia
heatmap(pssm(dtm), colormap=Reverse(:viridis))
```

![](fsqvtnr.png){width=600px height=450px}

:::

We can also overlay the `pssm` visualization on top of a colored `aspect` map.

```julia
f = heatmap(aspect(dtm); colormap=:curl)
heatmap!(pssm(dtm); colormap=Reverse(:greys), alpha=0.5)
f
```

![](cizbbuj.png){width=600px height=450px}

## Derivatives {#Derivatives}

Common derivatives are implemented in Geomorphometry.jl. These include [`slope`](/reference#Geomorphometry.slope-Tuple{AbstractMatrix{<:Real}}), [`aspect`](/reference#Geomorphometry.aspect-Tuple{AbstractMatrix{<:Real}}), and `curvature`. The latter is ill-defined, here we provide [`plan_curvature`](/reference#Geomorphometry.plan_curvature-Tuple{AbstractMatrix{<:Real}}) (also called _projected contour curvature_), [`profile_curvature`](/reference#Geomorphometry.profile_curvature-Tuple{AbstractMatrix{<:Real}}) (also called _normal slope line curvature_), and [`tangential_curvature`](/reference#Geomorphometry.tangential_curvature-Tuple{AbstractMatrix{<:Real}}) (also called _normal contour curvature_). Note that functions here allow for a custom radius (but fixed positions, see X), as demonstrated for `profile_curvature`.

### Slope {#Slope}

:::tabs

== Horn

```julia
heatmap(slope(dtm; method=Horn()); colormap=:matter, colorrange=(0, 60))
```

![](oophaho.png){width=600px height=450px}

== ZevenbergenThorne

```julia
heatmap(slope(dtm; method=ZevenbergenThorne()); colormap=:matter, colorrange=(0, 60))
```

![](eskqgbl.png){width=600px height=450px}

== MaximumDownwardGradient

```julia
heatmap(slope(dtm; method=Geomorphometry.MaximumDownwardGradient()); colormap=:matter, colorrange=(0, 60))
```

![](mjsipux.png){width=600px height=450px}

:::

We can also determine the slope of the terrain along a specific direction.

:::tabs

== Horn (0°)

```julia
heatmap(slope(dtm; method=Horn(), direction=0); colormap=:matter, colorrange=(-45, 45))
```

![](vmxlpye.png){width=600px height=450px}

== ZevenbergenThorne (0°)

```julia
heatmap(slope(dtm; method=ZevenbergenThorne(), direction=0); colormap=:matter, colorrange=(-45, 45))
```

![](iiuevgo.png){width=600px height=450px}

== Horn (90°)

```julia
heatmap(slope(dtm; method=Horn(), direction=90); colormap=:matter, colorrange=(-45, 45))
```

![](juvawiv.png){width=600px height=450px}

== ZevenbergenThorne (90°))

```julia
heatmap(slope(dtm; method=ZevenbergenThorne(), direction=90); colormap=:matter, colorrange=(-45, 45))
```

![](wenswih.png){width=600px height=450px}

:::

### Aspect {#Aspect}

:::tabs

== Horn

```julia
heatmap(aspect(dtm; method=Horn()); colormap=:romaO)
```

![](ugjlckj.png){width=600px height=450px}

== ZevenbergenThorne

```julia
heatmap(aspect(dtm; method=ZevenbergenThorne()); colormap=:romaO)
```

![](paxgodg.png){width=600px height=450px}

== MaximumDownwardGradient

```julia
heatmap(aspect(dtm; method=Geomorphometry.MaximumDownwardGradient()); colormap=:romaO)
```

![](fkrgjub.png){width=600px height=450px}

:::

### Curvature {#Curvature}

:::tabs

== Profile curvature

```julia
heatmap(profile_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

![](dsrdxfa.png){width=600px height=450px}

== Plan curvature

```julia
heatmap(plan_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

![](iervqor.png){width=600px height=450px}

== Tangential curvature

```julia
heatmap(tangential_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

![](adxdamr.png){width=600px height=450px}

== Laplacian

```julia
heatmap(laplacian(dtm); colorrange=(-1,1), colormap=:tarn)
```

![](hawrtft.png){width=600px height=450px}

:::

We can also determine the curvature of the terrain along a specific direction.

:::tabs

== Laplacian in 90° direction

```julia
heatmap(laplacian(dtm, radius=1, direction=90); colorrange=(-1,1), colormap=:tarn)
```

![](xmepwvp.png){width=600px height=450px}

== Laplacian in 0° direction

```julia
heatmap(laplacian(dtm, radius=1, direction=0); colorrange=(-1,1), colormap=:tarn)
```

![](ccvbvww.png){width=600px height=450px}

:::

## Relative position {#Relative-position}

There are several terrain descriptors that can be used to analyze the relative position of a point with respect to its neighbors. These include [`TPI`](/reference#Geomorphometry.TPI), [`TRI`](/reference#Geomorphometry.TRI-Tuple{AbstractMatrix{<:Real}}), [`RIE`](/reference#Geomorphometry.RIE), [`BPI`](/reference#Geomorphometry.BPI), [`rugosity`](/reference#Geomorphometry.rugosity-Tuple{AbstractMatrix{<:Real}}) and [`roughness`](/reference#Geomorphometry.roughness). Here we use `BPI`, but with a custom sized `Window` (from Stencils.jl). All these functions can be used with a custom window size.

:::tabs

== BPI

```julia
heatmap(BPI(dtm, Geomorphometry.Annulus(5, 3)); colormap=:delta, colorrange=(-10,10))
```

![](bizovmh.png){width=600px height=450px}

== TPI

```julia
heatmap(TPI(dtm); colormap=:delta, colorrange=(-10,10))
```

![](vuwswht.png){width=600px height=450px}

== TRI

```julia
heatmap(TRI(dtm); colormap=:speed)
```

![](selwrso.png){width=600px height=450px}

== RIE

```julia
heatmap(RIE(dtm); colormap=:speed)
```

![](rrarnac.png){width=600px height=450px}

== rugosity

```julia
heatmap(rugosity(dtm); colormap=:speed)
```

![](gnmskoh.png){width=600px height=450px}

== roughness

```julia
heatmap(roughness(dtm); colormap=:speed)
```

![](mxacuau.png){width=600px height=450px}

:::

## Hydrology {#Hydrology}

Hydrological operations are used to analyze the flow of water on the terrain. We provide [`filldepressions`](/reference#Geomorphometry.filldepressions) to fill depressions, and [`flowaccumulation`](/reference#Geomorphometry.flowaccumulation) to calculate the flow accumulation. Here we use `flowaccumulation` to calculate the flow accumulation. Note that the local drainage direction is also returned. By default the FD8 algorithm is used, but you can also use the D∞ or D8 algorithm by setting the `method` keyword argument to `DInf()` or `D8()`.

:::tabs

== Flow accumulation with FD8

```julia
acc, ldd = flowaccumulation(dtm; method=FD8(2))
heatmap(log10.(acc); colormap=:rain)
```

![](lvwnabj.png){width=600px height=450px}

== Flow accumulation with D∞

```julia
acc, ldd = flowaccumulation(dtm; method=DInf())
heatmap(log10.(acc); colormap=:rain)
```

![](amptzpj.png){width=600px height=450px}

== Flow accumulation with D8

```julia
acc, ldd = flowaccumulation(dtm; method=D8())
heatmap(log10.(acc); colormap=:rain)
```

![](boigrvf.png){width=600px height=450px}

== Underlying method

We use the Priority Flood method by [Barnes _et al._ (2014)](/bibliography#barnesPriorityFloodOptimalDepressionFilling2014).
<video autoplay loop muted playsinline controls src="./test.mp4" />


:::

Combined with the previous analysis, we can calculate [`TWI`](/reference#Geomorphometry.TWI-Tuple{AbstractMatrix}) and [`SPI`](/reference#Geomorphometry.SPI-Tuple{AbstractMatrix}) indices. These indices are used to analyze the terrain&#39;s ability to accumulate water and the terrain&#39;s ability to store water respectively.

:::tabs

== TWI

```julia
twi = TWI(dtm; method=FD8())
heatmap(twi; colormap=:tempo)
```

![](nebbgoi.png){width=600px height=450px}

== SPI

```julia
twi = SPI(dtm; method=FD8())
heatmap(twi; colormap=:tempo)
```

![](erctftq.png){width=600px height=450px}

:::

## Terrain filters {#Terrain-filters}

While filters seem unrelated to the previously discussed methods, these are used to filter (rasterized) pointclouds i.e. surface models (DSM) to arrive at a terrain model (DTM). We provide [`pmf`](/reference#Geomorphometry.pmf-Tuple{AbstractMatrix{<:Real}}), [`smf`](/reference#Geomorphometry.smf-Tuple{AbstractMatrix{<:Real}}), and [`skb`](/reference#Geomorphometry.skb-Tuple{AbstractArray}) filters. Here we use the `pmf` filter to remove elevations that do not pass the slope filter (normally vegetation and buildings).

:::tabs

== PMF

```julia
B, flags = pmf(dsm, ωₘ = 15.0, slope = 0.2, dhₘ = 5.0, dh₀ = 0.3)
A = copy(dsm)
A[A .> B] .= NaN
heatmap(A)
```

![](rbtbwzs.png){width=600px height=450px}

== SMF

```julia
A = smf(dsm; ω=15.0, slope=0.2)
heatmap(A)
```

![](vvbycfz.png){width=600px height=450px}

== SKB

```julia
B, flags = skb(dsm)
A = copy(dsm)
A[A .> B] .= NaN
heatmap(A)
```

![](atvzeoo.png){width=600px height=450px}

:::

## Cost distance {#Cost-distance}

Finally there are so-called cost distance functions. Instead of calculating a distance in the Euclidean sense, these functions calculate the cost of moving from one point to another, given a spatial varying cost. We provide a single [`spread`](/reference#Geomorphometry.spread-Tuple{AbstractMatrix{<:Real},%20AbstractMatrix{<:Real},%20Real}) function with multiple `method`s. The function takes starting points, an initial cost, and a cost (friction) map.

Here we use the top of the volcano as the starting point, and the squared slope as the cost map. We compare the `Tomlin` method with the `Eikonal` and `Eastman` methods. `Eastman` requires a high number of iterations to converge (even at 25 it&#39;s not enough), and is thus slower.

:::tabs

== Tomlin

```julia
heatmap(spread(dtm .> 850, 0, slope(dtm).^2), colormap=:dense)
```

![](iwckowj.png){width=600px height=450px}

== Eikonal

```julia
heatmap(spread(dtm .> 850, 0, slope(dtm).^2, method=FastSweeping()), colormap=:dense)
```

![](xaslzxb.png){width=600px height=450px}

== Eastman

```julia
heatmap(spread(dtm .> 850, 0, slope(dtm).^2, method=Eastman(iterations=25)), colormap=:dense)
```

![](loqnurr.png){width=600px height=450px}

:::
