


## Index {#Index}
- [`Geomorphometry.D8`](#Geomorphometry.D8)
- [`Geomorphometry.DInf`](#Geomorphometry.DInf)
- [`Geomorphometry.FD8`](#Geomorphometry.FD8)
- [`Geomorphometry.Horn`](#Geomorphometry.Horn)
- [`Geomorphometry.MaximumDownwardGradient`](#Geomorphometry.MaximumDownwardGradient)
- [`Geomorphometry.SpreadMethod`](#Geomorphometry.SpreadMethod)
- [`Geomorphometry.ZevenbergenThorne`](#Geomorphometry.ZevenbergenThorne)
- [`Geomorphometry.BPI`](#Geomorphometry.BPI)
- [`Geomorphometry.RIE`](#Geomorphometry.RIE)
- [`Geomorphometry.SPI`](#Geomorphometry.SPI-Tuple{AbstractMatrix})
- [`Geomorphometry.TPI`](#Geomorphometry.TPI)
- [`Geomorphometry.TRI`](#Geomorphometry.TRI-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.TWI`](#Geomorphometry.TWI-Tuple{AbstractMatrix})
- [`Geomorphometry.aspect`](#Geomorphometry.aspect-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.cellsize`](#Geomorphometry.cellsize-Tuple{Any})
- [`Geomorphometry.entropy`](#Geomorphometry.entropy-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.entropy!`](#Geomorphometry.entropy!-Union{Tuple{T},%20Tuple{AbstractMatrix{<:Real},%20AbstractMatrix{T}}}%20where%20T<:Real)
- [`Geomorphometry.filldepressions`](#Geomorphometry.filldepressions)
- [`Geomorphometry.flowaccumulation`](#Geomorphometry.flowaccumulation)
- [`Geomorphometry.hillshade`](#Geomorphometry.hillshade-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.multihillshade`](#Geomorphometry.multihillshade-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.opening!`](#Geomorphometry.opening!-Union{Tuple{T},%20Tuple{AbstractMatrix{T},%20Integer,%20AbstractMatrix{T}}}%20where%20T<:Real)
- [`Geomorphometry.opening_circ!`](#Geomorphometry.opening_circ!-Union{Tuple{T},%20Tuple{AbstractMatrix{T},%20Integer,%20AbstractMatrix{T}}}%20where%20T<:Real)
- [`Geomorphometry.pitremoval`](#Geomorphometry.pitremoval-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.plan_curvature`](#Geomorphometry.plan_curvature-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.pmf`](#Geomorphometry.pmf-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.profile_curvature`](#Geomorphometry.profile_curvature-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.prominence`](#Geomorphometry.prominence)
- [`Geomorphometry.prominence`](#Geomorphometry.prominence-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.pssm`](#Geomorphometry.pssm-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.roughness`](#Geomorphometry.roughness)
- [`Geomorphometry.round_odd`](#Geomorphometry.round_odd-Tuple{Any})
- [`Geomorphometry.rugosity`](#Geomorphometry.rugosity-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.skb`](#Geomorphometry.skb-Tuple{AbstractArray})
- [`Geomorphometry.skbr`](#Geomorphometry.skbr-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.skbr`](#Geomorphometry.skbr)
- [`Geomorphometry.slope`](#Geomorphometry.slope-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.smf`](#Geomorphometry.smf-Tuple{AbstractMatrix{<:Real}})
- [`Geomorphometry.spread`](#Geomorphometry.spread-Tuple{AbstractMatrix{<:Real},%20AbstractMatrix{<:Real},%20Real})
- [`Geomorphometry.spread`](#Geomorphometry.spread-Tuple{Vector{CartesianIndex{2}},%20AbstractMatrix{<:Real},%20AbstractMatrix{<:Real}})
- [`Geomorphometry.spread`](#Geomorphometry.spread-Union{Tuple{T},%20Tuple{Tomlin,%20Vector{CartesianIndex{2}},%20AbstractMatrix{T},%20AbstractMatrix{<:Real}}}%20where%20T<:Real)
- [`Geomorphometry.spread`](#Geomorphometry.spread-Tuple{Eastman,%20Vector{CartesianIndex{2}},%20AbstractMatrix{<:Real},%20AbstractMatrix{<:Real}})
- [`Geomorphometry.tangential_curvature`](#Geomorphometry.tangential_curvature-Tuple{AbstractMatrix{<:Real}})


## Reference - Exported functions {#Reference-Exported-functions}
<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.D8' href='#Geomorphometry.D8'><span class="jlbinding">Geomorphometry.D8</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



D8 Flow Direction method by [Jenson (1988)](/bibliography#jensonExtractingTopographicStructure1988).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/hydrology.jl#LL33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.DInf' href='#Geomorphometry.DInf'><span class="jlbinding">Geomorphometry.DInf</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



DInf Flow Direction method by [Tarboton (1997)](/bibliography#tarbotonNewMethodDetermination1997).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/hydrology.jl#LL36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.FD8' href='#Geomorphometry.FD8'><span class="jlbinding">Geomorphometry.FD8</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



FD8 Flow Direction method by [Quin (1991)](/bibliography#quinnPredictionHillslopeFlow1991).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/hydrology.jl#LL39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.Horn' href='#Geomorphometry.Horn'><span class="jlbinding">Geomorphometry.Horn</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Third order finite difference estimator using all 8 neighbors by [Horn, (1981)](/bibliography#hornHillShadingReflectance1981).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.ZevenbergenThorne' href='#Geomorphometry.ZevenbergenThorne'><span class="jlbinding">Geomorphometry.ZevenbergenThorne</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Second order finite difference estimator using all 4 neighbors by [Zevenbergen and Thorne, (1987)](/bibliography#zevenbergen1987quantitative).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.BPI' href='#Geomorphometry.BPI'><span class="jlbinding">Geomorphometry.BPI</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
BPI(dem::AbstractMatrix{<:Real})
```


BPI stands for Bathymetric Position Index (Lundblad et al., 2006), which is defined as the difference between a central pixel and the mean of the cells in an annulus around it.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL26-L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.RIE' href='#Geomorphometry.RIE'><span class="jlbinding">Geomorphometry.RIE</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
RIE(dem::AbstractMatrix{<:Real})
```


RIE stands for Roughness Index Elevation, which quantifies the standard deviation of residual topography (Cavalli et al., 2008)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL33-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.SPI-Tuple{AbstractMatrix}' href='#Geomorphometry.SPI-Tuple{AbstractMatrix}'><span class="jlbinding">Geomorphometry.SPI</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
SPI(dem::AbstractMatrix; method=D8(), cellsize=cellsize(dem))
```


Computes the Stream Power Index (SPI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/hydrology.jl#LL371-L375" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.TPI' href='#Geomorphometry.TPI'><span class="jlbinding">Geomorphometry.TPI</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
TPI(dem::AbstractMatrix{<:Real})
```


TPI stands for Topographic Position Index, which is defined as the difference between a central pixel and the mean of its surrounding cells (see Wilson et al 2007, Marine Geodesy 30:3-35).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL17-L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.TRI-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.TRI-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.TRI</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
TRI(dem::AbstractMatrix{<:Real})
```


TRI stands for Terrain Ruggedness Index, which measures the difference between a central pixel and its surrounding cells. This algorithm uses the square root of the sum of the square of the absolute difference between a central pixel and its surrounding cells. This is recommended for terrestrial use cases.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL43-L49" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.TWI-Tuple{AbstractMatrix}' href='#Geomorphometry.TWI-Tuple{AbstractMatrix}'><span class="jlbinding">Geomorphometry.TWI</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
TWI(dem::AbstractMatrix; method=D8(), cellsize=cellsize(dem))
```


Computes the Topographic Wetness Index (TWI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/hydrology.jl#LL360-L364" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.aspect-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.aspect-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.aspect</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
aspect(dem::Matrix{<:Real}, method=Horn())
```


Aspect is direction of [`slope`](/reference#Geomorphometry.slope-Tuple{AbstractMatrix{<:Real}}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL103-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.entropy-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.entropy-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.entropy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
entropy(dem::AbstractMatrix{<:Real}; step=0.5)
```


Entropy calculates the Shannon entropy of the surrounding cells of a central cell. `step` is the bin size for the histogram.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL93-L98" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.filldepressions' href='#Geomorphometry.filldepressions'><span class="jlbinding">Geomorphometry.filldepressions</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
filldepressions(dem::AbstractMatrix, mask::Matrix{Bool})
```


Performs the Priority-Flood algorithm ([Barnes _et al._, 2014](/bibliography#barnesPriorityFloodOptimalDepressionFilling2014)) on the given digital elevation model (DEM) `dem` with an optional `mask`.

**Arguments**
- `dem::AbstractMatrix`: A 2D array representing the digital elevation model (DEM).
  
- `mask::AbstractMatrix{Bool}`: A 2D boolean array representing the mask. Cells with `true` values are considered in the computation, while cells with `false` values are ignored.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/hydrology.jl#LL16-L24" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.flowaccumulation' href='#Geomorphometry.flowaccumulation'><span class="jlbinding">Geomorphometry.flowaccumulation</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
flowaccumulation(dem::AbstractMatrix, closed::Matrix{Bool}, method::FlowDirectionMethod)
```


Computes the flow accumulation of a digital elevation model (DEM) `dem` with an optional `closed` mask and a `method` for flow direction. Returns the flow accumulation and the flow direction (local drainage direction or ldd)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/hydrology.jl#LL185-L190" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.hillshade-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.hillshade-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.hillshade</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
hillshade(dem::Matrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=cellsize(dem))
```


hillshade is the simulated illumination of a surface based on its [`slope`](/reference#Geomorphometry.slope-Tuple{AbstractMatrix{<:Real}}) and [`aspect`](/reference#Geomorphometry.aspect-Tuple{AbstractMatrix{<:Real}}) given a light source with azimuth and zenith angles in °, as defined in [Burrough, P. A., and McDonell, R. A., (1998)](/bibliography#burroughPrinciplesGeographicalInformation2015).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/plot.jl#LL24-L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.multihillshade-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.multihillshade-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.multihillshade</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
multihillshade(dem::AbstractMatrix{<:Real}; cellsize=cellsize(dem))
```


multihillshade is the simulated illumination of a surface based on its [`slope`](/reference#Geomorphometry.slope-Tuple{AbstractMatrix{<:Real}}) and [`aspect`](/reference#Geomorphometry.aspect-Tuple{AbstractMatrix{<:Real}}). Like [`hillshade`](/reference#Geomorphometry.hillshade-Tuple{AbstractMatrix{<:Real}}), but now using multiple sources as defined in [Mark, R.K. (1992)](/bibliography#mark1992multidirectional), similar to GDALs -multidirectional.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/plot.jl#LL69-L75" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.pitremoval-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.pitremoval-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.pitremoval</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pitremoval(dem::AbstractMatrix{<:Real})
```


Remove pits from a DEM Array if the center cell of a 3x3 patch is `limit` lower or than the surrounding cells.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL207-L211" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.plan_curvature-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.plan_curvature-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.plan_curvature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plan_curvature(dem::AbstractMatrix{<:Real}; cellsize = cellsize(dem), radius=1)
```


Calculate projected contour curvature (plan curvature) as defined by [Minár et al., (2020)](/bibliography#minarComprehensiveSystemDefinitions2020).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL387-L391" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.pmf-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.pmf-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.pmf</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
B, flags = pmf(A; ωₘ, slope, dhₘ, dh₀, cellsize, adjust, erode)
```


Applies the progressive morphological filter by [Zhang (2003)](/bibliography#keqizhangProgressiveMorphologicalFilter2003) to `A`.

**Output**
- `B::Array{T,2}` Maximum allowable values
  
- `flags::Array{Float64,2}` A sized array with window sizes if filtered, zero if not filtered.
  

Afterwards, one can retrieve the resulting mask for `A` by `A .<= B` or `flags .== 0.`.

**Arguments**
- `A::Array{T,2}` Input Array
  
- `ωₘ::Real=20.` Maximum window size [m]
  
- `slope::Real=0.01` Terrain slope [m/m]
  
- `dhₘ::Real=2.5` Maximum elevation threshold [m]
  
- `dh₀::Real=0.2` Initial elevation threshold [m]
  
- `cellsize::Real=1.` Cellsize in [m]
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/pmf.jl#LL1-L19" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.profile_curvature-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.profile_curvature-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.profile_curvature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
profile_curvature(dem::AbstractMatrix{<:Real}; cellsize = cellsize(dem), radius=1)
```


Calculate normal slope line curvature (profile curvature) as defined by [Minár et al., (2020)](/bibliography#minarComprehensiveSystemDefinitions2020).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL337-L341" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.pssm-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.pssm-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.pssm</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
image = pssm(dem; exaggeration=2.3, resolution=1.0)
```


Perceptually Shaded Slope Map by [Pingel, Clarke., (2014)](/bibliography#pingelPerceptuallyShadedSlope2014a).

**Output**
- `image::Gray{T,2}` Grayscale image
  

**Arguments**
- `A::Array{Real,2}` Input Array
  
- `exaggeration::Real=2.3` Factor to exaggerate elevation
  
- `cellsize::Real=1.0` Size of cell to account for horizontal resolution if different from vertical resolution
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/plot.jl#LL2-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.roughness' href='#Geomorphometry.roughness'><span class="jlbinding">Geomorphometry.roughness</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
roughness(dem::AbstractMatrix{<:Real})
```


Roughness is the largest inter-cell difference of a central pixel and its surrounding cell, as defined in Wilson et al (2007, Marine Geodesy 30:3-35).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.rugosity-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.rugosity-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.rugosity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rugosity(dem::AbstractMatrix{<:Real})
```


Compute the rugosity of a DEM, which is the ratio between the  surface area divided by the planimetric area.

Jenness 2019 https://onlinelibrary.wiley.com/doi/abs/10.2193/0091-7648%282004%29032%5B0829%3ACLSAFD%5D2.0.CO%3B2


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL161-L168" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.skb-Tuple{AbstractArray}' href='#Geomorphometry.skb-Tuple{AbstractArray}'><span class="jlbinding">Geomorphometry.skb</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
mask = skb(A; mean=mean(A))
```


Applies skewness balancing by [Bartels e.a (2006)](/bibliography#bartelsDTMGenerationLIDAR2006) to `A`. Improved the performance by applying a binary search to find the threshold value.

**Output**
- `mask::BitMatrix` Mask of allowed values
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/skew.jl#LL1-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.skbr-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.skbr-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.skbr</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
mask = skbr(A; iterations=10)
```


Applies recursive skewness balancing by [Bartels e.a (2010)](/bibliography#bartelsThresholdfreeObjectGround2010) to `A`. Applies `skb` `iterations` times to the object (non-terrain) mask, as to include more (sloped) terrain.

**Output**
- `mask::BitMatrix` Mask of allowed values
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/skew.jl#LL84-L93" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.slope-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.slope-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.slope</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
slope(dem::Matrix{<:Real}; cellsize=cellsize(dem), method=Horn(), exaggeration=1.0)
```


Slope is the rate of change between a cell and its neighbors.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL24-L28" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.smf-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.smf-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.smf</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
B = smf(A; ω, slope, dhₘ, dh₀, cellsize)
```


Applies the simple morphological filter by [Pingel et al. (2013)](/bibliography#pingelImprovedSimpleMorphological2013a) to `A`.

**Output**
- `B::Array{Float64,2}` A filtered version of A
  

**Arguments**
- `A::Array{T,2}` Input Array
  
- `ω::Float64=18.` Maximum window size [m]
  
- `slope::Float64=0.01` Terrain slope [m/m]
  
- `cellsize::Float64=1.` Cellsize in [m]
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/smf.jl#LL1-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.spread-Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}, Real}' href='#Geomorphometry.spread-Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}, Real}'><span class="jlbinding">Geomorphometry.spread</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
spread(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Real; distance=Euclidean(), res=1.0)
spread(points::Matrix{<:Real}, initial::Real, friction::Real; distance=Euclidean(), res=1.0)
```


Optimized (and more accurate) function based on the same friction everywhere.

When the friction is the same everywhere, there&#39;s no need for searching the shortest cost path, as one can just take a direct line to the input points.

The calculated cost is more accurate, as there&#39;s no &#39;zigzag&#39; from cell center to cell center.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/spread.jl#LL167-L177" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.spread-Tuple{Eastman, Vector{CartesianIndex{2}}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}' href='#Geomorphometry.spread-Tuple{Eastman, Vector{CartesianIndex{2}}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.spread</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
spread(::Eastman, points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf, iterations=3)
```


Pushbroom method for friction costs as discussed by [Eastman (1989)](/bibliography#eastman1989pushbroom). This method should scale better (linearly) than the [Tomlin (1983)](/bibliography#tomlin1983digital) method, but can require more `iterations` than set by default (3) in the case of maze-like, uncrossable obstacles.

**Output**
- `Array{Float64,2}` Total friction distance
  

**Arguments**
- `points::Matrix{<:Real}` Input Array
  
- `initial::Matrix{<:Real}` Factor to exaggerate elevation
  
- `friction::Matrix{<:Real}` Resolution of cell size
  
- `res=1` Resolution or cell size
  
- `limit=Inf` Initial fill value
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/spread.jl#LL85-L101" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.spread-Tuple{Vector{CartesianIndex{2}}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}' href='#Geomorphometry.spread-Tuple{Vector{CartesianIndex{2}}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.spread</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
spread(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; cellsize=(1,1), limit=Inf, method=Tomlin())
```


Total friction distance spread from `points` from `initial` with `friction`. By default uses Tomlin, see SpreadMethod for other algorithms.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/spread.jl#LL209-L214" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.spread-Union{Tuple{T}, Tuple{Tomlin, Vector{CartesianIndex{2}}, AbstractMatrix{T}, AbstractMatrix{<:Real}}} where T<:Real' href='#Geomorphometry.spread-Union{Tuple{T}, Tuple{Tomlin, Vector{CartesianIndex{2}}, AbstractMatrix{T}, AbstractMatrix{<:Real}}} where T<:Real'><span class="jlbinding">Geomorphometry.spread</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
spread(::Tomlin, points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf, method=Tomlin())
```


Total friction distance spread from `points` as by [Tomlin (1983)](/bibliography#tomlin1983digital). This is also the method implemented by [PCRaster](https://pcraster.geo.uu.nl/pcraster/4.0.2/doc/manual/op_spread.html).

**Output**
- `Array{Float64,2}` Total friction distance
  

**Arguments**
- `points::Matrix{<:Real}` Input Array
  
- `initial::Matrix{<:Real}` Initial values of the result
  
- `friction::Matrix{<:Real}` Resolution of cell size
  
- `res=1` Resolution or cell size
  
- `limit=Inf` Initial fill value
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/spread.jl#LL11-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.tangential_curvature-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.tangential_curvature-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.tangential_curvature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tangential_curvature(dem::AbstractMatrix{<:Real}; cellsize = cellsize(dem), radius=1)
```


Calculate normal contour curvature (tangential curvature) as defined by [Minár et al., (2020)](/bibliography#minarComprehensiveSystemDefinitions2020).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL360-L364" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Reference - Internal functions {#Reference-Internal-functions}
<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.MaximumDownwardGradient' href='#Geomorphometry.MaximumDownwardGradient'><span class="jlbinding">Geomorphometry.MaximumDownwardGradient</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Maximum Downward Gradient


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/terrain.jl#LL16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.SpreadMethod' href='#Geomorphometry.SpreadMethod'><span class="jlbinding">Geomorphometry.SpreadMethod</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Spread algorithms.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/spread.jl#LL6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.cellsize-Tuple{Any}' href='#Geomorphometry.cellsize-Tuple{Any}'><span class="jlbinding">Geomorphometry.cellsize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cellsize(dem)
```


Return an Tuple with the x and y length of each cell of the dem. Set them negatively to flip the image.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/utils.jl#LL200-L205" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.entropy!-Union{Tuple{T}, Tuple{AbstractMatrix{<:Real}, AbstractMatrix{T}}} where T<:Real' href='#Geomorphometry.entropy!-Union{Tuple{T}, Tuple{AbstractMatrix{<:Real}, AbstractMatrix{T}}} where T<:Real'><span class="jlbinding">Geomorphometry.entropy!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
entropy!(dem::AbstractMatrix{<:Real})
```


Entropy calculates the Shannon entropy of the surrounding cells of a central cell.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL107-L111" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.opening!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Integer, AbstractMatrix{T}}} where T<:Real' href='#Geomorphometry.opening!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Integer, AbstractMatrix{T}}} where T<:Real'><span class="jlbinding">Geomorphometry.opening!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Apply the opening operation to `A` with window size `ω`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/utils.jl#LL8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.opening_circ!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Integer, AbstractMatrix{T}}} where T<:Real' href='#Geomorphometry.opening_circ!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Integer, AbstractMatrix{T}}} where T<:Real'><span class="jlbinding">Geomorphometry.opening_circ!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Apply the opening operation to `A` with window size `ω`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/utils.jl#LL19" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.prominence-Tuple{AbstractMatrix{<:Real}}' href='#Geomorphometry.prominence-Tuple{AbstractMatrix{<:Real}}'><span class="jlbinding">Geomorphometry.prominence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
prominence(dem::AbstractMatrix{<:Real})
```


Prominence calculates the number of cells that are lower or equal than the central cell. Thus, 8 is a local maximum (peak), while 0 is a local minimum (pit).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/relative.jl#LL75-L80" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.round_odd-Tuple{Any}' href='#Geomorphometry.round_odd-Tuple{Any}'><span class="jlbinding">Geomorphometry.round_odd</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
round_odd(x)
```


Rounds `x` to the nearest odd number.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/733916e662f48a41cbad9468380ceae1d1f6ca6b/src/pmf.jl#LL174-L178" target="_blank" rel="noreferrer">source</a></Badge>

</details>

