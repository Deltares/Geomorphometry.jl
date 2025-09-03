
# Experimental {#Experimental}
> 
> [!WARNING] Methods here are experimental, not yet stable and may change or even be removed in future releases.
> 



<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.prominence' href='#Geomorphometry.prominence'><span class="jlbinding">Geomorphometry.prominence</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
prominence(dem::AbstractMatrix{<:Real})
```


Prominence calculates the number of cells that are lower or equal than the central cell. Thus, 8 is a local maximum (peak), while 0 is a local minimum (pit).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/b965e5cc3140d915fca5ca2b1f042e71a9732ba2/src/relative.jl#LL75-L80" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Geomorphometry.skbr' href='#Geomorphometry.skbr'><span class="jlbinding">Geomorphometry.skbr</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
mask = skbr(A; iterations=10)
```


Applies recursive skewness balancing by [Bartels e.a (2010)](/bibliography#bartelsThresholdfreeObjectGround2010) to `A`. Applies `skb` `iterations` times to the object (non-terrain) mask, as to include more (sloped) terrain.

**Output**
- `mask::BitMatrix` Mask of allowed values
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/Deltares/Geomorphometry.jl/blob/b965e5cc3140d915fca5ca2b1f042e71a9732ba2/src/skew.jl#LL84-L93" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```docs
Geomorphometry.pmf2
```

