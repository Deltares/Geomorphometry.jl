module Geomorphometry
using StatsBase: skewness
using Distances: Euclidean, euclidean, evaluate
using OffsetArrays: centered
using PaddedViews: PaddedView
using FillArrays: Fill
using StaticArrays: @SMatrix, @MMatrix, @MVector
import DataStructures
using LocalFilters: LocalFilters, dilate, localfilter!
using Stencils: Stencils, Annulus, Moore, NamedStencil, Stencil, Window, center, mapstencil
using Statistics: mean, std
using LocalFilters
using QuickHeaps: FastPriorityQueue, PriorityQueue, enqueue!, dequeue!

include("utils.jl")
include("relative.jl")
include("pmf.jl")
include("smf.jl")
include("plot.jl")
include("spread.jl")
include("terrain.jl")
include("skew.jl")
include("hydrology.jl")

export ZevenbergenThorne, Horn, MDG
export D8, DInf, FD8
export pmf, smf, psf
export pssm, hillshade, multihillshade
export pitremoval
export spread, Eastman, FastSweeping, Tomlin
export roughness, TRI, TPI, BPI, RIE, rugosity, entropy
export slope, aspect, curvature, plan_curvature, profile_curvature, tangential_curvature
export skb, skbr
export filldepressions, flowaccumulation, TWI, SPI

end # module
