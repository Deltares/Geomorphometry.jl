module Geomorphometry
using StatsBase: skewness
using OffsetArrays: OffsetMatrix
using ImageFiltering: mapwindow, sobel, imgradients
using Distances: Euclidean, euclidean, evaluate, colwise!
using PaddedViews: PaddedView
using FillArrays: Fill
using StaticArrays: @SMatrix, @MMatrix, SMatrix, MMatrix, MVector
# using DataStructures: Deque, PriorityQueue, enqueue!, dequeue!
import DataStructures
using Stencils
using Statistics: median, mean, std
# using ImageCore: scaleminmax, Gray
using LocalFilters
using QuickHeaps: FastPriorityQueue, PriorityQueue, enqueue!, dequeue!
import Eikonal

include("utils.jl")
include("relative.jl")
include("pmf.jl")
include("smf.jl")
include("psf.jl")
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
export roughness, TRI, TPI, BPI, RIE
export slope, aspect, curvature, plan_curvature, profile_curvature, tangential_curvature
export skb, skbr
export filldepressions, flowaccumulation, TWI, SPI

end # module
