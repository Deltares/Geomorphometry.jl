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
using Statistics: median, mean
using ImageCore: scaleminmax, Gray
using LocalFilters
using QuickHeaps: FastPriorityQueue, PriorityQueue, enqueue!, dequeue!

include("utils.jl")
include("pmf.jl")
include("smf.jl")
include("psf.jl")
include("plot.jl")
include("spread.jl")
include("terrain.jl")
include("skew.jl")
include("hydrology.jl")

export ZevenbergenThorne, Horn, MDG
export D8, Dinf
export pmf, smf, psf
export pssm
export pitremoval
export spread, spread2
export roughness, TRI, TPI, BPI, slope, aspect, curvature, hillshade, multihillshade
export skb, skbr
export priorityflood, streamflow

end # module
