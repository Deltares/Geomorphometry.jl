"""
    Geomorphometry

Geospatial terrain analysis for digital elevation models (DEMs). Provides terrain
derivatives such as [`slope`](@ref), [`aspect`](@ref) and curvature; relative relief
metrics like [`roughness`](@ref), [`topographic_position_index`](@ref) and
[`terrain_ruggedness_index`](@ref); hydrological algorithms including
[`filldepressions`](@ref), [`flowaccumulation`](@ref) and
[`height_above_nearest_drainage`](@ref); ground filters such as
[`progressive_morphological_filter`](@ref), [`simple_morphological_filter`](@ref) and
[`skewness_balancing`](@ref); friction-distance [`spread`](@ref) algorithms; and
visualization helpers like [`pssm`](@ref) and [`hillshade`](@ref).

Functions operate on plain `AbstractMatrix` DEMs. Extensions add support for `GeoArray`
and `Raster` inputs, from which the cell size is derived automatically.
"""
module Geomorphometry
using StatsBase: skewness
using Distances: Euclidean, euclidean, evaluate
using OffsetArrays: OffsetArray
using PaddedViews: PaddedView
using FillArrays: Fill
using StaticArrays: @SMatrix, @MVector
import DataStructures
using LocalFilters: LocalFilters, dilate, localfilter!
using Stencils: Stencils, Annulus, Moore, NamedStencil, Stencil, Window, center, mapstencil
using Statistics: mean, std
using LocalFilters
using QuickHeaps: FastPriorityQueue, PriorityQueue, enqueue!, dequeue!
using KernelAbstractions: KernelAbstractions, @kernel, @index, @Const, get_backend

include("utils.jl")
include("relative.jl")
include("pmf.jl")
include("smf.jl")
include("plot.jl")
include("spread.jl")
include("terrain.jl")
include("skew.jl")
include("flowdir.jl")
include("hydrology.jl")
include("horizon.jl")

export ZevenbergenThorne, Horn, MDG
export D8, DInf, FD8
export progressive_morphological_filter, simple_morphological_filter
export pssm, hillshade, multihillshade
export depression_depth, depression_volume, drainage_potential
export pitremoval, percentile_elevation
export spread, Eastman, FastSweeping, Tomlin
export roughness,
    terrain_ruggedness_index,
    topographic_position_index,
    bathymetric_position_index,
    roughness_index_elevation,
    rugosity,
    entropy
export slope,
    aspect, laplacian, plan_curvature, profile_curvature, tangential_curvature
export skewness_balancing
export filldepressions,
    flowaccumulation,
    topographic_wetness_index,
    stream_power_index,
    height_above_nearest_drainage
export FlowDirection, LDD, D8D
export horizon_angle, sky_view_factor, viewshed, total_viewshed

end # module
