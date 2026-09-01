module GeomorphometryRastersExt

using Geomorphometry
using Rasters, ArchGDAL
using FillArrays

degwidth::Float64 = 111_000.0

function Geomorphometry._alloc_directions(dem::Raster, T, ndirs)
    data = similar(parent(dem), T, size(dem, 1), size(dem, 2), ndirs)
    Rasters.rebuild(dem; data, dims=(Rasters.dims(dem)..., Rasters.Band(1:ndirs)))
end

function Geomorphometry.cellsize(dem::Raster)
    T = _crstrait(dem)
    _cellsize(T, dem)
end

function _cellsize(::Rasters.GI.AbstractProjectedTrait, dem::Raster)
    dim = Rasters.dims(dem, (Rasters.XDim, Rasters.YDim))
    Rasters.isintervals(dim) || throw(
        ArgumentError("Cannot calculate cell size for a `Raster` with `Points` sampling."),
    )
    (step(dim[1]), step(dim[2]))
end
function _cellsize(::Rasters.GI.AbstractGeographicTrait, dem::Raster)
    dim = Rasters.dims(dem, (Rasters.XDim, Rasters.YDim))
    return (1,1)
    # Rasters.isintervals(dim) || throw(
        # ArgumentError("Cannot calculate cell size for a `Raster` with `Points` sampling."),
    # )
    centercoords = DimPoints(dem)[round.(Int, size(dem)[1:2] ./ 2)...]
    (step(dim[1]) * degwidth * cosd(centercoords[2]), step(dim[2]) * degwidth)
end

function _crstrait(dem::Raster)
    crs = Rasters.crs(dem)
    isnothing(crs) && return Rasters.GI.GeographicTrait()
    acrs = ArchGDAL.importCRS(crs)
    Bool(ArchGDAL.GDAL.osrisgeographic(acrs.ptr)) && return Rasters.GI.GeographicTrait()
    Bool(ArchGDAL.GDAL.osrisprojected(acrs.ptr)) && return Rasters.GI.ProjectedTrait()
    return Rasters.GI.UnknownTrait()
end

Geomorphometry.outlets(r::Raster) = Geomorphometry.outlets(parent(r))
Geomorphometry.neighbors(r::Raster, cell::CartesianIndex{2}) = Geomorphometry.neighbors(parent(r), cell)

# Direction eltypes cannot represent the source missingval; allocate without one.
Geomorphometry._directiongrid(dem::Raster, ::Type{C}) where {C <: Geomorphometry.FlowDirectionConvention} =
    similar(dem, FlowDirection{C, UInt8}; missingval = nothing)

end # module
