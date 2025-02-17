module GeomorphometryRastersExt

using Geomorphometry
using Rasters, ArchGDAL
using FillArrays

degwidth::Float64 = 111_000.0

function Geomorphometry.cellsize(dem::Raster)
    T = _crstrait(dem)
    _cellsize(T, dem)
end

function _cellsize(::Rasters.GI.AbstractProjectedTrait, dem::Raster)
    dim = Rasters.dims(dem, (Rasters.XDim, Rasters.YDim))
    isintervals(dim) || throw(
        ArgumentError("Cannot calculate cell size for a `Raster` with `Points` sampling."),
    )
    (step(dim[1]), step(dim[2]))
end
function _cellsize(::Rasters.GI.AbstractGeographicTrait, dem::Raster)
    dim = Rasters.dims(dem, (Rasters.XDim, Rasters.YDim))
    Rasters.isintervals(dim) || throw(
        ArgumentError("Cannot calculate cell size for a `Raster` with `Points` sampling."),
    )
    centercoords = DimPoints(dem)[round.(Int, size(dem)[1:2] ./ 2)...]
    (step(dim[1]) * degwidth * cosd(centercoords[2]), step(dim[2]) * degwidth)
end

function _crstrait(dem::Raster)
    crs = Rasters.crs(dem)
    acrs = ArchGDAL.importCRS(crs)
    Bool(ArchGDAL.GDAL.osrisgeographic(acrs.ptr)) && return Rasters.GI.GeographicTrait()
    Bool(ArchGDAL.GDAL.osrisprojected(acrs.ptr)) && return Rasters.GI.ProjectedTrait()
    return Rasters.GI.UnknownTrait()
end

end # module
