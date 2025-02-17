module GeomorphometryRastersExt

using Geomorphometry
using Rasters
using FillArrays

function Geomorphometry.cellsize(dem::Raster)
    dim = Rasters.dims(dem, (Rasters.XDim, Rasters.YDim));
    isintervals(dim) || throw(ArgumentError("Cannot calculate cell size for a `Raster` with `Points` sampling."))
    xbnds, ybnds = Rasters.DD.intervalbounds(dims)
    broadcast(xb -> xb[2] - xb[1], xbnds) ./ broadcast(yb -> yb[2] - yb[1], ybnds)'
end

end # module
