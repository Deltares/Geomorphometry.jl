module GeomorphometryRastersExt

using Geomorphometry
using Rasters
using FillArrays

function Geomorphometry.xyratio(dem::Raster)
    return Fill(1.0, size(dem))
end

end # module
