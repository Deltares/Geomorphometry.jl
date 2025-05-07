module GeomorphometryGeoArraysExt

using Geomorphometry
using GeoArrays
using FillArrays

degwidth::Float64 = 111_000.0

function Geomorphometry.cellsize(dem::GeoArray)
    T = GeoArrays.GI.crstrait(dem)
    _cellsize(T, dem)
end

_cellsize(::GeoArrays.GI.AbstractProjectedTrait, dem::GeoArray) =
    (dem.f.linear[1], dem.f.linear[4])
function _cellsize(::GeoArrays.GI.AbstractGeographicTrait, dem::GeoArray)
    centercoords = GeoArrays.coords(dem, round.(Int, size(dem)[1:2] ./ 2))
    (dem.f.linear[1] * degwidth * cosd(centercoords[2]), dem.f.linear[4] * degwidth)
end

function Geomorphometry.slope(
    ::Geomorphometry.GDAL,
    dem::GeoArray,
    cellsize = Geomorphometry.cellsize(dem),
    method = Geomorphometry.Horn();
    kwargs...,
)
    T = GeoArrays.GI.crstrait(dem)
    options = GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsnew(
        [
            "-alg",
            string(typeof(method)),
            "-s",
            T isa GeoArrays.GI.AbstractGeographicTrait ? "111120" : "1",
            "-of",
            "GTiff",
        ],
        C_NULL,
    )
    fn_out = tempname() * ".tif"
    GeoArrays.ArchGDAL.Dataset(dem) do ds
        ds_dempr = GeoArrays.ArchGDAL.GDAL.gdaldemprocessing(
            fn_out,
            ds.ptr,
            "slope",
            C_NULL,
            options,
            C_NULL,
        )
        GeoArrays.ArchGDAL.GDAL.gdalclose(ds_dempr)
    end
    GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsfree(options)
    GeoArrays.read(fn_out)
end

function Geomorphometry.aspect(
    ::Geomorphometry.GDAL,
    dem::GeoArray,
    cellsize = Geomorphometry.cellsize(dem),
    method = Geomorphometry.Horn();
    kwargs...,
)
    options = GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsnew(
        ["-alg", string(typeof(method)), "-of", "GTiff"],
        C_NULL,
    )
    fn_out = tempname() * ".tif"
    GeoArrays.ArchGDAL.Dataset(dem) do ds
        ds_dempr = GeoArrays.ArchGDAL.GDAL.gdaldemprocessing(
            fn_out,
            ds.ptr,
            "aspect",
            C_NULL,
            options,
            C_NULL,
        )
        GeoArrays.ArchGDAL.GDAL.gdalclose(ds_dempr)
    end
    GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsfree(options)
    GeoArrays.read(fn_out)
end

function Geomorphometry.hillshade(
    ::Geomorphometry.GDAL,
    dem::GeoArray,
    cellsize = Geomorphometry.cellsize(dem),
    method = Geomorphometry.Horn();
    kwargs...,
)
    T = GeoArrays.GI.crstrait(dem)
    options = GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsnew(
        [
            "-alg",
            string(typeof(method)),
            "-s",
            T isa GeoArrays.GI.AbstractGeographicTrait ? "111120" : "1",
            "-of",
            "GTiff",
        ],
        C_NULL,
    )
    fn_out = tempname() * ".tif"
    GeoArrays.ArchGDAL.Dataset(dem) do ds
        ds_dempr = GeoArrays.ArchGDAL.GDAL.gdaldemprocessing(
            fn_out,
            ds.ptr,
            "hillshade",
            C_NULL,
            options,
            C_NULL,
        )
        GeoArrays.ArchGDAL.GDAL.gdalclose(ds_dempr)
    end
    GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsfree(options)
    GeoArrays.read(fn_out)
end

function Geomorphometry.multihillshade(
    ::Geomorphometry.GDAL,
    dem::GeoArray;
    method = Geomorphometry.Horn(),
    kwargs...,
)
    options = GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsnew(
        ["-multidirectional", "-alg", string(typeof(method)), "-of", "GTiff"],
        C_NULL,
    )
    fn_out = tempname() * ".tif"
    GeoArrays.ArchGDAL.Dataset(dem) do ds
        ds_dempr = GeoArrays.ArchGDAL.GDAL.gdaldemprocessing(
            fn_out,
            ds.ptr,
            "hillshade",
            C_NULL,
            options,
            C_NULL,
        )
        GeoArrays.ArchGDAL.GDAL.gdalclose(ds_dempr)
    end
    GeoArrays.ArchGDAL.GDAL.gdaldemprocessingoptionsfree(options)
    GeoArrays.read(fn_out)
end

end # module
