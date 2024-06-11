from typing import Optional
import numpy as np
from pyproj import Transformer
import xarray as xr

try:
    import rioxarray
except ImportError:
    raise Exception("You need rioxarray. See https://corteva.github.io/rioxarray")


def handle_bathymetry(
    tif_file: str,
    verbose: bool = False,
) -> xr.DataArray:
    xda = rioxarray.open_rasterio(
        tif_file,
        engine="rasterio",
        default_name="bathymetry",
    )

    xda = -1.0 * xda
    xda = xda.fillna(-9999.0)
    xda.attrs["_FillValue"] = -9999.0

    if xda.ndim == 3 and xda.shape[0] == 1:
        # Squeeze out band dimension
        xda = xda[0, :, :]

    if verbose:
        print(f"{xda.rio.crs}")
        x = xda.coords[xda.rio.x_dim]
        y = xda.coords[xda.rio.y_dim]
        print(f"Grid size (nx, ny) = ({len(x)}, {len(y)})")
        print(f"X-axis: {x.min():.7} - {x.max():.7} m")
        print(f"Y-axis: {y.min():.7} - {y.max():.7} m")

    return xda


def create_lon_lat_grid(
    xds: xr.DataArray,
    epsg: str,
    verbose: bool = False,
):
    # -> xr.DataArray, xr.DataArray:

    to_crs = epsg
    if verbose:
        print(to_crs)

    transformer = Transformer.from_crs(xds.rio.crs, to_crs, always_xy=True)
    x_coords, y_coords = np.meshgrid(
        xds.coords[xds.rio.x_dim], xds.coords[xds.rio.y_dim]
    )
    lon, lat = transformer.transform(x_coords, y_coords)

    lon_xr = xr.DataArray(
        lon,
        name="lon",
        attrs=dict(
            units="degrees_east",
            description="lon",
        ),
        coords=xds.coords,
        dims=["y", "x"],
    )

    lat_xr = xr.DataArray(
        lat,
        name="lat",
        attrs=dict(
            units="degrees_north",
            description="lat",
        ),
        coords=xds.coords,
        dims=["y", "x"],
    )

    if xds.ndim == 3 and xds.shape[0] == 1:
        lon_xr = lon_xr[0, :, :]
        lat_xr = lat_xr[0, :, :]

    if verbose:
        print(f"Longitude: {lon.min():.7} - {lon.max():.7} degrees East")
        print(f"Latitude:  {lat.min():.7} - {lat.max():.7} degrees North")

    return lon_xr, lat_xr


def convert_ds_netcdf(
    components: list,
    nc_path: str,
    verbose: bool = False,
):
    ds = xr.merge(components)
    ds.to_netcdf(nc_path)
    if verbose:
        ds.info()
    return ds


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--tif_file",
        nargs="?",
        help="GeoTIFF file to read",
        default="grid_malaren_depth_epsg3006masked.tif",
    )
    parser.add_argument("--epsg", nargs="?", help="EPGS code", default="EPSG:4326")
    parser.add_argument(
        "--nc_path", nargs="?", help="Output NetCDF file", default="bathymetry.nc"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="debug output")
    args = parser.parse_args()

    bathy = handle_bathymetry(args.tif_file, verbose=args.verbose)

    lon, lat = create_lon_lat_grid(bathy, args.epsg, verbose=args.verbose)

    ds = convert_ds_netcdf([bathy, lon, lat], args.nc_path, verbose=args.verbose)


# https://stackoverflow.com/questions/58128824/calculationg-lat-lon-from-geostationary-projection
# https://pyproj4.github.io/pyproj/stable/gotchas.html#proj-not-a-generic-latitude-longitude-to-projection-converter
# from osgeo import gdal
# inputfile='grid_malaren_depth_epsg3006masked.tif'
# outputfile='grid_malaren_depth_epsg3006masked.nc'
# ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
