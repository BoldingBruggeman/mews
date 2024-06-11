from gettext import npgettext
import scipy.interpolate
import numpy as np
import xarray as xr

# CRS.from_epsg(4326).to_cf() from pyproj import CRS
# https://py.geocompx.org/
# https://corteva.github.io/rioxarray/stable/readme.html

fname = "bathyfinal-4326.csv"
fname = "bathyfinal-3397.csv"

CARTESIAN = 1
SPHERICAL = 2

grid_type = CARTESIAN


def read_data(fname: str):
    print(f"Reading data from {fname}")
    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    x, y, z = dat[:, 0], dat[:, 1], dat[:, 2]

    if -360.0 < x[0] and x[0] < 360.0:
        grid_type = SPHERICAL
    else:
        grid_type = CARTESIAN

    return x, y, z


def interp_grid(x, y, z, xlen=50, ylen=50, zref=520, method="linear"):
    from matplotlib import pyplot

    print(f"Interpolating irregular data to regular grid")

    xgrid = np.linspace(x.min(), x.max(), xlen)
    ygrid = np.linspace(y.min(), y.max(), ylen)

    zgrid = scipy.interpolate.griddata(
        (x, y), z, (xgrid[np.newaxis, :], ygrid[:, np.newaxis]), method=method
    )

    # define surface as reference
    zgrid = zgrid - zref

    # mask out land
    zgrid[zgrid > 0.0] = np.nan

    # GETM require depths to be positive
    zgrid = -1.0 * zgrid

    # convert to an XArray data array
    if grid_type == CARTESIAN:
        data_xr = xr.DataArray(
            zgrid,
            name="bathymetry",
            attrs=dict(
                units="m",
                description="depth data",
                interpolation=method,
                project="MEWS",
            ),
            coords=dict(
                x=(["x"], xgrid, {"units": "m"}),
                y=(["y"], ygrid, {"units": "m"}),
            ),
            dims=["y", "x"],
        )
    else:
        data_xr = xr.DataArray(
            zgrid,
            name="bathymetry",
            attrs=dict(
                units="m",
                description="depth data",
                interGpolation=method,
                project="MEWS",
            ),
            coords=dict(
                lon=(["lon"], xgrid, {"units": "degrees_east"}),
                lat=(["lat"], ygrid, {"units": "degrees_north"}),
            ),
            dims=["lat", "lon"],
        )

    data_xr = data_xr.fillna(-9999.0)
    data_xr.attrs["_FillValue"] = -9999.0

    fig, ax = pyplot.subplots()
    pc = ax.pcolormesh(xgrid, ygrid, zgrid)
    cb = fig.colorbar(pc)
    if grid_type == CARTESIAN:
        fig.savefig("bathymetry_cartesian.png")
    else:
        fig.savefig("bathymetry_spherical.png")

    return data_xr


def create_lon_lat(
    xda: xr.DataArray,
    epsg: str,
    verbose: bool = False,
):
    from pyproj import Transformer

    print(f"Creating lon, lat grids:")

    to_crs = epsg
    if verbose:
        print(to_crs)

    transformer = Transformer.from_crs("EPSG:3397", to_crs, always_xy=True)
    x_coords, y_coords = np.meshgrid(xda.coords["x"], xda.coords["y"])
    lon, lat = transformer.transform(x_coords, y_coords)

    lon_xr = xr.DataArray(
        lon,
        name="lon",
        attrs=dict(
            units="degrees_east",
            description="lon",
        ),
        coords=xda.coords,
        dims=["y", "x"],
    )

    lat_xr = xr.DataArray(
        lat,
        name="lat",
        attrs=dict(
            units="degrees_north",
            description="lat",
        ),
        coords=xda.coords,
        dims=["y", "x"],
    )

    if xda.ndim == 3 and xda.shape[0] == 1:
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


def print_summary():

    print(f"Summary:")
    print(f"(xlen,ylen) = ({xlen}, {ylen})")
    if grid_type == CARTESIAN:
        print(f"(xmin,xmax) = ({x.min():.2f}, {x.max():.2f})")
        print(f"(ymin,ymax) = ({y.min():.2f}, {y.max():.2f})")
        dx = (x.max() - x.min()) / xlen
        dy = (y.max() - y.min()) / ylen
        print(f"(dx,dy) = ({dx:.2f}, {dy:.2f})")
        data_xr.to_netcdf("bathymetry_cartesian.nc")
    else:
        print(f"(xmin,xmax) = ({x.min():.4f}, {x.max():.4f})")
        print(f"(ymin,ymax) = ({y.min():.4f}, {y.max():.4f})")
        dx = (x.max() - x.min()) / xlen
        dy = (y.max() - y.min()) / ylen
        print(f"(dx,dy) =     ({dx:.5e}, {dy:.5e})")
        data_xr.to_netcdf("bathymetry_spherical.nc")



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--data_file",
        nargs="?",
        help="CSV separated data file to read",
        default="bathyfinal-3397.csv",
    )
    parser.add_argument("--epsg", nargs="?", help="EPGS code", default="EPSG:4326")
    parser.add_argument(
        "--nc_path", nargs="?", help="Output NetCDF file", default="bathymetry.nc"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="debug output")
    args = parser.parse_args()

    x, y, z = read_data(fname)

    xda = interp_grid(x, y, z)

    if grid_type == CARTESIAN:
        xlon, xlat = create_lon_lat(xda, args.epsg, args.verbose)
        components = [xda, xlon, xlat]
    else:
        components = [xda]

    convert_ds_netcdf(components, args.nc_path)

    # print_summary()
