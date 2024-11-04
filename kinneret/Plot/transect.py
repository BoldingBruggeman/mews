from matplotlib import pyplot
import netCDF4
import numpy as np
import cmocean

import xarray as xr

def temperature(itime, i=20):
    with netCDF4.Dataset('../20220101/kinneret_3d.nc') as nc:
    #with xr.open_mfdataset("20220101/kinneret_3d.nc") as nc:
        temp = nc['temp'][itime, :, :, i]
        z = nc['zct'][itime, :, :, i]
        y = nc['yt'][:, i]
        nctime = nc['time']
        timestr = netCDF4.num2date(nctime[itime], units=nctime.units).strftime('%Y-%m-%d')
    y = np.broadcast_to(y, z.shape)

    fig, ax = pyplot.subplots()
    pc = ax.contourf(y, z, temp, cmap=cmocean.cm.thermal)
    cb = fig.colorbar(pc)
    cb.set_label('temperature (°C)')
    ax.set_ylim(-45., 1.)
    ax.set_title(timestr)
    fig.savefig(f'temp_{itime:03}.png', dpi=300)
    pyplot.close(fig)


for itime in range(33):
    temperature(itime)
