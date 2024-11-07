from matplotlib import pyplot
import netCDF4
import numpy as np
import cmocean

def temperature(itime, i=20):
    with netCDF4.Dataset('kinneret_3d.nc') as nc:
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

def oxygen(itime, i=20):
    with netCDF4.Dataset('kinneret_3d.nc') as nc:
        o2 = nc['selma_o2'][itime, :, :, i]
        z = nc['zct'][itime, :, :, i]
        y = nc['yt'][:, i]
        nctime = nc['time']
        timestr = netCDF4.num2date(nctime[itime], units=nctime.units).strftime('%Y-%m-%d')
    y = np.broadcast_to(y, z.shape)

    fig, ax = pyplot.subplots()
    pc = ax.contourf(y, z, o2, 15, cmap=cmocean.cm.haline_r)
    cb = fig.colorbar(pc)
    cb.set_label('oxygen concentration (µmol/L)')
    ax.set_ylim(-45., 1.)
    ax.set_title(timestr)
    fig.savefig(f'o2_{itime:03}.png', dpi=300)
    pyplot.close(fig)

def cyanobacteria(itime, k=-1):
    with netCDF4.Dataset('kinneret_3d.nc') as nc:
        cyano = nc['cyanobacteria_c'][itime, k, :, :]
        u = nc['uk'][itime, -1, :, :]
        v = nc['vk'][itime, -1, :, :]
        x2 = nc['xt'][:, :]
        y2 = nc['yt'][:, :]
        nctime = nc['time']
        timestr = netCDF4.num2date(nctime[itime], units=nctime.units).strftime('%Y-%m-%d')

    fig, ax = pyplot.subplots()
    pc = ax.contourf(x2, y2, cyano, np.arange(0., 120.01, 10.), cmap=cmocean.cm.algae, extend='max')
    q = ax.quiver(x2, y2, u, v, scale=2)
    path = f'cyano_{itime:03}.png'
    print(f'Saving {path}')
    ax.quiverkey(q, 0.25, 0.1, .1, '0.25 m/s', labelpos='W')
    cb = fig.colorbar(pc, ax=ax)
    cb.set_label('cyanobacteria (µmol C/L)')
    ax.set_title(timestr)
    ax.axis('equal')
    fig.savefig(path, dpi=96)
    pyplot.close(fig)

temperature(105)
oxygen(105)
#for itime in range(2*365):
#    cyanobacteria(itime)