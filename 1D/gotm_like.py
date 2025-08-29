import datetime

import numpy as np

import pygetm
import pygetm.input.igotm

setup = "ohra"
setup = "kinneret"
setup = "ekoln"

# setup specifics
if setup == "ekoln":
    lon, lat, depth = 17.61, 59.75, 30.0
    timestep, split_factor = 360.0, 10
    tinit = 0.1
    sinit = 0.1
    output = "ekoln.nc"

if setup == "kinneret":
    lon, lat, depth = 35.6, 32.8, 45.0
    timestep, split_factor = 360.0, 10
    tinit = 16.0
    sinit = 0.1
    output = "kinneret.nc"

if setup == "ohra":
    lon, lat, depth = 10.7, 50.76, 42.0
    timestep, split_factor = 360.0, 10
    tinit = 0.1
    sinit = 0.1
    output = "ohra.nc"


jerlov = True
if jerlov:
    jerlov_type = pygetm.Jerlov.Type_III
else:
    A = 0.78
    kc1 = 1.0 / 1.5
    kc2 = 1.0 / 40.0

x = np.linspace(-1000, 1000, 3)
y = np.linspace(-1000, 1000, 3)
domain = pygetm.domain.create_cartesian(
    x,
    y,
    interfaces=True,
    H=depth,
    periodic_x=True,
    periodic_y=True,
    lat=lat,
    lon=lon,
    z0=0.03,
)

if jerlov:
    radiation = pygetm.radiation.TwoBand(jerlov_type=jerlov_type)
else:
    radiation = pygetm.radiation.TwoBand()

airsea = pygetm.airsea.FluxesFromMeteo(
    humidity_measure=pygetm.HumidityMeasure.DEW_POINT_TEMPERATURE,
)
vertical_coordinates = pygetm.vertical_coordinates.Sigma(int(depth))
sim = pygetm.Simulation(
    domain,
    airsea=airsea,
    radiation=radiation,
    vertical_coordinates=vertical_coordinates,
    gotm="gotm.yaml",
)
if not jerlov:
    sim.radiation.A.fill(A)
    sim.radiation.kc1.fill(kc1)
    sim.radiation.kc2.fill(kc2)
    sim.logger.info(
        f"User provided attunuation:A={A:.2f}, kc1={kc1:.4f} m-1, kc2={kc2:.4f} m-1"
    )
    quit()

era5 = pygetm.input.igotm.download_era5(lon, lat, 2020, logger=sim.logger)
sim.airsea.t2m.set(era5["t2m"])
sim.airsea.d2m.set(era5["d2m"])
sim.airsea.u10.set(era5["u10"])
sim.airsea.v10.set(era5["v10"])
sim.airsea.sp.set(era5["sp"])
sim.airsea.tcc.set(era5["tcc"])
start = era5.time.values[0]
stop = era5.time.values[-1]

sim.temp.set(tinit)
sim.salt.set(sinit)

out = sim.output_manager.add_netcdf_file(
    output,
    interval=datetime.timedelta(hours=1),
    sync_interval=None,
)

out.request("temp", "salt", "sst", "uk", "vk", "SS", "NN", "rad", "nuh")
sim.start(
    start,
    timestep,
    split_factor,
    report=datetime.timedelta(hours=24),
    report_totals=datetime.timedelta(days=10),
)
while sim.time < stop:
    sim.advance()
sim.finish()
