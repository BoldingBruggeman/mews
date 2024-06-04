#!/usr/bin/env python

import threading
import multiprocessing
import cdsapi

c = cdsapi.Client()
dk = [
    7.75,
    58.0,
    15.25,
    54.5,
]
europa = [
    70.0,
    -20.0,
    25.0,
    45.0,
]
north_atlantic = [
    70.0,
    -45.0,
    25.0,
    45.0,
]

kinneret = [
        33.0,
        35.5,
        32.5,
        35.75,
        ]

ohra = [
        50.0,
        10.5,
        49.,
        11.,
        ]

maleren = [
        59.75,
        15.75,
        59.0,
        18.25,
        ]

area = ohra


def download_era5_year(year: int, path: str):
    request = {
        "variable": [
             "10m_u_component_of_wind",
             "10m_v_component_of_wind",
             "2m_temperature",
             "2m_dewpoint_temperature",
             "surface_pressure",
             "total_cloud_cover",
            # "total_precipitation",
        ],
        "product_type": "reanalysis",
        "format": "netcdf",
        "area": area,
        "year": "%04i" % year,
        "month": ["%02i" % m for m in range(1, 13)],
        "day": ["%02i" % d for d in range(1, 32)],
        "time": ["%02i:00" % h for h in range(0, 24)],
        "grid": ["0.25/0.25"],
    }
    r = c.retrieve("reanalysis-era5-single-levels", request, "era5_%04i.nc" % year)
    r.download(path)
    return path


start_year = 2021
stop_year = 2023
pool = multiprocessing.Pool(processes=(stop_year - start_year + 1))
results = []
for year in range(start_year, stop_year + 1):
    path = "era5_{}.nc".format(year)
    print("Queuing download of %s..." % path)
    results.append(pool.apply_async(download_era5_year, args=(year, path)))

for res in results:
    print("%s done." % res.get())
