#!/usr/bin/env python

from datetime import date, timedelta
import calendar
import subprocess

start_date = date(2022, 1, 1)
#stop_date = date(2022, 4, 1)
stop_date = date(2023, 1, 1)

setup="maleren"
script="run_maleren.py"

start = date(2022, 4, 1)
start = start_date
while start < stop_date:
    days_in_month = calendar.monthrange(start.year, start.month)[1]
    stop = start + timedelta(days=days_in_month)
    x = start.strftime("%Y-%m-%d %H:%M:%S")
    y = stop.strftime("%Y-%m-%d %H:%M:%S")
    print(f"Simulation from {x} to {y}")

    x = start.strftime("%Y%m%d")
    y = stop.strftime("%Y%m%d")
    if start == start_date:
        option = '--initial'
        restart_in = ''
    else:
        option = ''
        x = start.strftime("%Y%m%d")
        restart_in = f"--load_restart restart_{setup}_{x}.nc"

    restart_out = f"--save_restart restart_{setup}_{y}.nc"

    command = "mpiexec -np 4 python %s \"%s\" \"%s\" %s %s %s" % (script, start.strftime('%Y-%m-%d %H:%M:%S'), stop.strftime('%Y-%m-%d %H:%M:%S'), option, restart_in, restart_out )
    #print(command)
#    subprocess.run(["echo", command])
    subprocess.run([command], shell=True)
    subprocess.run(["mkdir", x])
    subprocess.run(["mv", "meteo.nc", "ohra_2d.nc", "ohra_3d.nc", x])
    #subprocess.run(["mv", "getm-\*.log", x])
    #print(command)

    start = stop

