#!/usr/bin/env python

from datetime import date, timedelta
import calendar
from pathlib import Path
import subprocess

start_date = date(1995, 1, 1)
# stop_date = date(2022, 4, 1)
stop_date = date(1995, 3, 31)

setup = "malaren"
script = "run_malaren.py"
np = 25

start = date(1995, 1, 1)
start = start_date

while start < stop_date:
    command = [
        "mpiexec",
        "-n",
        str(np),
        "python",
        script,
    ]
    
    days_in_month = calendar.monthrange(start.year, start.month)[1]
    stop = start + timedelta(days=days_in_month)
    x = start.strftime("%Y-%m-%d %H:%M:%S")
    y = stop.strftime("%Y-%m-%d %H:%M:%S")
    print(f"Simulation from {x} to {y}")
    
    command.extend([x, y])
    
    x = start.strftime("%Y%m%d")
    y = stop.strftime("%Y%m%d")
    if start == start_date:
        print("Initial simulation")
    else:
        command.extend([f"--load_restart", f"restart_{setup}_{x}.nc"])
    
    
    command.extend([f"--save_restart", f"restart_{setup}_{y}.nc"])

    p = Path(x)
    if p.is_dir():
        print(f"Folder {x} already exists - move/delete and run again")
        exit()
    else:
        p.mkdir(parents=True)
    
    # command.extend(["--output_dir", "x"])
    command.extend([f"--output_dir", f"{x}"])
    
    subprocess.run(command, shell=True)
    
    # if True:
        # print(command)
    # else:
        # subprocess.run(command)
    start = stop
