#!/usr/bin/env python

from datetime import date, timedelta
import calendar
from pathlib import Path
import subprocess

setup = "ohra"
script = "run_ohra.py"
np = 4

start_date = date(2022, 1, 1)
stop_date = date(2023, 1, 1)

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

    command.extend([f"--output_dir", f"{x}"])

    if True:
        print(command)
    else:
        subprocess.run(command)

    start = stop
