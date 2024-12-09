#!/usr/bin/env python

from datetime import date, timedelta
import calendar
from pathlib import Path
import subprocess

start_date = date(1995, 1, 1)
# stop_date = date(2022, 4, 1)
stop_date = date(1995, 3, 31)

setup = "ekoln"
script = "run_ekoln.py"
np = 16

start = date(1995, 1, 1)
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
        option = "--initial"
        restart_in = ""
    else:
        option = ""
        restart_in = f"--load_restart restart_{setup}_{x}.nc"

    restart_out = f"--save_restart restart_{setup}_{y}.nc"

    p = Path(x)
    if p.is_dir():
        print(f"Folder {x} already exists - move/delete and run again")
        exit()
    else:
        p.mkdir(parents=True)
    output_dir = f"--output_dir {x}"

    command = 'mpiexec -np %d python %s "%s" "%s" %s %s %s %s' % (
        np,
        script,
        start.strftime("%Y-%m-%d %H:%M:%S"),
        stop.strftime("%Y-%m-%d %H:%M:%S"),
        option,
        restart_in,
        restart_out,
        output_dir,
    )
    
    # Full command to run inside Anaconda
    # full_command = f'cmd.exe /K "{activate_script} && conda activate {env_name} && {command}"'
    
    subprocess.run(command, shell=True)
    # subprocess.run(full_command, shell=True)
    # subprocess.run(["mv", "getm-\*.log", x])

    start = stop
