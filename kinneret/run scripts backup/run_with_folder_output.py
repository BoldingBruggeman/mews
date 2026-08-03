from datetime import date, timedelta
from pathlib import Path
import subprocess

setup = "kinneret"
script = "run_kinneret.py"
np = 10

start_date = date(2010, 1, 1)
stop_date = date(2013, 1, 1)

start = start_date

while start < stop_date:
    stop = date(start.year + 1, 1, 1)  

    # Format date strings
    start_str = start.strftime("%Y-%m-%d %H:%M:%S")
    stop_str = stop.strftime("%Y-%m-%d %H:%M:%S")
    start_tag = start.strftime("%Y%m%d")
    stop_tag = stop.strftime("%Y%m%d")

    print(f"Running simulation for full year: {start.year}")
    command = [
        "mpiexec",
        "-n",
        str(np),
        "python",
        script,
        start_str,
        stop_str,
    ]

    # Load restart if not first year
    if start.year > 2010:
        command.extend(["--load_restart", f"restart_{setup}_{start_tag}.nc"])

    # Save restart at end of the year (Jan 1 of next year)
    command.extend(["--save_restart", f"restart_{setup}_{stop_tag}.nc"])

    # Output folder is one per year
    output_dir = Path(f"{start.year}")
    if output_dir.exists():
        print(f"Folder {output_dir} already exists - move/delete and run again")
        exit()
    else:
        output_dir.mkdir(parents=True)

    command.extend(["--output_dir", str(output_dir)])

    print("Running command:", " ".join(command))
    subprocess.run(command, check=True)

    # Advance to next year
    start = stop
