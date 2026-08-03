import os, sys
import parsac.job.gotm as gotm
from parsac.util import TextFormat
from parsac.sensitivity import Sobol  # or: from parsac.sensitivity import Morris

WORKDIR   = "."
EXECUTABLE = "gotm_sp.exe"
GOTM_YAML = "gotm.yaml"
OBS_FILE  = r".\OBS\SST.obs"

SOBOL_NSAMPLES = 200  # start small; raise later

def main():
    sim = gotm.Simulation(WORKDIR, executable=EXECUTABLE)

    # Attach target directly to the Simulation (newer API)
    sim.request_comparison(
        variable="temperature",
        filename=OBS_FILE,
        fmt=TextFormat.DEPTH_INDEPENDENT,
    )

    exp = Sobol(sim, nsamples=SOBOL_NSAMPLES)
    # For screening:
    # from parsac.sensitivity import Morris
    # exp = Morris(sim, ntraj=15, levels=8)

    exp.add_parameter(sim.get_parameter(GOTM_YAML, "turbulence/turb_param/k_min"),
                      1e-8, 1e-5, logscale=True)
    exp.add_parameter(sim.get_parameter(GOTM_YAML, "turbulence/turb_param/galp"),
                      0.1, 8.0, logscale=False)

    print("[parsac] Starting sensitivity run...")
    exp.run()
    db = os.path.splitext(os.path.basename(__file__))[0] + ".results.db"
    print("[parsac] Done. Results in:", db)
    print("Plot:", "python -m parsac.optimize.plot", db)
    print("Dump:", "python -m parsac.record", db, "--dump results.txt")

if __name__ == "__main__":
    main()
