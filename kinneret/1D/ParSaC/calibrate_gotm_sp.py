import parsac.job.gotm
import parsac.optimize
from parsac.util import TextFormat

if __name__ == "__main__":
    experiment = parsac.optimize.Optimization(max_workers=8)

    sim = parsac.job.gotm.Simulation(".", executable="gotm_sp.exe")

    experiment.add_parameter(
        sim.get_parameter("gotm.yaml", "turbulence/turb_param/k_min"),
        1e-8,
        1e-5,
        logscale=True,
    )
    experiment.add_parameter(
        sim.get_parameter("gotm.yaml", "turbulence/turb_param/galp"),
        0.1,
        8,
        logscale=False,
    )
    experiment.add_target(
        sim.request_comparison(
            "output.nc",
            "temp[:,-1]",
            "./OBS/SST.obs",
            obs_file_format=TextFormat.DEPTH_INDEPENDENT,
        )
    )

    p = experiment.run(reltol=0.00001, maxgen=10)
