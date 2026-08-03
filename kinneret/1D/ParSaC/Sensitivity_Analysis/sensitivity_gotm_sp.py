import parsac.job.gotm
import parsac.sensitivity
import pandas as pd

if __name__ == "__main__":
    # experiment = parsac.sensitivity.MVR(20)
    # experiment = parsac.sensitivity.Morris(8)
    # experiment = parsac.sensitivity.CV(100)
    experiment = parsac.sensitivity.Sobol(8)

    sim = parsac.job.gotm.Simulation(".", executable="gotm_sp.exe")

    experiment.add_parameter(
        sim.get_parameter("gotm.yaml", "turbulence/turb_param/k_min"),
        1e-7,
        1e-5,
        logscale=True,
    )
    experiment.add_parameter(
        sim.get_parameter("gotm.yaml", "turbulence/turb_param/galp"),
        0.1,
        8,
        logscale=False,
    )
    # Targets for sensitivity analysis
    sim.record_output("output.nc", "temp[:,0].max()")
    sim.record_output("output.nc", "temp[:,0].min()")
    sim.record_output("output.nc", "temp[:,-1].max()")
    sim.record_output("output.nc", "temp[:,-1].min()")

    experiment.add_job(sim)

    p = experiment.run()

    print(type(p))
    print(p.shape)
    print(p)

    # Your parameter names in the same order as defined
    param_names = ["k_min", "galp"]
    
    # Your output targets in the same order as you added them
    output_names = [
        "temp[:,0].max()",
        "temp[:,0].min()",
        "temp[:,-1].max()",
        "temp[:,-1].min()"
    ]
    
    # Reshape p: currently (2,4) = params × outputs
    # Transpose to (4,2) = outputs × params
    df = pd.DataFrame(p.T, columns=param_names, index=output_names)
    
    # Save to CSV
    df.to_csv("sobol_indices.csv")
    print("Saved sensitivity indices to sobol_indices.csv")
    print(df)