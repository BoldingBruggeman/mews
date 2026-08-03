# This code reads parameters from sensitivity_config.yaml file, run parsac sensitivity 
# and write results to sobol_indices.csv. In the code you should adjust the path,
# adjust for gotm.yaml or fabm.yaml, adjust executable name (gotm_sp1_2.exe)
# remember to delete (or rename) ..results.db file before running sensitivity.

import parsac.job.gotm
import parsac.sensitivity
import yaml
import pandas as pd
import numpy as np
import os
import sys

# ----------------------------
# Helper function
# ----------------------------
def clean_name(expr):
    """
    Replace characters that can break CSV or Pandas index
    """
    return expr.replace("[", "_").replace("]", "").replace(":", "_").replace(",", "_")

# ----------------------------
# R-style signif
# ----------------------------
def signif(x, digits=3):
    """Round to significant digits like R's signif()."""
    if x == 0:
        return 0
    return round(x, digits - int(np.floor(np.log10(abs(x)))) - 1)

# ----------------------------
# Main function
# ----------------------------
def main():
    # Paths
    folder0 = r"C:/Users/mestr/OneDrive - IOLR/MEWS"
    #folder0 = r"C:/Users/shaja/OneDrive - IOLR/MEWS"
    path_DB = os.path.join(folder0, "Git_MEWS", "mews", "kinneret", "1D")
    config_file = os.path.join(path_DB, "sensitivity_config.yaml")
    output_csv = os.path.join(path_DB, "sobol_indices.csv")

    # Load YAML configuration
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    parameters = config["parameters"]
    targets = config["targets"]
    N = config.get("sobol", {}).get("N", 256)   # default = 256 if missing

    # Setup experiment
    experiment = parsac.sensitivity.Sobol(N)
    sim = parsac.job.gotm.Simulation(".", executable="gotm_sp1_2.exe")

    # Add parameters
    for p in parameters:
        min_val = float(p["min"])
        max_val = float(p["max"])
        print(p["name"], min_val, max_val)
        experiment.add_parameter(
            sim.get_parameter("fabm.yaml", p["name"]),
            min_val,
            max_val,
            logscale=p.get("logscale", False)
        )

    # Add targets
    for t in targets:
        filename, expr = t.split(":", 1)  # split on first colon only
        sim.record_output(filename, expr)

    # Add job and run
    experiment.add_job(sim)
    p = experiment.run()

    # Prepare names for DataFrame
    param_names = [p["name"].split("/")[-1] for p in parameters]  # k_min, galp, etc.
    output_names = [clean_name(t.split(":",1)[1]) for t in targets]

    # Round values to 3 significant digits
    df = pd.DataFrame(p.T, columns=param_names, index=output_names)
    df = df.applymap(lambda v: signif(v, 3))  # 3 significant digits

    # Save CSV with first column labeled "Targets"
    df.to_csv(output_csv, index_label="Targets")
    print(f"Sobol indices saved to {output_csv}")
    print(df)


if __name__ == "__main__":
    main()
