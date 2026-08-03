# remmember to run conda activate fabmws2025
import argparse
from pathlib import Path

import yaml
import parsac.job.gotm
import parsac.optimize
from parsac.util import TextFormat

def load_config(cfg_path: Path) -> dict:
    with cfg_path.open("r") as f:
        return yaml.safe_load(f) or {}


# ----------------------------
# Parameter selection by mode
# ----------------------------
def parameter_belongs_to_mode(key: str, mode: str) -> bool:
    """
    Which PARAMETERS to include for a given mode.

    Modes:
      - all    : base + Dino + Green parameters
      - base   : base parameters only
      - do     : base parameters only (targets will be DO only)
      - po4    : base parameters only (targets will be PO4 only)
      - dino   : Dino parameters only
      - green  : Green parameters only
    """
    if mode == "all":
        return True

    if mode in ("base", "do", "po4"):
        return key.startswith("instances/selmaprotbas/")

    if mode == "dino":
        return key.startswith("instances/Dino/")

    if mode == "green":
        return key.startswith("instances/Green/")

    return False


# ----------------------------
# Target selection by mode
# ----------------------------
def add_targets(experiment: parsac.optimize.Optimization, sim, mode: str) -> None:
    """
    Which TARGETS to include for a given mode.

    Modes:
      - do     : DO only
      - po4    : PO4 only
      - base   : DO + PO4
      - dino   : Dino only
      - green  : Green only
      - all    : DO + PO4 + Dino + Green
    """

    def add_do():
        experiment.add_target(
            sim.request_comparison(
                "output.nc",
                "selmaprotbas_DO_mg",        # model: mg O2 / m³
                "./OBS/DO_mgm3.obs",         # obs:   mg O2 / m³
                obs_file_format=TextFormat.DEPTH_EXPLICIT,
            )
        )

    def add_po4():
        experiment.add_target(
            sim.request_comparison(
                "output.nc",
                "selmaprotbas_po",           # model: mmol P / m³
                "./OBS/PO4_mmolm3.obs",      # obs:   mmol P / m³
                obs_file_format=TextFormat.DEPTH_EXPLICIT,
            )
        )

    def add_dino():
        experiment.add_target(
            sim.request_comparison(
                "output.nc",
                "Dino_c",
                "./OBS/Dino.obs",
                obs_file_format=TextFormat.DEPTH_EXPLICIT,
                mindepth=-10.0,
                maxdepth=0.0,
            )
        )

    def add_green():
        experiment.add_target(
            sim.request_comparison(
                "output.nc",
                "Green_c",
                "./OBS/Green.obs",
                obs_file_format=TextFormat.DEPTH_EXPLICIT,
                mindepth=-10.0,
                maxdepth=0.0,
            )
        )

    if mode == "do":
        add_do()
        return

    if mode == "po4":
        add_po4()
        return

    if mode == "base":
        add_do()
        add_po4()
        return

    if mode == "dino":
        add_dino()
        return

    if mode == "green":
        add_green()
        return

    if mode == "all":
        add_do()
        add_po4()
        add_dino()
        add_green()
        return

    raise ValueError(f"Unknown mode: {mode}")

def build_experiment(
    mode: str,
    cfg_path: Path = Path("calibration_config.yaml"),
    executable: str = "gotm_sp1_2.exe",
    max_workers: int = 7,
) -> parsac.optimize.Optimization:
    """
    Build the parsac Optimization object for a given mode.
    """
    experiment = parsac.optimize.Optimization(max_workers=max_workers)
    sim = parsac.job.gotm.Simulation(".", executable=executable)

    # --- Load parameter definitions from YAML ---
    cfg = load_config(cfg_path)

    for p in cfg.get("parameters", []):
        key = str(p["key"])
        if not parameter_belongs_to_mode(key, mode):
            continue

        param = sim.get_parameter(p["file"], key)
        experiment.add_parameter(
            param,
            float(p["lower"]),
            float(p["upper"]),
            logscale=bool(p.get("logscale", False)),
        )

    # --- Add targets based on mode ---
    add_targets(experiment, sim, mode)

    return experiment


def main():
    parser = argparse.ArgumentParser(
        description="GOTM-FABM calibration with parsac (DO/PO4 separately or grouped; Dino/Green separately; or all)."
    )
    parser.add_argument(
        "--mode",
        choices=["all", "base", "do", "po4", "dino", "green"],
        default="all",
        help=(
            "Calibration mode:\n"
            "  do    : base params; DO target only\n"
            "  po4   : base params; PO4 target only\n"
            "  base  : base params; DO + PO4 targets\n"
            "  dino  : Dino params; Dino target only\n"
            "  green : Green params; Green target only\n"
            "  all   : base + Dino + Green params; DO + PO4 + Dino + Green targets"
        ),
    )
    parser.add_argument(
        "--reltol",
        type=float,
        default=1e-4,
        help="Relative tolerance for parsac.optimize.Optimization.run (default: 1e-4).",
    )
    parser.add_argument(
        "--maxgen",
        type=int,
        default=50,
        help="Maximum number of generations for the optimizer (default: 20).",
    )
    parser.add_argument(
        "--exe",
        type=str,
        default="gotm_sp1_2.exe",
        help="Name/path of GOTM executable (default: gotm_sp.exe).",
    )
    parser.add_argument(
        "--max_workers",
        type=int,
        default=7,
        help="Max parallel workers (default: 10).",
    )
    args = parser.parse_args()

    print(f"Building experiment in '{args.mode}' mode...")
    experiment = build_experiment(
        mode=args.mode,
        cfg_path=Path("calibration_config.yaml"),
        executable=args.exe,
        max_workers=args.max_workers,
    )

    print("Starting calibration...")
    best_params = experiment.run(reltol=args.reltol, maxgen=args.maxgen)
    print("Calibration finished.")
    print("Best parameter set:")
    print(best_params)


if __name__ == "__main__":
    main()
