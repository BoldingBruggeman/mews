import argparse
import sqlite3
import re
from typing import Optional, Dict, Any

import pandas as pd

# Pattern for FABM parameter paths in "wide" DB columns
PAT_ANY = re.compile(
    r"(?:^|[:.])fabm(?:\.yaml)?[:./]?(?:fabm[:./])?instances[:./]([^/:.]+)[:./]parameters[:./](.+)$",
    re.IGNORECASE,
)

OK_STATUSES = {"ok", "success", "done", "completed", "1", "true"}


def pick_best_row(df: pd.DataFrame, fitness: Optional[str], maximize: bool) -> pd.Series:
    df = df.copy()

    # Prefer rows that look "successful"
    if "status" in df.columns:
        ok = df["status"].astype(str).str.strip().str.lower().isin(OK_STATUSES)
        if ok.any():
            df = df.loc[ok]

    # Auto-pick a fitness column if none provided
    if fitness is None or fitness not in df.columns:
        for c in ["objective", "fitness", "lnl", "nll", "rmse", "mae"]:
            if c in df.columns:
                fitness = c
                break

    # Pick best by fitness if possible
    if fitness and fitness in df.columns:
        f = pd.to_numeric(df[fitness], errors="coerce")
        if f.notna().any():
            return df.loc[f.idxmax() if maximize else f.idxmin()]

    # Fallback: last row by a reasonable ordering
    order_cols = [c for c in ["generation", "iteration", "id", "rowid", "time", "created_at"] if c in df.columns]
    if order_cols:
        return df.sort_values(order_cols).iloc[-1]

    return df.iloc[-1]


def extract_from_wide_columns(row: pd.Series) -> Dict[str, Any]:
    """
    Extract parameters from wide-format column names like:
      fabm/instances/<inst>/parameters/<param>  (many variants)
    """
    out: Dict[str, Any] = {"instances": {}}

    for col in row.index:
        m = PAT_ANY.search(str(col))
        if not m:
            continue

        inst, pname = m.group(1), m.group(2)
        val = row[col]

        if pd.isna(val):
            continue

        # Only accept numeric values
        try:
            val_f = float(val)
        except Exception:
            continue

        out["instances"].setdefault(inst, {}).setdefault("parameters", {})[pname] = val_f

    return out


def to_4dp_float(x: float) -> float:
    """
    Round to at most 4 decimals (stored as float). Example: 0.000123 -> 0.0001
    """
    return float(f"{x:.3f}")


def apply_params_to_fabm(
    fabm_in: str,
    fabm_out: str,
    params: Dict[str, Any],
    instance: Optional[str],
    decimals: int = 3,
):
    """
    Apply extracted parameters onto fabm.yaml while preserving comments/formatting
    using ruamel.yaml round-trip mode. Also formats numbers to max N decimals.
    """
    from ruamel.yaml import YAML

    y = YAML(typ="rt")  # round-trip: preserves comments
    y.preserve_quotes = True
    # Keep a stable, readable indentation
    y.indent(mapping=2, sequence=4, offset=2)

    with open(fabm_in, "r") as f:
        base = y.load(f)

    if not isinstance(base, dict) or "instances" not in base:
        raise SystemExit("ERROR: 'instances' section missing in fabm.yaml")

    updated_any = False

    for inst, block in params.get("instances", {}).items():
        if instance and inst != instance:
            continue

        pset = block.get("parameters", {})
        if not pset:
            continue

        base_inst = base["instances"].setdefault(inst, {})
        base_params = base_inst.setdefault("parameters", {})

        for p, v in pset.items():
            # Round to max N decimals
            v = float(f"{float(v):.{decimals}f}")
            base_params[p] = v
            updated_any = True
            print(f"Updated {inst}.{p} = {v:.{decimals}f}")

    if not updated_any:
        raise SystemExit("ERROR: No matching parameters were applied (check --instance name).")

    with open(fabm_out, "w") as f:
        y.dump(base, f)

    print(f"\n✔ Wrote calibrated FABM file (comments preserved): {fabm_out}")


def main():
    ap = argparse.ArgumentParser(description="Extract & apply best parsac parameters to fabm.yaml (preserve comments)")
    ap.add_argument("db", help="SQLite database file")
    ap.add_argument("--fabm-in", required=True, help="Path to original fabm.yaml")
    ap.add_argument("--fabm-out", required=True, help="Output calibrated fabm.yaml")
    ap.add_argument("--instance", required=True, help="Instance name (e.g. Dino or selmaprotbas)")
    ap.add_argument("--fitness", default=None, help="Fitness column name (optional)")
    ap.add_argument("--maximize", action="store_true", help="Maximize fitness (default: minimize)")
    ap.add_argument("--table", default="results", help="SQLite table name (default: results)")
    ap.add_argument("--decimals", type=int, default=3, help="Max decimals to write for parameters (default: 3)")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    try:
        df = pd.read_sql(f"SELECT * FROM {args.table}", con)
    finally:
        con.close()

    if df.empty:
        raise SystemExit("No rows in results table.")

    best = pick_best_row(df, args.fitness, args.maximize)
    params = extract_from_wide_columns(best)

    if not params.get("instances"):
        raise SystemExit("ERROR: No parameters extracted from DB (column names didn’t match FABM pattern).")

    apply_params_to_fabm(
        fabm_in=args.fabm_in,
        fabm_out=args.fabm_out,
        params=params,
        instance=args.instance,
        decimals=args.decimals,
    )


if __name__ == "__main__":
    main()
