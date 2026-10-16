#!/usr/bin/env python3
"""Regenerate v4 initial states from recovered Group I parameters."""
import argparse
import csv
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", type=int, help="Run one state; default: all missing states")
    args = parser.parse_args()
    from iloscar_backend import parameters, init_start, model_run
    base = Path(__file__).resolve().parent
    os.chdir(base)
    states = base / "steady_state_fsed"
    states.mkdir(exist_ok=True)
    with open("mc_accepted_runs.csv") as f:
        rows = list(csv.DictReader(f))
    if args.run_id is not None and not any(int(r["run_id"]) == args.run_id for r in rows):
        parser.error("run-id is not in the supplied design")
    for row in rows:
        run_id = int(row["run_id"])
        if args.run_id is not None and run_id != args.run_id:
            continue
        state = states / f"mc_{run_id}_fsed_steady.dat"
        if state.exists():
            continue
        p = parameters(pcycle=0, FTYS=1, printflag=1)
        p.update(RUN_TYPE=1, LOADFLAG=0, coeff_flag=2,
                 exp_name=f"spinup_{run_id}", svresults=0, svstart=1,
                 svfile=str(state), tfinal=2e7, cinpflag=0,
                 cac=0.02, mgc=0.03, d13cin=2, d13cvc=-1.5,
                 pref=float(row["pCO2_init"]), pcsi=float(row["pCO2_init"]),
                 fsi0=float(row["silicate_weathering_init"]),
                 finc0=float(row["carbonate_weathering_init"]),
                 rrain=float(row["rain_ratio"]),
                 fsh=float(row["shelf_basin_calcite_burial"]),
                 thc=float(row["thc_flux"]), po4_scaler=float(row["po4_scaler"]))
        init_start(p)
        model_run()
        print(f"Saved {state}")


if __name__ == "__main__":
    main()
