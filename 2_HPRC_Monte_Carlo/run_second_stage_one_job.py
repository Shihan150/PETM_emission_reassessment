import os
import traceback
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from iloscar_backend import parameters, init_start, model_run


# =========================
# Settings
# =========================
BATCH_SIZE = 50 

group1_file = "mc_accepted_runs.csv"
group2_file = "mc_group2_samples.csv"

steady_dir = Path(__file__).resolve().parent / "steady_state_fsed"
output_dir = Path(__file__).resolve().parent / "generated_results"
output_dir.mkdir(exist_ok=True)

target_file = "mc_inv_ph.csv"
BASE_DIR = Path(__file__).resolve().parent
os.chdir(BASE_DIR)


# =========================
# Read parameter files
# =========================
df_group1 = pd.read_csv(group1_file)
df_group1 = df_group1[df_group1["accepted"] == True].reset_index(drop=True)

df_group2 = pd.read_csv(group2_file).reset_index(drop=True)

n_group1 = len(df_group1)
n_group2 = len(df_group2)

TOTAL_RUNS = n_group1 * n_group2
N_BATCHES = int(np.ceil(TOTAL_RUNS / BATCH_SIZE))

print(f"Accepted Group I runs: {n_group1}")
print(f"Group II samples: {n_group2}")
print(f"Total inversions: {TOTAL_RUNS}")
print(f"Batch size: {BATCH_SIZE}")
print(f"Number of batches: {N_BATCHES}")


# =========================
# Read SLURM array task ID
# =========================
task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])

start_id = task_id * BATCH_SIZE
end_id = min(start_id + BATCH_SIZE, TOTAL_RUNS)

if start_id >= TOTAL_RUNS:
    raise ValueError(
        f"task_id={task_id} starts at {start_id}, "
        f"but TOTAL_RUNS={TOTAL_RUNS}"
    )

print(f"Running task_id={task_id}")
print(f"Global IDs: {start_id} to {end_id - 1}")


results = []
failed = []


# =========================
# Batch loop
# =========================
for global_id in range(start_id, end_id):

    group1_index = global_id // n_group2
    group2_index = global_id % n_group2

    row1 = df_group1.iloc[group1_index]
    row2 = df_group2.iloc[group2_index]

    group1_id = int(row1["run_id"])
    sample_id = int(row2["sample_id"])

    steady_file = steady_dir / f"mc_{group1_id}_fsed_steady.dat"

    exp_name = f"MC_G1_{group1_id}_G2_{sample_id}"

    print("=" * 60)
    print(f"Running global_id={global_id}")
    print(f"exp_name={exp_name}")
    print(f"steady_file={steady_file}")

    work_dir = Path(tempfile.mkdtemp(prefix=f"petm_{global_id}_"))
    try:
        if not steady_file.exists():
            raise FileNotFoundError(f"Missing steady-state file: {steady_file}")

        shutil.copy2(BASE_DIR / "myami_ca_mg_coeff_table.npz", work_dir)
        target = pd.read_csv(BASE_DIR / target_file)
        target.iloc[:, 1] *= float(row2["ph_drop_target"])
        target.to_csv(work_dir / "target_ph.csv", index=False)
        os.chdir(work_dir)

        p = parameters(
            pcycle=0,
            LOADFLAG=1,
            FTYS=1,
            printflag=0
        )

        # =========================
        # Load accepted steady state
        # =========================
        p["id"] = global_id
        p["RUN_TYPE"] = 2
        p["exp_name"] = exp_name
        p["LOADFLAG"] = 1
        p["initfile"] = str(steady_file)

        # =========================
        # Fixed parameters
        # =========================
        p["cac"] = 0.02
        p["mgc"] = 0.03
        p["d13cin"] = 2
        p["d13cvc"] = -1.5
        p["coeff_flag"] = 2

        # =========================
        # Accepted Group I parameters
        # =========================
        p["pref"] = row1["pCO2_init"]
        p["pcsi"] = row1["pCO2_init"]

        p["fsi0"] = row1["silicate_weathering_init"]
        p["finc0"] = row1["carbonate_weathering_init"]

        p["rrain"] = row1["rain_ratio"]
        p["fsh"] = row1["shelf_basin_calcite_burial"]

        p["thc"] = row1["thc_flux"]
        p["po4_scaler"] = row1["po4_scaler"]
        p["atlantic_extra_stagnation"] = row2["atlantic_extra_stagnation"]

        # =========================
        # Sampled Group II parameters
        # =========================
        p["nsi"] = row2["nsi"]
        p["ncc"] = row2["ncc"]

        p["sclim"] = row2["climate_sensitivity"]
        p["thcweak"] = row2["thc_weakening"]
        p["r_deep_atlantic"] = row2["deep_atlantic_carbon_injection"]

        # =========================
        # Save transient results
        # =========================
        p["svresults"] = 1
        p["svstart"] = 0

        # =========================
        # Inversion settings
        # =========================
        p["target"] = "GSpH"
        p["target_file"] = "target_ph.csv"
        p["toms_low"] = -0.2
        p["toms_high"] = 1

        p["tfinal"] = 5e4

        init_start(p)

        (
            emi_onset,
            emi_total,
            wt_drop,
            pacific_ccd_shoaling,
            atlantic_ccd_shoaling,
            ph_drop,
            rco2,
            delta_t
        ) = model_run()

        ccd = pd.read_csv(Path(exp_name) / "ccd.csv")
        fcp = pd.read_csv(Path(exp_name) / "fcp.csv")
        results.append({
            "ph_drop_target": row2["ph_drop_target"],
            "ph_drop_achieved": ph_drop,
            "ph_drop_error": ph_drop - row2["ph_drop_target"],
            "po4_scaler": row1["po4_scaler"],
            "atlantic_extra_stagnation": row2["atlantic_extra_stagnation"],
            "pacific_ccd_min": ccd["Pacific"].min(),
            "atlantic_ccd_min": ccd["Atlantic"].min(),
            "pacific_2900_min": fcp["2900.0"].min(),
            "pacific_3200_min": fcp["3200.0"].min(),
            "global_id": global_id,
            "task_id": task_id,
            "group1_index": group1_index,
            "group2_index": group2_index,
            "group1_id": group1_id,
            "sample_id": sample_id,
            "status": "success",

            "emi_onset": emi_onset,
            "emi_total": emi_total,
            "wt_drop": wt_drop,
            "pacific_ccd_shoaling": pacific_ccd_shoaling,
            "atlantic_ccd_shoaling": atlantic_ccd_shoaling,
            "ph_drop": ph_drop,
            "RCO2": rco2,
            "delta_T": delta_t,

            "pCO2_init": row1["pCO2_init"],
            "silicate_weathering_init": row1["silicate_weathering_init"],
            "carbonate_weathering_init": row1["carbonate_weathering_init"],
            "rain_ratio": row1["rain_ratio"],
            "shelf_basin_calcite_burial": row1["shelf_basin_calcite_burial"],
            "thc_flux": row1["thc_flux"],

            "climate_sensitivity": row2["climate_sensitivity"],
            "nsi": row2["nsi"],
            "ncc": row2["ncc"],
            "thc_weakening": row2["thc_weakening"],
            "deep_atlantic_carbon_injection": row2["deep_atlantic_carbon_injection"],
        })

        print(
            f"SUCCESS global_id={global_id}, "
            f"emi_total={emi_total}, ph_drop={ph_drop}, "
            f"RCO2={rco2}, delta_T={delta_t}"
        )

    except Exception as e:
        failed.append({
            "global_id": global_id,
            "task_id": task_id,
            "group1_index": group1_index,
            "group2_index": group2_index,
            "group1_id": group1_id,
            "sample_id": sample_id,
            "status": "failed",
            "error": repr(e),
            "traceback": traceback.format_exc(),

            "pCO2_init": row1["pCO2_init"],
            "silicate_weathering_init": row1["silicate_weathering_init"],
            "carbonate_weathering_init": row1["carbonate_weathering_init"],
            "rain_ratio": row1["rain_ratio"],
            "shelf_basin_calcite_burial": row1["shelf_basin_calcite_burial"],
            "thc_flux": row1["thc_flux"],

            "climate_sensitivity": row2["climate_sensitivity"],
            "nsi": row2["nsi"],
            "ncc": row2["ncc"],
            "thc_weakening": row2["thc_weakening"],
            "deep_atlantic_carbon_injection": row2["deep_atlantic_carbon_injection"],
        })

        print(f"FAILED global_id={global_id}")
        print(repr(e))
        traceback.print_exc()

    finally:
        os.chdir(BASE_DIR)
        shutil.rmtree(work_dir)


# =========================
# Save one file per batch
# =========================
result_file = output_dir / f"batch_{task_id:04d}_results.csv"
failed_file = output_dir / f"batch_{task_id:04d}_failed.csv"

if results:
    pd.DataFrame(results).to_csv(result_file, index=False)
    print(f"Saved successful results to {result_file}")

if failed:
    pd.DataFrame(failed).to_csv(failed_file, index=False)
    print(f"Saved failed results to {failed_file}")

print("Batch finished.")