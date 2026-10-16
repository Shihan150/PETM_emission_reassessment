# PETM carbon-cycle and ocean-acidification simulations

Model simulation code, input data, results, and plotting tools for investigating carbon emissions, ocean pH, carbonate compensation depth (CCD), sediment carbonate preservation, and alkalinity feedbacks during the Paleocene–Eocene Thermal Maximum (PETM).

This package retains the current `code/iloscar_backend.py` implementation and the latest archived Monte Carlo ensemble (**v4**). Older ensembles, analysis versions, figure variants, short-circuit experiments, and exploratory notebooks have been removed. Each workflow contains its own identical copy of the current model backend and precomputed MyAMI coefficient table so it can be used independently.

## Contents

| Folder | Simulation code | Data and results |
| --- | --- | --- |
| `1_intermodel_comparison/` | `run_simulation.py`, `iloscar_backend.py` | Comparison workbook, proxy targets, initial states, emission forcing, selected scenario outputs, and the finalized comparison figure |
| `2_HPRC_Monte_Carlo/` | `spinup_states.py`, `run_second_stage_one_job.py`, `iloscar_backend.py`, Slurm scripts | Recovered v4 parameter designs, archived v4 batches, and v4 analysis |
| `3_Illustrative_examples/` | `run_simulation.py`, `iloscar_backend.py` | Conventional/modified feedback initial state, forcing, and paired output tables |

## Setup

Use Python 3.10 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Windows activation: `.venv\Scripts\activate`. Dependencies are not pinned because the source workspace does not provide a complete validated environment specification. Slurm is required only for cluster submission.

## 1. Intermodel comparison

The simulation driver contains three independent scenarios extracted from the current finite-sediment experiment notebook:

```bash
python 1_intermodel_comparison/run_simulation.py gutjahr2017
python 1_intermodel_comparison/run_simulation.py zeebe2009
python 1_intermodel_comparison/run_simulation.py this-study
```

- `gutjahr2017`: joint pH–carbon-isotope inversion using the supplied proxy target files and Gutjahr initial state.
- `zeebe2009`: forward simulation using the supplied Zeebe emission history.
- `this-study`: modified-weathering inversion using the current study's initial state and proxy targets.

The driver sets its working directory automatically. Model outputs are written to named experiment directories beside it. The original source scenario parameter values are retained; the Gutjahr/Zeebe scenarios use their configured default carbonate chemistry, while the current-study scenario uses MyAMI (`coeff_flag=2`). These are implementations of published scenarios in this model, not the original external model source packages.

To regenerate the comparison plot:

```bash
python 1_intermodel_comparison/plot_fig1_grl_style.py
```

The plot reads the processed workbook `data/fig1_fsed.xlsx` and writes a PDF/PNG beside the script. Archived selected outputs are in `model_outputs/`, and the finalized supplied figure is in `figures/`. Running simulations does not automatically replace the processed workbook; rebuilding the comparison from new simulations requires exporting the relevant series into the workbook sheets.

## 2. HPRC Monte Carlo — v4 only

### Analyze the latest archived ensemble

`second_stage_results_v4/` contains 12,000 successful result rows in 60 batch CSVs. Existing v4 analysis tables and figures are in `analysis/petm_mc_analysis_v4/`.

```bash
cd 2_HPRC_Monte_Carlo
python analyze_mc_results.py
```

The analysis reads v4 and updates `analysis/petm_mc_analysis_v4/`. Its `joint_ok` criterion combines warming of 4.5 < ΔT < 6.5 °C, Pacific minimum CCD between 2900 and 3300 m, and Atlantic CCD shoaling ≥2000 m. Depth-specific carbonate criteria are calculated separately. Consult the script for exact thresholds.

### Regenerate initial states and run simulations

`mc_accepted_runs.csv` and `mc_group2_samples.csv` were recovered from parameter columns in the archived v4 batches. Repeated parameter combinations were checked for consistency. The design comprises **40 Group I initial states × 300 Group II samples = 12,000 inversions**. Group II includes sampled pH decline and additional Atlantic stagnation; Group I includes phosphate scaling.

**The original saved v4 initial-state files are absent from the source workspace.** The nine saved states previously packaged belong to another design and have been removed. Before rerunning Stage 2, regenerate states from the recovered Group I parameters:

```bash
python spinup_states.py
# Or regenerate one state:
python spinup_states.py --run-id 2
```

This performs 20 Myr integrations using the current backend, skips existing state files, and writes `steady_state_fsed/mc_<run_id>_fsed_steady.dat`. The spin-up duration/settings are a reconstruction, not verified original v4 spin-up metadata. Regenerated states and resulting numerical outputs must be checked before treating a rerun as a reproduction of the archived ensemble.

For Slurm submission, activate the Python environment, edit cluster account/partition/resource settings in `main_job.sh`, then run:

```bash
bash submit.sh
```

The submission helper verifies all state files exist, creates `logs/`, and calculates the array size from the parameter tables: currently 240 tasks, each handling up to 50 inversions. This batching differs from the 60 archived files; scientific combinations are indexed from the recovered design.

A single batch can run locally:

```bash
SLURM_ARRAY_TASK_ID=0 python run_second_stage_one_job.py
```

The updated runner uses the current backend and MyAMI coefficients, applies phosphate scaling and Atlantic stagnation, scales the normalized target shape in `mc_inv_ph.csv` by each sample's pH decline, and records CCD/carbonate-depth diagnostics. Each inversion runs in its own temporary directory to prevent collisions between cluster tasks. New batch CSVs go to `generated_results/`, preserving the archived v4 data. Model failures are recorded in separate CSVs.

The analysis script currently targets the archived v4 directory; inspecting a new ensemble requires deliberately updating its `SOURCE` and output path. Numerical equivalence of the reconstructed runner to the unavailable original v4 execution environment has not been established.

## 3. Illustrative examples

Run the paired feedback scenarios:

```bash
python 3_Illustrative_examples/run_simulation.py conventional
python 3_Illustrative_examples/run_simulation.py modified
```

These use the supplied `petm_steady_fsed.dat` initial state and `emi_rw.dat` emission forcing. Scenario values were extracted from the current finite-sediment notebook. Outputs are written to `petm_test_no_rw/` and `petm_test_rw/`, replacing the supplied output tables if those directories already exist.

To plot the paired outputs and calculate feedback metrics:

```bash
python 3_Illustrative_examples/plot_alkalinity_feedback_figure3.py
```

Outputs go to `3_Illustrative_examples/figures/`.

## Data, provenance, and validation

CSV files hold samples, batch results, or time series; `.dat` files hold model states or forcing histories; `.npz` stores precomputed carbonate chemistry coefficients; `.xlsx` contains processed comparison/proxy data. Use individual headers and source code for units. Sediment carbonate fractions and wt% require conversion before comparison.

`FILE_MANIFEST.csv` lists the final distribution files and their provenance. Model equations are retained from the current workspace backend. Packaging changes include independent simulation drivers, v4 parameter-table recovery, path updates, a reconstructed state generator, the updated batch runner, and portable Slurm submission helpers.

The upstream MyAMI source checkout is not included. The supplied coefficient table supports the configured scenarios; regenerating the table requires upstream MyAMI software. Manuscripts, private notes, reference PDFs, environments, caches, and older versions are excluded. This package provides no additional license grant for third-party model code or data; preserve applicable attribution.

Python syntax, shell syntax, simulation CLI help, scenario input references, and v4 design consistency were checked. Full numerical integrations, plotting, and Slurm execution remain unvalidated because the available local environments lack required dependencies.

Upload the contents of this directory as the GitHub repository root. The adjacent `github_upload.zip` contains the same cleaned package.
