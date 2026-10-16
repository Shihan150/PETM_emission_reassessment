from pathlib import Path
import glob
import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
VERSION = "v4"
SOURCE = ROOT / f"second_stage_results_{VERSION}"
OUT = ROOT / f"analysis/petm_mc_analysis_{VERSION}"
OUT.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="ticks", context="paper")
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load_results():
    files = sorted(glob.glob(str(SOURCE / "batch_*_results.csv")))
    frames = [pd.read_csv(f).assign(source_batch=Path(f).name) for f in files]
    d = pd.concat(frames, ignore_index=True)
    d["weathering_ratio"] = d["carbonate_weathering_init"] / d["silicate_weathering_init"]
    d["warming_ok"] = d["delta_T"].gt(4.5) & d["delta_T"].lt(6.5)
    d["pacific_2900_ok"] = d["pacific_2900_min"].ge(0.20)
    d["pacific_depths_ok"] = d["pacific_2900_ok"] & d["pacific_3200_min"].between(0.01, 0.05)
    d["pacific_ccd_ok"] = d["pacific_ccd_min"].between(2900, 3300)
    d["atlantic_ok"] = d["atlantic_ccd_shoaling"].ge(2000)
    d["joint_ok"] = d["warming_ok"] & d["pacific_ccd_ok"] & d["atlantic_ok"]
    return d


d = load_results()
f = d[d["warming_ok"]].copy()

headline = {
    "source_files": int(d["source_batch"].nunique()),
    "total_rows": int(len(d)),
    "successful_model_runs": int(d["status"].eq("success").sum()),
    "initial_states": int(d["group1_id"].nunique()),
    "group2_designs": int(d["group2_index"].nunique()),
    "warming_filtered": int(len(f)),
    "pacific_ccd_success": int(f["pacific_ccd_ok"].sum()),
    "pacific_2900_success_strict_test": int(f["pacific_2900_ok"].sum()),
    "pacific_two_depth_success_strict_test": int(f["pacific_depths_ok"].sum()),
    "atlantic_success": int(f["atlantic_ok"].sum()),
    "joint_success": int(f["joint_ok"].sum()),
    "joint_fraction_filtered": float(f["joint_ok"].mean()),
    "joint_fraction_given_pacific": float(f.loc[f["pacific_ccd_ok"], "atlantic_ok"].mean()),
}
(OUT / "headline.json").write_text(json.dumps(headline, indent=2))

d.to_csv(OUT / f"combined_results_{VERSION}.csv", index=False)
f.to_csv(OUT / f"combined_results_{VERSION}_deltaT_4p5_to_6p5.csv", index=False)
f[f["joint_ok"]].to_csv(OUT / "joint_successes.csv", index=False)

group1 = [
    "pCO2_init", "silicate_weathering_init", "carbonate_weathering_init",
    "weathering_ratio", "rain_ratio", "shelf_basin_calcite_burial",
    "thc_flux", "po4_scaler",
]
group2 = [
    "climate_sensitivity", "nsi", "ncc", "thc_weakening",
    "atlantic_extra_stagnation", "deep_atlantic_carbon_injection",
    "ph_drop_target",
]
params = group1 + group2

# Distribution shifts. SMD is descriptive because rows share state/design values.
rows = []
for p in params:
    yes = f.loc[f["joint_ok"], p].dropna()
    no = f.loc[~f["joint_ok"], p].dropna()
    pooled = np.sqrt(((len(yes)-1)*yes.var() + (len(no)-1)*no.var()) / (len(yes)+len(no)-2))
    rho_p, p_p = spearmanr(f[p], f["pacific_ccd_shoaling"])
    rho_a, p_a = spearmanr(f[p], f["atlantic_ccd_shoaling"])
    rows.append({
        "parameter": p,
        "joint_n": len(yes),
        "joint_p05": yes.quantile(.05), "joint_p25": yes.quantile(.25),
        "joint_median": yes.median(), "joint_p75": yes.quantile(.75),
        "joint_p95": yes.quantile(.95), "filtered_median": f[p].median(),
        "standardized_mean_difference": (yes.mean()-no.mean())/pooled,
        "spearman_pacific_ccd": rho_p, "spearman_pacific_p": p_p,
        "spearman_atlantic_ccd": rho_a, "spearman_atlantic_p": p_a,
    })
profile = pd.DataFrame(rows)
profile.to_csv(OUT / "parameter_profiles.csv", index=False)

# Rates by initial state and repeated Group II design.
state_params = group1
state_summary = f.groupby(["group1_index", "group1_id"]).agg(
    eligible=("global_id", "size"),
    pacific_successes=("pacific_depths_ok", "sum"),
    joint_successes=("joint_ok", "sum"),
    **{p: (p, "first") for p in state_params},
).reset_index()
state_summary["pacific_rate"] = state_summary["pacific_successes"] / state_summary["eligible"]
state_summary["joint_rate"] = state_summary["joint_successes"] / state_summary["eligible"]
state_summary.to_csv(OUT / "state_summary.csv", index=False)

design_summary = f.groupby("group2_index").agg(
    eligible=("global_id", "size"),
    pacific_successes=("pacific_depths_ok", "sum"),
    joint_successes=("joint_ok", "sum"),
    **{p: (p, "first") for p in group2},
).reset_index()
design_summary["pacific_rate"] = design_summary["pacific_successes"] / design_summary["eligible"]
design_summary["joint_rate"] = design_summary["joint_successes"] / design_summary["eligible"]
design_summary.to_csv(OUT / "group2_summary.csv", index=False)

# Threshold sensitivity.
sens = []
for a in [1800, 1900, 2000, 2100, 2200, 2500]:
    for p29 in [0.15, 0.18, 0.20, 0.22, 0.25]:
        ok = f["pacific_2900_min"].ge(p29) & f["pacific_3200_min"].between(.01, .05) & f["atlantic_ccd_shoaling"].ge(a)
        sens.append({"atlantic_threshold_m": a, "pacific_2900_threshold": p29, "n": int(ok.sum()), "fraction": ok.mean()})
pd.DataFrame(sens).to_csv(OUT / "threshold_sensitivity.csv", index=False)

# Figure 1: interbasin CCD response.
fig, ax = plt.subplots(figsize=(6.5, 5.2))
sc = ax.scatter(f["pacific_ccd_shoaling"], f["atlantic_ccd_shoaling"],
                c=f["atlantic_extra_stagnation"], cmap="viridis", s=12, alpha=.45, linewidth=0)
j = f[f["joint_ok"]]
ax.scatter(j["pacific_ccd_shoaling"], j["atlantic_ccd_shoaling"],
           facecolors="none", edgecolors="#d62728", s=38, linewidth=1.0, label=f"Joint success (n={len(j)})")
ax.axhline(2000, color="black", ls="--", lw=1, label="Atlantic proxy threshold")
ax.set(xlabel="Pacific CCD shoaling (m)", ylabel="Atlantic CCD shoaling (m)")
cb = fig.colorbar(sc, ax=ax, pad=.02); cb.set_label("Additional Atlantic stagnation")
ax.legend(frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(OUT / "fig1_interbasin_ccd.png"); fig.savefig(OUT / "fig1_interbasin_ccd.pdf"); plt.close(fig)

# Figure 2: direct Pacific depth constraints.
fig, ax = plt.subplots(figsize=(6.5, 5.2))
sc = ax.scatter(f["pacific_2900_min"]*100, f["pacific_3200_min"]*100,
                c=f["atlantic_ccd_shoaling"], cmap="magma", s=13, alpha=.5, linewidth=0)
ax.axvline(20, color="black", ls="--", lw=1)
ax.axhspan(1, 5, xmin=(20-ax.get_xlim()[0])/(ax.get_xlim()[1]-ax.get_xlim()[0]), color="#2ca02c", alpha=.10)
ax.axhline(1, color="black", ls=":", lw=.8); ax.axhline(5, color="black", ls=":", lw=.8)
ax.scatter(j["pacific_2900_min"]*100, j["pacific_3200_min"]*100,
           facecolors="none", edgecolors="#00a6d6", s=38, linewidth=1)
ax.set(xlabel="Minimum CaCO$_3$ at 2900 m (wt%)", ylabel="Minimum CaCO$_3$ at 3200 m (wt%)")
cb = fig.colorbar(sc, ax=ax, pad=.02); cb.set_label("Atlantic CCD shoaling (m)")
fig.tight_layout(); fig.savefig(OUT / "fig2_pacific_depth_constraints.png"); fig.savefig(OUT / "fig2_pacific_depth_constraints.pdf"); plt.close(fig)

# Figure 3: parameter distribution shifts, standardized to sampled percentile.
plot_params = ["ph_drop_target", "climate_sensitivity", "thc_weakening",
               "atlantic_extra_stagnation", "deep_atlantic_carbon_injection",
               "weathering_ratio", "pCO2_init", "silicate_weathering_init"]
labels = ["pH decline", "Climate sensitivity", "THC weakening", "Atlantic stagnation",
          "Deep Atlantic injection", "Carbonate/silicate weathering", "Initial pCO$_2$", "Silicate weathering"]
long = []
for p, label in zip(plot_params, labels):
    lo, hi = f[p].min(), f[p].max()
    vals = (f.loc[f["joint_ok"], p]-lo)/(hi-lo)
    long.extend({"Parameter": label, "Sampled-range position": v} for v in vals)
long = pd.DataFrame(long)
fig, ax = plt.subplots(figsize=(7.2, 5.4))
sns.boxplot(data=long, y="Parameter", x="Sampled-range position", color="#7fcdbb", width=.6,
            showfliers=False, ax=ax)
sns.stripplot(data=long, y="Parameter", x="Sampled-range position", color="#253494", size=2.2,
              alpha=.35, jitter=.22, ax=ax)
ax.axvline(.5, color="0.35", ls="--", lw=.8)
ax.set_xlim(-.03, 1.03); ax.set_xlabel("Position within the sampled range (0 = minimum, 1 = maximum)")
ax.set_ylabel("")
fig.tight_layout(); fig.savefig(OUT / "fig3_joint_parameter_distributions.png"); fig.savefig(OUT / "fig3_joint_parameter_distributions.pdf"); plt.close(fig)

# Figure 4: binned joint success rates for the most relevant controls.
bin_params = ["ph_drop_target", "climate_sensitivity", "atlantic_extra_stagnation",
              "deep_atlantic_carbon_injection", "weathering_ratio", "pCO2_init"]
bin_labels = ["pH decline", "Climate sensitivity (°C)", "Atlantic stagnation",
              "Deep Atlantic injection", "Carbonate/silicate weathering", "Initial pCO$_2$ (ppm)"]
fig, axes = plt.subplots(2, 3, figsize=(10, 6.2))
for ax, p, label in zip(axes.flat, bin_params, bin_labels):
    tmp = f[[p, "joint_ok"]].copy()
    tmp["bin"] = pd.qcut(tmp[p], 8, duplicates="drop")
    g = tmp.groupby("bin", observed=True).agg(x=(p,"median"), rate=("joint_ok","mean"), n=("joint_ok","size"), k=("joint_ok","sum"))
    se = np.sqrt(g.rate*(1-g.rate)/g.n)
    ax.errorbar(g.x, 100*g.rate, yerr=100*1.96*se, marker="o", color="#2166ac", lw=1.3, capsize=2)
    ax.set(xlabel=label, ylabel="Joint success (%)")
fig.tight_layout(); fig.savefig(OUT / "fig4_binned_success_rates.png"); fig.savefig(OUT / "fig4_binned_success_rates.pdf"); plt.close(fig)

# Figure 5: initial-state dependence.
ss = state_summary.sort_values("joint_rate")
fig, ax = plt.subplots(figsize=(7.2, 6.2))
ax.barh(ss["group1_id"].astype(str), 100*ss["joint_rate"], color=np.where(ss.joint_successes.gt(0), "#3182bd", "#d9d9d9"))
ax.set(xlabel="Joint success among warming-filtered runs (%)", ylabel="Initial-state ID")
fig.tight_layout(); fig.savefig(OUT / "fig5_initial_state_success.png"); fig.savefig(OUT / "fig5_initial_state_success.pdf"); plt.close(fig)

print(json.dumps(headline, indent=2))
print("\nStrongest distribution shifts (absolute SMD):")
print(profile.reindex(profile.standardized_mean_difference.abs().sort_values(ascending=False).index)[
    ["parameter", "joint_median", "filtered_median", "standardized_mean_difference"]
].head(10).to_string(index=False))
