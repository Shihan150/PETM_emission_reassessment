#!/usr/bin/env python3
"""Create GRL Figure 3 from the paired PETM alkalinity-feedback experiments."""

from pathlib import Path
import os
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
CODE = ROOT
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

RUNS = {
    "Conventional": CODE / "petm_test_no_rw",
    "Modified": CODE / "petm_test_rw",
}
COLORS = {"Conventional": "#6b6b6b", "Modified": "#2b8cbe"}
MM = 1 / 25.4


def load_run(path):
    summary = pd.read_csv(path / "summary.csv")
    # summary.csv contains duplicate Age headings; retain one explicit time column.
    summary = summary.iloc[:, 1:].copy()
    summary.rename(columns={summary.columns[0]: "time"}, inplace=True)
    surface = pd.read_csv(path / "Surface_dic_alk_d13c_ph.csv").rename(
        columns={"Age": "time"}
    )
    return summary, surface


def carbonate_ph_grid(dic_grid, alk_grid):
    """Evaluate MyAMI pH contours at representative PETM surface conditions."""
    old_cwd = Path.cwd()
    os.chdir(CODE)
    try:
        sys.path.insert(0, str(CODE))
        import iloscar_backend as model

        # Only the carbonate-system constants are needed here; define the model
        # geometry flags expected by load_const without initializing a new run.
        model.FTYS = 1
        model.FSED = 1
        model.FSEDCC = 1
        model.load_const()
        table = model.read_myami_ca_mg_table("myami_ca_mg_coeff_table.npz")
        ck0, ck1, ck2, ckb, ckw, _, _ = model.get_myami_coeffs_from_table(
            0.02, 0.03, table
        )
        shape = dic_grid.shape
        n = dic_grid.size
        _, _, _, _, ph, _, _ = model.fcsys(
            dic_grid.ravel() / model.RHO,
            alk_grid.ravel() / model.RHO,
            np.full(n, model.HGUESS),
            np.full(n, 28.0),
            np.full(n, 35.0),
            np.zeros(n),
            0.02,
            0.03,
            2,
            ck0,
            ck1,
            ck2,
            ckb,
            ckw,
        )
        return np.asarray(ph).reshape(shape)
    finally:
        os.chdir(old_cwd)


def interp(frame, column, time):
    return float(np.interp(time, frame["time"], frame[column]))


def cumulative_trapezoid(x, y):
    result = np.zeros_like(x, dtype=float)
    result[1:] = np.cumsum(np.diff(x) * (y[:-1] + y[1:]) / 2)
    return result


def run_metrics(summary, surface):
    p0 = float(summary.surface_pH.iloc[0])
    pac0 = float(summary.Pacific_CCD.iloc[0])
    atl0 = float(summary.Atlantic_CCD.iloc[0])
    wt0 = float(summary["wt%"].iloc[0])
    return {
        "initial_pH": p0,
        "minimum_pH": float(summary.surface_pH.min()),
        "pH_decline": p0 - float(summary.surface_pH.min()),
        "pacific_ccd_min": float(summary.Pacific_CCD.min()),
        "pacific_shoaling": pac0 - float(summary.Pacific_CCD.min()),
        "atlantic_ccd_min": float(summary.Atlantic_CCD.min()),
        "atlantic_shoaling": atl0 - float(summary.Atlantic_CCD.min()),
        "minimum_sediment_caco3": float(summary["wt%"].min() * 100),
        "sediment_caco3_drop": (wt0 - float(summary["wt%"].min())) * 100,
        "peak_pco2": float(summary.pco2.max()),
        "peak_warming": float(summary.Temp.max() - summary.Temp.iloc[0]),
        "dic_5k": interp(surface, "surface_dic", 5000),
        "alk_5k": interp(surface, "surface_alk", 5000),
        "pH_5k": interp(surface, "surface_pH", 5000),
    }


def main():
    data = {name: load_run(path) for name, path in RUNS.items()}
    metrics = {
        name: run_metrics(summary, surface)
        for name, (summary, surface) in data.items()
    }
    pd.DataFrame(metrics).T.to_csv(OUT / "Figure_3_alkalinity_feedback_metrics.csv")

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.0,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.3,
            "ytick.labelsize": 6.3,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
        }
    )

    fig = plt.figure(figsize=(180 * MM, 137 * MM))
    gs = fig.add_gridspec(2, 2)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 0])
    axd = fig.add_subplot(gs[1, 1])

    # (a) Alkalinity-flux anomalies in the modified experiment. Weathering
    # anomalies are relative to their initial fluxes; sediment compensation is
    # the reduction in the net sedimentary alkalinity sink relative to its
    # initial value.
    fluxes = {}
    for name, path in RUNS.items():
        flux = pd.read_csv(path / "flux.csv")
        flux.rename(columns={flux.columns[0]: "time"}, inplace=True)
        fluxes[name] = flux
    flux = fluxes["Modified"]
    t_flux = flux.time.to_numpy()
    ca_sil = flux.alk_sil.to_numpy() - float(flux.alk_sil.iloc[0])
    sediment = float(flux.alk_burial.iloc[0]) - flux.alk_burial.to_numpy()
    srw = data["Modified"][0]
    weather2 = 10.0 * (srw.pco2.to_numpy() / 1000.0) ** 0.2
    reverse = 10.0 * (srw.surface_pH.to_numpy() / srw.surface_pH.iloc[0]) ** 22.4
    added_summary = 2.0 * (weather2 - reverse)
    added = np.interp(t_flux, srw.time.to_numpy(), added_summary)
    modified_processes = [
        (ca_sil, "Ca-silicate weathering", ":"),
        (sediment, r"Sedimentary CaCO$_3$ compensation", "-"),
        (added, "Added alkalinity feedback", "-."),
    ]
    conventional_flux = fluxes["Conventional"]
    conventional_processes = [
        (conventional_flux.alk_sil.to_numpy() - float(conventional_flux.alk_sil.iloc[0]),
         ":"),
        (float(conventional_flux.alk_burial.iloc[0]) - conventional_flux.alk_burial.to_numpy(),
         "-"),
    ]
    for values, linestyle in conventional_processes:
        axa.plot(conventional_flux.time / 1000, values, lw=1.25, ls=linestyle,
                 color=COLORS["Conventional"], alpha=0.95)
    for values, label, linestyle in modified_processes:
        axa.plot(t_flux / 1000, values, lw=1.45, ls=linestyle,
                 color=COLORS["Modified"], label=label)
    axa.axhline(0, color="0.65", lw=0.55)
    axa.axvline(5, color="0.65", lw=0.55, ls=":")
    axa.set(xlabel="Time after onset (kyr)",
            ylabel=r"Alkalinity flux change (Tmol eq yr$^{-1}$)", xlim=(0, 50))
    axa.legend(handles=[
        Line2D([0], [0], color="0.25", lw=1.3, ls=":",
               label="Ca-silicate weathering"),
        Line2D([0], [0], color="0.25", lw=1.3, ls="-",
               label=r"Sedimentary CaCO$_3$ compensation"),
        Line2D([0], [0], color="0.25", lw=1.3, ls="-.",
               label="Added alkalinity feedback"),
    ], frameon=False, fontsize=5.0, loc="upper right", handlelength=2.2)

    # (b) DIC-alkalinity phase space with MyAMI pH contours.
    xv = np.linspace(2.16, 3.08, 150)
    yv = np.linspace(2.27, 3.12, 150)
    xx, yy = np.meshgrid(xv, yv)
    ph_grid = carbonate_ph_grid(xx, yy)
    levels_fill = np.arange(7.25, 7.86, 0.05)
    cf = axb.contourf(xx, yy, ph_grid, levels=levels_fill, cmap="YlGnBu", alpha=0.72)
    cs = axb.contour(xx, yy, ph_grid, levels=np.arange(7.3, 7.81, 0.1),
                     colors="0.35", linewidths=0.45)
    axb.clabel(cs, fmt="%.1f", fontsize=5.2, inline_spacing=2)
    for name, (_, surface) in data.items():
        axb.plot(surface.surface_dic, surface.surface_alk, color=COLORS[name], lw=1.5,
                 label=name, zorder=4)
        for t, marker in [(0, "o"), (5000, "s"), (50000, "^")]:
            axb.scatter(interp(surface, "surface_dic", t), interp(surface, "surface_alk", t),
                        marker=marker, s=22, color=COLORS[name], edgecolor="white",
                        linewidth=0.45, zorder=5)
    axb.set(xlabel=r"Mean surface DIC (mol m$^{-3}$)",
            ylabel=r"Mean surface alkalinity (mol m$^{-3}$)",
            xlim=(2.16, 3.08), ylim=(2.27, 3.12))
    cbar = fig.colorbar(cf, ax=axb, pad=0.02, fraction=0.047)
    cbar.set_label("pH at representative conditions", fontsize=6.2)
    cbar.ax.tick_params(labelsize=5.5)

    # (c) Atmospheric CO2 and global-mean warming.
    axt = axc.twinx()
    for name, (summary, _) in data.items():
        axc.plot(summary.time / 1000, summary.pco2, color=COLORS[name], lw=1.5)
        warming = summary.Temp - float(summary.Temp.iloc[0])
        axt.plot(summary.time / 1000, warming, color=COLORS[name], lw=1.25, ls="--")
    axc.set(xlabel="Time after onset (kyr)",
            ylabel=r"Atmospheric pCO$_2$ (ppm)", xlim=(0, 50))
    axt.set_ylabel(r"Global-mean warming ($^{\circ}$C)")
    axt.spines["top"].set_visible(False)
    axc.legend(handles=[
        Line2D([0], [0], color="0.25", lw=1.4, label=r"pCO$_2$"),
        Line2D([0], [0], color="0.25", lw=1.2, ls="--", label="Warming"),
    ], frameon=False, fontsize=5.8, loc="lower right")

    # (d) Basin CCD histories; depth axis increases downward.
    axd.axhspan(2900, 3300, color="#deebf7", alpha=0.65, lw=0)
    axd.axhline(4448 - 2000, color="#b2182b", ls=":", lw=0.9)
    for name, (summary, _) in data.items():
        axd.plot(summary.time / 1000, summary.Pacific_CCD, color=COLORS[name], lw=1.5)
        axd.plot(summary.time / 1000, summary.Atlantic_CCD, color=COLORS[name], lw=1.15,
                 ls="--")
    axd.set(xlabel="Time after onset (kyr)", ylabel="CCD depth (m)",
            xlim=(0, 50), ylim=(5100, 1900))
    axd.text(49, 3070, "Pacific constraint", ha="right", va="center", fontsize=5.7,
             color="#2166ac")
    axd.text(30, 2448, "Atlantic constraint", ha="left", va="bottom", fontsize=5.7,
             color="#b2182b")

    experiment_handles = [
        Line2D([0], [0], color=COLORS["Conventional"], lw=1.6,
               label="Baseline"),
        Line2D([0], [0], color=COLORS["Modified"], lw=1.6,
               label=r"+ non-CaCO$_3$ alkalinity feedback"),
    ]
    phase_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="0.35",
               markeredgecolor="white", markersize=5, label="0 kyr"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="0.35",
               markeredgecolor="white", markersize=5, label="5 kyr"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="0.35",
               markeredgecolor="white", markersize=5, label="50 kyr"),
    ]
    axb.legend(handles=phase_handles, frameon=False, fontsize=5.5,
               loc="upper left", ncol=3, columnspacing=0.8, handlelength=1.2)
    axd.legend(handles=[
        Line2D([0], [0], color="0.25", lw=1.4, label="Pacific"),
        Line2D([0], [0], color="0.25", lw=1.2, ls="--", label="Atlantic"),
    ], frameon=False, fontsize=5.8, loc="upper right", ncol=1)

    for ax, label in zip([axa, axb, axc, axd], "abcd"):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.15, 1.05, label, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=8.5, fontweight="bold", clip_on=False)

    fig.legend(handles=experiment_handles, frameon=False, fontsize=6.5,
               loc="upper center", bbox_to_anchor=(0.5, 1.005), ncol=2,
               columnspacing=2.0, handlelength=2.5)

    fig.subplots_adjust(left=0.095, right=0.91, bottom=0.10, top=0.93,
                        wspace=0.34, hspace=0.34)
    pdf = OUT / "Main_Figure_3_Alkalinity_feedback.pdf"
    png = OUT / "Main_Figure_3_Alkalinity_feedback.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    plt.close(fig)

    caption = (
        "Figure 3. Illustrative effect of a transient weathering-reverse-weathering "
        "imbalance on PETM carbonate chemistry. Both experiments receive the same "
        "9,500 Gt C input over 50 kyr (5,000 Gt C during the first 5 kyr and 4,500 "
        "Gt C during the subsequent 45 kyr). (a) Alkalinity-flux anomalies from "
        "Ca-silicate weathering and sedimentary carbonate "
        "compensation in both experiments, together with the added feedback in the "
        "modified experiment. Colors identify the experiments throughout the figure; "
        "line styles in (a) identify alkalinity processes. Weathering "
        "anomalies are relative to their initial fluxes; sedimentary compensation denotes "
        "the reduction in the net sedimentary alkalinity sink. "
        "(b) Area-weighted surface-ocean DIC and alkalinity trajectories. Background "
        "contours show pH calculated with MyAMI "
        "at representative surface conditions (28 degrees C, salinity 35, Ca = 20 "
        "mmol kg-1, and Mg = 30 mmol kg-1); symbols mark 0, 5, and 50 kyr. (c) "
        "Atmospheric pCO2 (solid) and global-mean warming (dashed). (d) Pacific "
        "(solid) and Atlantic "
        "(dashed) CCD histories. Blue shading marks the reconstructed Pacific event-CCD "
        "interval, and the dotted line marks 2,000 m of Atlantic shoaling. The experiments "
        "are mechanistic sensitivity tests and do "
        "not constitute a reconstruction of reverse-weathering fluxes."
    )
    (OUT / "Main_Figure_3_Alkalinity_feedback_caption.txt").write_text(caption + "\n")
    print(pdf)
    print(png)


if __name__ == "__main__":
    main()
