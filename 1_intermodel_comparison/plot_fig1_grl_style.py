import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "data" / "fig1_fsed.xlsx"
OUTPUT_FILE = BASE_DIR / "fig1_scenario_comparison_0_50kyr_fsed_grl.pdf"

MM_TO_INCH = 1 / 25.4


def set_grl_style():
    """Apply compact AGU/GRL-oriented Matplotlib defaults."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "legend.fontsize": 5.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "lines.linewidth": 1.0,
            "lines.markersize": 3.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
        }
    )


def keep_0_50k(df):
    return df[(df["Age"] >= 0) & (df["Age"] <= 50000)].copy()


def load_data():
    gutjahr = pd.read_excel(DATA_FILE, sheet_name="gutjahr2017")
    gutjahr["Age"] = gutjahr["Age"] - 40999.5
    gutjahr["wt%"] = gutjahr["wt%"] / 100

    this_study = pd.read_excel(DATA_FILE, sheet_name="this_study")
    zeebe = pd.read_excel(DATA_FILE, sheet_name="zeebe2009")

    gutjahr_em = pd.read_excel(DATA_FILE, sheet_name="gutjahr_emission")
    this_study_em = pd.read_excel(DATA_FILE, sheet_name="this_study_emission")
    zeebe_em = pd.read_excel(DATA_FILE, sheet_name="zeebe2009_emission")
    target = pd.read_excel(DATA_FILE, sheet_name="Target")

    frames = [gutjahr, this_study, zeebe, gutjahr_em, this_study_em, zeebe_em, target]
    gutjahr, this_study, zeebe, gutjahr_em, this_study_em, zeebe_em, target = [
        keep_0_50k(df) for df in frames
    ]

    for df in [gutjahr, this_study, zeebe]:
        df["delta_Temp"] = df["Temp"] - df["Temp"].iloc[0]

    return gutjahr, this_study, zeebe, gutjahr_em, this_study_em, zeebe_em, target


def plot_line(ax, df, x, y, style, label=None, transform=None):
    values = df[y]
    if transform is not None:
        values = transform(values)
    ax.plot(df[x] / 1000, values, label=label, **style)


def main():
    set_grl_style()
    gutjahr, this_study, zeebe, gutjahr_em, this_study_em, zeebe_em, target = load_data()

    styles = {
        "gutjahr": {"ls": "-", "lw": 1.0, "color": "0.35"},
        "this_study": {"ls": "--", "lw": 1.0, "color": "0.05"},
        "zeebe": {"ls": ":", "lw": 1.2, "color": "0.05"},
    }

    labels = {
        "gutjahr": "cGENIE inversion (Gutjahr et al., 2017)",
        "this_study": "iLOSCAR inversion",
        "zeebe": "iLOSCAR forward (Zeebe et al., 2009)",
        "target": "Inversion target",
    }

    fig, axes = plt.subplots(
        4,
        2,
        figsize=(170 * MM_TO_INCH, 220 * MM_TO_INCH),
        sharex=True,
        constrained_layout=True,
    )
    axes = axes.ravel()

    plot_line(
        axes[0],
        gutjahr_em,
        "Age",
        "total_carbon_emission_Gt",
        styles["gutjahr"],
        labels["gutjahr"],
    )
    plot_line(
        axes[0],
        this_study_em,
        "Age",
        "total_carbon_emission_Gt",
        styles["this_study"],
        labels["this_study"],
    )
    plot_line(
        axes[0],
        zeebe_em,
        "Age",
        "total_carbon_emission_Gt",
        styles["zeebe"],
        labels["zeebe"],
    )
    axes[0].set_ylabel("Cumulative emission (Gt C)")

    plot_line(axes[1], this_study, "Age", "surface_d13c", styles["this_study"])
    plot_line(axes[1], zeebe, "Age", "surface_d13c", styles["zeebe"])
    axes[1].scatter(
        target["Age"] / 1000,
        target["d13c"],
        marker="o",
        s=10,
        facecolors="white",
        edgecolors="black",
        linewidths=0.6,
        label=labels["target"],
        zorder=5,
    )
    axes[1].set_ylabel(r"Mean surface $\delta^{13}$C (‰)")

    plot_line(
        axes[2],
        this_study,
        "Age",
        "surface_pH",
        styles["this_study"],
        transform=lambda x: x - x.iloc[0],
    )
    plot_line(
        axes[2],
        zeebe,
        "Age",
        "surface_pH",
        styles["zeebe"],
        transform=lambda x: x - x.iloc[0],
    )
    axes[2].scatter(
        target["Age"] / 1000,
        target["ph"] - target["ph"].iloc[0],
        marker="o",
        s=10,
        facecolors="white",
        edgecolors="black",
        linewidths=0.6,
        zorder=5,
    )
    axes[2].set_ylabel(r"$\Delta$ mean surface pH")

    plot_line(axes[3], gutjahr, "Age", "pco2", styles["gutjahr"])
    plot_line(axes[3], this_study, "Age", "pco2", styles["this_study"])
    plot_line(axes[3], zeebe, "Age", "pco2", styles["zeebe"])
    axes[3].set_ylabel(r"pCO$_2$ (ppm)")

    plot_line(axes[4], gutjahr, "Age", "delta_Temp", styles["gutjahr"])
    plot_line(axes[4], this_study, "Age", "delta_Temp", styles["this_study"])
    plot_line(axes[4], zeebe, "Age", "delta_Temp", styles["zeebe"])
    axes[4].set_ylabel(r"$\Delta$ temperature (°C)")

    plot_line(
        axes[5],
        gutjahr,
        "Age",
        "wt%",
        styles["gutjahr"],
        transform=lambda x: x * 100,
    )
    plot_line(
        axes[5],
        this_study,
        "Age",
        "wt%",
        styles["this_study"],
        transform=lambda x: x * 100,
    )
    plot_line(
        axes[5],
        zeebe,
        "Age",
        "wt%",
        styles["zeebe"],
        transform=lambda x: x * 100,
    )
    axes[5].set_ylabel(r"Surface sediment CaCO$_3$ (wt%)")

    plot_line(axes[6], this_study, "Age", "Pacific_CCD", styles["this_study"])
    plot_line(axes[6], zeebe, "Age", "Pacific_CCD", styles["zeebe"])
    axes[6].set_ylabel("Pacific CCD (m)")
    axes[6].invert_yaxis()

    plot_line(axes[7], this_study, "Age", "Atlantic_CCD", styles["this_study"])
    plot_line(axes[7], zeebe, "Age", "Atlantic_CCD", styles["zeebe"])
    axes[7].set_ylabel("Atlantic CCD (m)")
    axes[7].invert_yaxis()

    for ax, label in zip(axes, "abcdefgh"):
        ax.text(
            0.02,
            0.96,
            label,
            transform=ax.transAxes,
            fontsize=7,
            fontweight="bold",
            va="top",
            ha="left",
        )
        ax.set_xlim(0, 50)
        ax.set_xticks([0, 10, 20, 30, 40, 50])
        ax.tick_params(direction="out")
        ax.grid(False)
        ax.set_xlabel("")

    fig.supxlabel("Age from PETM onset (kyr)", fontsize=6)

    handles0, labels0 = axes[0].get_legend_handles_labels()
    handles1, labels1 = axes[1].get_legend_handles_labels()
    fig.legend(
        handles0 + handles1[-1:],
        labels0 + labels1[-1:],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.015),
        ncol=2,
        frameon=False,
        handlelength=2.8,
        columnspacing=1.0,
    )

    fig.savefig(OUTPUT_FILE)
    fig.savefig(OUTPUT_FILE.with_suffix(".png"), dpi=600)


if __name__ == "__main__":
    main()
