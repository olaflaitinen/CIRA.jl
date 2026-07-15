import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullLocator, NullFormatter
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.0,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#9400D3", "#E69F00", "#000000"]
TWO = ["#0072B2", "#D55E00"]
DISP = ["SIR", "SEIR", "Goodwin", "JAK-STAT", "PK", "MM", "LV", "HIV"]


def canon(name):
    n = name.lower()
    if "seir" in n:
        return "SEIR"
    if "goodwin" in n:
        return "Goodwin"
    if "jak" in n:
        return "JAK-STAT"
    if "pk" in n or "compartment" in n or "pharmac" in n:
        return "PK"
    if "michaelis" in n or n.strip() == "mm":
        return "MM"
    if "lotka" in n or n.strip() == "lv":
        return "LV"
    if "hiv" in n:
        return "HIV"
    return "SIR"


def read(name):
    with open(os.path.join(DATA, name)) as fh:
        rd = csv.reader(fh)
        header = next(rd)
        rows = [r for r in rd if r]
    return header, rows


def panel_label(ax, text):
    ax.text(-0.16, 1.04, text, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="right")


def save(fig, stem):
    fig.savefig(os.path.join(FIG, stem + ".png"), dpi=400)
    fig.savefig(os.path.join(FIG, stem + ".pdf"))
    plt.close(fig)


def figure1():
    _, rows = read("scaling.csv")
    ks = np.array([float(r[0]) for r in rows])
    ts = np.array([float(r[2]) for r in rows])
    slope, intercept = np.polyfit(np.log(ks), np.log(ts), 1)
    fig, ax = plt.subplots(figsize=(3.4, 2.7))
    xf = np.linspace(ks.min(), ks.max(), 100)
    ax.plot(xf, np.exp(intercept) * xf ** slope, "--", color="0.45", lw=0.9,
            label="Power-law fit (exponent %.2f)" % slope, zorder=1)
    ax.plot(ks, ts, "o-", color=PALETTE[0], lw=1.1, ms=3.5,
            label="CIRA reference pipeline", zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([2, 5, 10, 20, 50])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(1.8, 56)
    ax.set_xlabel("Number of parameters, q")
    ax.set_ylabel("Analysis runtime (s)")
    ax.legend(frameon=False, loc="upper left", handlelength=1.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "Fig1")


def figure2():
    _, prows = read("practical.csv")
    _, rprows = read("reparameterization.csv")
    _, nrows = read("robustness_noise.csv")
    fig = plt.figure(figsize=(7.2, 2.5))
    gs = GridSpec(1, 3, width_ratios=[1.0, 1.55, 1.0], wspace=0.42, figure=fig)

    ax = fig.add_subplot(gs[0, 0])
    sir = [r for r in prows if canon(r[0]) == "SIR"][0]
    vals = [float(sir[4]) * 100, float(sir[5]) * 100]
    ax.bar(["beta", "gamma"], vals, color=TWO, width=0.6, edgecolor="black", linewidth=0.4)
    ax.set_ylabel("Relative standard error (%)")
    ax.set_ylim(0, max(vals) * 1.25)
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals) * 0.03, "%.1f" % v, ha="center", va="bottom", fontsize=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    panel_label(ax, "a")

    ax2 = fig.add_subplot(gs[0, 1])
    nonid = [r for r in rprows if int(r[2]) < int(r[1])]
    labels = [canon(r[0]) for r in nonid]
    cf = np.array([float(r[4]) for r in nonid])
    cr = np.array([float(r[5]) for r in nonid])
    x = np.arange(len(labels))
    w = 0.4
    ax2.bar(x - w / 2, cf, width=w, color=PALETTE[1], edgecolor="black", linewidth=0.4, label="Full parameterization")
    ax2.bar(x + w / 2, cr, width=w, color=PALETTE[2], edgecolor="black", linewidth=0.4, label="CIRA reparameterization")
    ax2.set_yscale("log")
    ax2.set_ylim(1e0, 1e44)
    ax2.set_yticks([1e0, 1e8, 1e16, 1e24, 1e32, 1e40])
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=35, ha="right")
    ax2.set_ylabel("Fisher condition number")
    ax2.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=2, handlelength=1.4, columnspacing=1.2)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    panel_label(ax2, "b")

    ax3 = fig.add_subplot(gs[0, 2])
    xs = np.array([float(r[0]) * 100 for r in nrows])
    ax3.plot(xs, [float(r[1]) * 100 for r in nrows], "o-", color=TWO[0], lw=1.1, ms=3.5, label="beta")
    ax3.plot(xs, [float(r[2]) * 100 for r in nrows], "s-", color=TWO[1], lw=1.1, ms=3.5, label="gamma")
    ax3.set_xlabel("Measurement noise (%)")
    ax3.set_ylabel("Relative standard error (%)")
    ax3.set_ylim(0, None)
    ax3.legend(frameon=False, loc="upper left", handlelength=1.6)
    for s in ("top", "right"):
        ax3.spines[s].set_visible(False)
    panel_label(ax3, "c")

    save(fig, "Fig2")


def figure3():
    _, rrows = read("robustness_rank.csv")
    _, grows = read("robustness_gap.csv")
    factors = sorted(set(float(r[1]) for r in rrows))
    fac_labels = [("%.2f" % f).rstrip("0").rstrip(".") + "x" for f in factors]

    mat = np.full((len(DISP), len(factors)), np.nan)
    correct = np.zeros((len(DISP), len(factors)))
    for r in rrows:
        m = canon(r[0])
        i = DISP.index(m)
        j = factors.index(float(r[1]))
        mat[i, j] = int(r[3])
        correct[i, j] = 1.0 if int(r[3]) == int(r[4]) else 0.0

    fig = plt.figure(figsize=(7.2, 3.0))
    gs = GridSpec(1, 2, width_ratios=[1.1, 1.25], wspace=0.32, figure=fig)

    ax = fig.add_subplot(gs[0, 0])
    cmap = matplotlib.colors.ListedColormap(["#E8807F", "#7FB98A"])
    ax.imshow(correct, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(factors)))
    ax.set_xticklabels(fac_labels)
    ax.set_yticks(range(len(DISP)))
    ax.set_yticklabels(DISP)
    ax.set_xlabel("Integration grid density (relative to baseline)")
    for i in range(len(DISP)):
        for j in range(len(factors)):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, "%d" % int(mat[i, j]), ha="center", va="center", fontsize=6.5, color="black")
    ax.set_xticks(np.arange(-0.5, len(factors), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(DISP), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    panel_label(ax, "a")

    ax2 = fig.add_subplot(gs[0, 1])
    for i, name in enumerate(DISP):
        pts = [(float(r[1]), float(r[3])) for r in grows if canon(r[0]) == name]
        pts.sort()
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax2.plot(xs, ys, "o-", color=PALETTE[i], lw=1.0, ms=3.0, label=name)
    ax2.axhspan(1e-17, 1e-13, color="0.85", zorder=0)
    ax2.text(factors[-1], 3e-14, "numerical null space", ha="right", va="bottom", fontsize=5.5, color="0.35")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-17, 3e0)
    ax2.set_xlabel("Integration grid density (relative to baseline)")
    ax2.set_ylabel("Smallest identifiable singular value\n(normalized)")
    ax2.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), handlelength=1.4)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    panel_label(ax2, "b")

    save(fig, "Fig3")


def main():
    figure1()
    figure2()
    figure3()
    print("figures written to", FIG)
    print(sorted(os.listdir(FIG)))


if __name__ == "__main__":
    main()
