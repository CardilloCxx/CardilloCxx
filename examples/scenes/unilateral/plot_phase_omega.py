#!/usr/bin/env python3
from pathlib import Path
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

p = Path(__file__).resolve().parents[3] / "vtk_out" / "unilateral_tracked.csv"
out_pdf = Path(__file__).resolve().parent / "unilateral_phase_wy_integral.pdf"
groups = defaultdict(lambda: {"t": [], "wy": []})
with p.open() as f:
    for r in csv.DictReader(f):
        name = r.get("name", "")
        groups[name]["t"].append(float(r["t"]))
        groups[name]["wy"].append(float(r["wy"]))

# Prefer explicit a-list, but fall back to any tracked unilateral rods.
a_tags = ["0p10", "0p20", "0p30"]
names = [f"unilateral_rod_a{a}" for a in a_tags if f"unilateral_rod_a{a}" in groups]
if not names:
    names = sorted([n for n in groups.keys() if n.startswith("unilateral_rod")])

if not names:
    raise SystemExit(f"No tracked unilateral rods found in {p}")

fig, axes = plt.subplots(
    nrows=1,
    ncols=len(names),
    figsize=(2 * len(names), 3.5),
    constrained_layout=True,
    sharey=True,
)
if len(names) == 1:
    axes = [axes]

for ax, name in zip(axes, names):
    t = groups[name]["t"]
    wy = groups[name]["wy"]
    if not t:
        continue

    # Ensure monotone time for integration.
    order = sorted(range(len(t)), key=lambda k: t[k])
    t_sorted = [t[k] for k in order]
    wy_sorted = [wy[k] for k in order]

    I = cumulative_trapezoid(wy_sorted, t_sorted, initial=0.0)
    I = [v - I[-1] for v in I]

    # Ensure 0 is centered on the x-axis by using symmetric limits.
    if I:
        max_abs = max(abs(min(I)), abs(max(I)))
        if max_abs > 0:
            ax.set_xlim(-max_abs, +max_abs)

    # Title from name, e.g. unilateral_rod_a0p10 -> a=0.10
    title = name
    if "_a" in name:
        a_part = name.split("_a", 1)[1]
        title = "a=" + a_part.replace("p", ".")

    ax.plot(I, wy_sorted, color="black", linewidth=1.25)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r"$\varphi$", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axes[0].set_ylabel(r"$\dot{\varphi}$", fontsize=12, rotation=0, labelpad=15)
fig.savefig(out_pdf)
