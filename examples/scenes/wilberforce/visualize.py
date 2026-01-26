import os
import sys
import numpy as np
import matplotlib

# Use non-interactive backend when no DISPLAY is available (headless test runs)
if "DISPLAY" not in os.environ:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt


if __name__ == "__main__":
    file = "../../../vtk_out/wilberforce_tracked.csv"
    usecols = (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    if not os.path.exists(file):
        print(f"File not found: {file}")
        sys.exit(1)

    data = np.loadtxt(file, delimiter=",", skiprows=1, usecols=usecols)
    print(f"data.shape: {data.shape}")
    t, x, y, z, ux, uy, uz, omegax, omegay, omegaz = data.T

    dt = t[1] - t[0]
    alpha = np.cumsum(omegax) * dt * 180 / np.pi
    beta = np.cumsum(omegay) * dt * 180 / np.pi
    gamma = np.cumsum(omegaz) * dt * 180 / np.pi

    # Clip time window for plotting
    start_time = 0.3
    end_time = 8.0
    mask = (t >= start_time) & (t <= end_time)
    if not np.any(mask):
        print(f"Warning: no data in [{start_time}, {end_time}]s — plotting full range")
    else:
        t = t[mask]
        alpha = alpha[mask]
        z = z[mask]

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    # Plot alpha in black on left y-axis
    ax.plot(t, alpha, color="black", label="alpha (deg)")
    ax.set_xlabel("time")
    ax.set_ylabel("alpha (deg)", color="black")
    ax.tick_params(axis='y', labelcolor='black')
    ax.set_xlim(t[0], t[-1])
    ax.grid(True)

    # Create a secondary y-axis for z (red)
    ax2 = ax.twinx()
    ax2.plot(t, z, color="red", label="z")
    ax2.set_ylabel("z", color="red")
    ax2.tick_params(axis='y', labelcolor='red')

    # Combined legend: gather handles from both axes
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    # ax.legend(handles1 + handles2, labels1 + labels2, loc="upper right")

    plt.tight_layout()

    # Always save a PNG (useful for headless runs); also show interactively if DISPLAY is set
    out = "visualize_wilberforce.png"
    fig.savefig(out, dpi=150)
    print(f"Saved figure to {out}")