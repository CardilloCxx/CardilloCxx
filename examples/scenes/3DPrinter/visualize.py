import os
import sys
import numpy as np
import matplotlib

# Use non-interactive backend when no DISPLAY is available
if "DISPLAY" not in os.environ:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt


if __name__ == "__main__":
    file = "../../../vtk_out/3DPrinter_tracked.csv"

    if not os.path.exists(file):
        print(f"File not found: {file}")
        sys.exit(1)

    # Robust loading (skip incomplete lines)
    data = np.genfromtxt(
        file,
        delimiter=",",
        names=True,
        dtype=None,
        encoding="utf-8",
        invalid_raise=False
    )

    # Remove rows with NaNs (from partially written lines)
    valid = (
        ~np.isnan(data["t"]) &
        ~np.isnan(data["px"]) &
        ~np.isnan(data["py"])
    )
    data = data[valid]

    # Filter gantry
    mask = data["name"] == "Gantry"
    if not np.any(mask):
        print("No entries found for 'gantry'")
        sys.exit(1)

    d = data[mask]

    x = d["px"]
    y = d["py"]

    # ---- 2D trajectory plot ----
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot(x, y, color="blue", linewidth=1.5)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (py)")
    ax.set_title("Gantry trajectory (x-y)")
    ax.axis("equal")
    ax.grid(True)

    # Mark start and end
    ax.scatter(x[0], y[0], color="green", label="start", zorder=3)
    ax.scatter(x[-1], y[-1], color="red", label="end", zorder=3)
    ax.legend()

    plt.tight_layout()

    out = "trajectory_xy.png"
    fig.savefig(out, dpi=150)
    print(f"Saved figure to {out}")

    if "DISPLAY" in os.environ:
        plt.show()