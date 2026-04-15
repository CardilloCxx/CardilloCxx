import os
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


if __name__ == "__main__":
    csv_file = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "../../../vtk_out/tube_tracked.csv")
    )

    if not os.path.exists(csv_file):
        raise SystemExit(f"Tracked CSV not found: {csv_file}")

    data = np.genfromtxt(csv_file, delimiter=",", names=True, dtype=None, encoding="utf-8")

    if data.size == 0:
        raise SystemExit("Tracked CSV is empty.")

    if data.ndim == 0:
        data = np.array([data], dtype=data.dtype)

    names = data["name"]
    if names.dtype.kind in ("U", "S", "O"):
        mask = names == "tube_sphere"
    else:
        mask = np.ones(len(data), dtype=bool)

    if np.any(mask):
        d = data[mask]
    else:
        d = data

    t = d["t"]
    x = d["px"]
    y = d["py"]
    z = d["pz"]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4))

    ax0.plot(t, x, label="x")
    ax0.plot(t, y, label="y")
    ax0.plot(t, z, label="z")
    ax0.set_xlabel("time [s]")
    ax0.set_ylabel("position [m]")
    ax0.set_title("Tube sphere position")
    ax0.grid(True)
    ax0.legend()

    ax1.plot(y, z, color="black")
    ax1.set_xlabel("y [m]")
    ax1.set_ylabel("z [m]")
    ax1.set_title("Cross-section path (y-z)")
    ax1.set_aspect("equal", adjustable="box")
    ax1.grid(True)

    plt.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "visualize_tube.png")
    fig.savefig(out, dpi=150)
    print(f"Saved figure to {out}")