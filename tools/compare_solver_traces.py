#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from scipy.ndimage import gaussian_filter1d
import numpy as np

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Case:
    scene: str
    label: str
    base_config: Path
    solver_type: str
    nesterov: bool
    overrides: dict[str, str]


TRACE_RE = re.compile(r"Simulating\s+\d+%\|.*?\|(\s*)(\d+)/(\d+)\s*\[[^\]]*?([\d.]+)\s*it/s\].*?iters=(\d+)")
TOTAL_WALLCLOCK_RE = re.compile(r"Total \(inclusive\)\s+([\d.]+)")


# Case.scene is used for display labels/output-dir grouping and doesn't always match the CMake
# target name registered in examples/CMakeLists.txt's cardillo_add_scene_example() calls.
SCENE_BINARY_OVERRIDES = {"3DPrinter": "three_d_printer"}


def resolve_bin_path(repo_root: Path, scene: str) -> Path:
    """Each scene now builds its own executable (commit 061ebe8 removed the single shared
    `build/bin/main` binary in favor of one target per scene via cardillo_add_scene_example()).
    There's no CMAKE_RUNTIME_OUTPUT_DIRECTORY override, so CMake's default applies: the binary
    lands under build/examples/ (the directory of the CMakeLists.txt that defines the target),
    not build/bin/."""
    target = SCENE_BINARY_OVERRIDES.get(scene, scene)
    p = repo_root / "build" / "examples" / target
    if not p.exists():
        raise FileNotFoundError(f"Missing per-scene binary: {p} (build target '{target}' first, e.g. `cmake --build build --target {target}`)")
    return p


def rewrite_config(base_config: Path, solver_type: str, nesterov: bool, overrides: dict[str, str], tmp_dir: Path) -> Path:
    text = base_config.read_text()

    def replace_or_append(key: str, value: str, src: str) -> str:
        pattern = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
        replacement = f"{key}={value}"
        if pattern.search(src):
            return pattern.sub(replacement, src)
        return src.rstrip() + "\n" + replacement + "\n"

    # Always enforce output.interval_steps=1000000 by default for all runs
    text = replace_or_append("output.interval_steps", "1000000", text)
    text = replace_or_append("solver.type", solver_type, text)
    text = replace_or_append("pj.nesterov", "true" if nesterov else "false", text)
    
    for key, value in overrides.items():
        text = replace_or_append(key, value, text)

    out_path = tmp_dir / f"{base_config.stem}_{solver_type}_{'nest' if nesterov else 'plain'}.config"
    out_path.write_text(text)
    return out_path


@dataclass(frozen=True)
class CaseResult:
    steps: list[int]
    iters: list[int]
    it_per_sec: list[float]
    total_wallclock_s: float | None
    log_path: Path


def run_case(repo_root: Path, case: Case, run_dir: Path, timeout_s: int) -> CaseResult:
    case_dir = run_dir / case.scene / case.label
    case_dir.mkdir(parents=True, exist_ok=True)

    bin_path = resolve_bin_path(repo_root, case.scene)
    tmp_config = rewrite_config(case.base_config, case.solver_type, case.nesterov, case.overrides, case_dir)
    log_path = case_dir / "console.log"

    try:
        proc = subprocess.run(
            [str(bin_path), "tools/" + str(tmp_config)],
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout_s,
        )
        output = proc.stdout
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode() if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        output = stdout + stderr

    log_path.write_text(output)

    steps: list[int] = []
    iters: list[int] = []
    it_per_sec: list[float] = []
    for line in output.splitlines():
        match = TRACE_RE.search(line)
        if match:
            steps.append(int(match.group(2)))
            it_per_sec.append(float(match.group(4)))
            iters.append(int(match.group(5)))

    wallclock_match = TOTAL_WALLCLOCK_RE.search(output)
    total_wallclock_s = float(wallclock_match.group(1)) if wallclock_match else None

    return CaseResult(steps, iters, it_per_sec, total_wallclock_s, log_path)


def worker(repo_root: Path, case: Case, log_root: Path, timeout_s: int) -> tuple[Case, CaseResult]:
    print(f"[START] Running {case.scene} / {case.label} ...")
    result = run_case(repo_root, case, log_root, timeout_s)
    print(f"[FINISH] {case.scene} / {case.label} captured {len(result.steps)} points, total_wallclock_s={result.total_wallclock_s}.")
    return case, result


import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.ndimage import gaussian_filter1d

def plot_combined_results(
    results: dict[tuple[str, str], tuple[list[int], list[int]]], 
    unique_scenes: list[str], 
    out_path: Path, 
    smoothness_fraction: float = 0.013
) -> None:
    colors = {
        ("pj", False): "#a0cbe8",
        ("pj", True): "#1f77b4",
        ("pgs", False): "#ff9896",
        ("pgs", True): "#d62728",
        ("clarabel", False): "#c5b0d5",
        ("clarabel", True): "#9467bd",
        ("condensed", False): "#98df8a",
        ("condensed", True): "#2ca02c",
    }
    styles = {False: "-", True: "--"}
    markers = {False: "o", True: "s"}

    num_scenes = len(unique_scenes)
    fig, axes = plt.subplots(num_scenes, 1, figsize=(13, 4.5 * num_scenes), sharex=False)

    if num_scenes == 1:
        axes = [axes]

    condensed_cycle = plt.get_cmap("tab10").colors

    def plot_series(ax: plt.Axes, steps: list[int], iters: list[int], label: str, color, linestyle: str, marker: str, max_steps: float) -> None:
        steps_arr = np.array(steps, dtype=float)
        iters_arr = np.array([it if it > 0 else 1 for it in iters], dtype=float)

        # Calculate sigma relative to the percentage of the dataset's total span
        sigma = 0.0
        if smoothness_fraction > 0.0 and len(steps_arr) > 1:
            avg_spacing = np.mean(np.diff(steps_arr))
            if avg_spacing > 0:
                target_window = max_steps * smoothness_fraction
                sigma = target_window / avg_spacing

        if sigma > 0.1:
            log_iters = np.log10(iters_arr)
            smoothed_log = gaussian_filter1d(log_iters, sigma=sigma, mode="nearest")
            line_iters = 10**smoothed_log
            ax.plot(steps_arr, iters_arr, linestyle="None", marker=marker, markersize=0.5, color=color, alpha=0.3)
        else:
            line_iters = iters_arr

        ax.plot(steps_arr, line_iters, linestyle=linestyle, marker=None, linewidth=2, color=color, label=label)
        ax.set_yscale("log")

    for ax, scene in zip(axes, unique_scenes, strict=True):

        max_steps = 0
        for (result_scene, label), (steps, _iters) in results.items():
            if result_scene == scene and steps:
                max_steps = max(max_steps, max(steps))

        for solver_type in ("pj", "pgs", "clarabel", "condensed"):
            for nesterov in (False, True):
                key = (scene, f"{solver_type}_{'nest' if nesterov else 'plain'}")
                if key not in results:
                    continue
                steps, iters = results[key]
                label = f"{solver_type.upper()} {'Nesterov' if nesterov else 'Plain'}"
                plot_series(ax, steps, iters, label, colors[(solver_type, nesterov)], styles[nesterov], markers[nesterov], max_steps)

        # Condensed cases don't fit the (solver_type, nesterov) grid -- one label per
        # sweep_mode x local_solve combination (see main()'s condensed case generation) -- so
        # plot every "condensed_*" result for this scene with its own color from a fixed cycle.
        condensed_labels = sorted(label for (result_scene, label) in results if result_scene == scene and label.startswith("condensed_"))
        for i, label in enumerate(condensed_labels):
            steps, iters = results[(scene, label)]
            plot_series(ax, steps, iters, label, condensed_cycle[i % len(condensed_cycle)], "-", "o", max_steps)

        display_title = scene.replace("_", " ").title() if scene != "3DPrinter" else "3DPrinter"

        ax.set_title(display_title, fontsize=12, fontweight="bold")
        ax.set_xlabel("Simulation step")
        ax.set_ylabel("Solver iterations")
        ax.grid(True, alpha=0.25)
        ax.legend(ncol=3, fontsize=9, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run solver comparisons and plot iteration traces into a combined plot.")
    parser.add_argument("--timeout", type=int, default=15, help="Seconds to allow each solver run before timing out.")
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_results/solver_compare"), help="Directory for logs and plots.")
    parser.add_argument("--workers", type=int, default=None, help="Number of concurrent threads. Defaults to CPU core count.")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        # Domino Scene Cases
        Case("domino", "pj_plain", repo_root / "examples/scenes/domino/scene_standard.config", "pj", False, {"sim.T": "2.0"}),
        Case("domino", "pj_nest", repo_root / "examples/scenes/domino/scene_standard.config", "pj", True, {"sim.T": "2.0"}),
        Case("domino", "pgs_plain", repo_root / "examples/scenes/domino/scene_standard.config", "pgs", False, {"sim.T": "2.0"}),
        Case("domino", "pgs_nest", repo_root / "examples/scenes/domino/scene_standard.config", "pgs", True, {"sim.T": "2.0"}),
        Case("domino", "clarabel_plain", repo_root / "examples/scenes/domino/scene_standard.config", "clarabel", False, {"sim.T": "2.0"}),
        
        # Slinky Scene Cases
        Case("slinky", "pj_plain", repo_root / "examples/scenes/slinky/scene_standard.config", "pj", False, {}),
        Case("slinky", "pj_nest", repo_root / "examples/scenes/slinky/scene_standard.config", "pj", True, {}),
        Case("slinky", "pgs_plain", repo_root / "examples/scenes/slinky/scene_standard.config", "pgs", False, {}),
        Case("slinky", "pgs_nest", repo_root / "examples/scenes/slinky/scene_standard.config", "pgs", True, {}),
        Case("slinky", "clarabel_plain", repo_root / "examples/scenes/slinky/scene_standard.config", "clarabel", False, {}),
        
        # 3DPrinter Scene Cases
        Case("3DPrinter", "pj_plain", repo_root / "examples/scenes/3DPrinter/scene.config", "pj", False, {}),
        Case("3DPrinter", "pj_nest", repo_root / "examples/scenes/3DPrinter/scene.config", "pj", True, {}),
        Case("3DPrinter", "pgs_plain", repo_root / "examples/scenes/3DPrinter/scene.config", "pgs", False, {}),
        Case("3DPrinter", "pgs_nest", repo_root / "examples/scenes/3DPrinter/scene.config", "pgs", True, {}),
        Case("3DPrinter", "clarabel_plain", repo_root / "examples/scenes/3DPrinter/scene.config", "clarabel", False, {}),

        # Stacked Spheres Scene Cases
        Case("stacked_spheres", "pj_plain", repo_root / "examples/scenes/stacked_spheres/scene.config", "pj", False, {}),
        Case("stacked_spheres", "pj_nest", repo_root / "examples/scenes/stacked_spheres/scene.config", "pj", True, {}),
        Case("stacked_spheres", "pgs_plain", repo_root / "examples/scenes/stacked_spheres/scene.config", "pgs", False, {}),
        Case("stacked_spheres", "pgs_nest", repo_root / "examples/scenes/stacked_spheres/scene.config", "pgs", True, {}),
        Case("stacked_spheres", "clarabel_plain", repo_root / "examples/scenes/stacked_spheres/scene.config", "clarabel", False, {}),
    ]

    # Condensed solver ablation: sweep_mode x local_solve, on the two scenes this solver was
    # specifically designed/validated against -- domino (pure rigid, contact-only) and slinky
    # (compliant chain + friction). Uses each scene's own scene_standard.config as the base (same
    # convention as the pj/pgs/clarabel cases above), overriding only the condensed.* keys; this
    # inherits scene_standard.config's already-tuned pj.alpha (domino's dense contact set needs
    # alpha=0.001 -- the generic default of 0.3 diverges in jacobi sweep mode on this scene).
    for scene, sim_t_override in (("domino", {"sim.T": "2.0"}), ("slinky", {})):
        base_config = repo_root / f"examples/scenes/{scene}/scene_standard.config"
        for sweep in ("gauss_seidel", "jacobi", "colored", "chaotic"):
            for local in ("projection", "newton"):
                cases.append(
                    Case(
                        scene,
                        f"condensed_{sweep}_{local}",
                        base_config,
                        "condensed",
                        False,
                        {**sim_t_override, "condensed.sweep_mode": sweep, "condensed.local_solve": local},
                    )
                )

    unique_scenes = sorted(list(set(case.scene for case in cases)))

    results: dict[tuple[str, str], tuple[list[int], list[int]]] = {}
    wallclock_rows: list[dict[str, str | float | None]] = []
    summary_rows: list[dict[str, str | int | float]] = []
    log_root = output_dir / "raw_logs"

    print(f"Starting parallel execution using up to {args.workers or 'default'} workers...")
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(worker, repo_root, case, log_root, args.timeout)
            for case in cases
        ]

        for future in futures:
            case, result = future.result()
            results[(case.scene, case.label)] = (result.steps, result.iters)
            wallclock_rows.append({"scene": case.scene, "case": case.label, "total_wallclock_s": result.total_wallclock_s})

            for step, it, itps in list(zip(result.steps, result.iters, result.it_per_sec, strict=True))[:10]:
                summary_rows.append({"scene": case.scene, "case": case.label, "step": step, "iters": it, "it_per_sec": itps})

    summary_csv = output_dir / "trace_summary.csv"
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["scene", "case", "step", "iters", "it_per_sec"])
        writer.writeheader()
        writer.writerows(summary_rows)

    wallclock_csv = output_dir / "wallclock_summary.csv"
    with wallclock_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["scene", "case", "total_wallclock_s"])
        writer.writeheader()
        writer.writerows(wallclock_rows)
    print(f"Saved wall-clock summary CSV to {wallclock_csv}")
    print(f"Saved summary CSV to {summary_csv}")

    plot_path = output_dir / "solver_comparison.png"
    plot_combined_results(results, unique_scenes, plot_path)
    print(f"Saved combined scene plots to {plot_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())