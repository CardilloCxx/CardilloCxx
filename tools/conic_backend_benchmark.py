#!/usr/bin/env python3
"""Compare the conic contact backends (Clarabel, QOCO, ConicXX, MOSEK, SCS) on real scenes.

For every scene this does two things:

1. End-to-end: runs the scene for --steps steps with each backend (solver.type=...) and the same
   tolerance (pj.tol_abs = pj.tol_rel = mosek.tol_* = --tol; SCS's eps defaults to pj.tol_*), and
   reports the time spent in the backend (sum of its Assembly/Setup/Update/Solve/Extract timers).
   Trajectories may diverge between backends, so these are "what a simulation costs" numbers.
2. Same-problem replay: the SCS run dumps every --dump-every-th cone program together with the warm
   start it built (scs.dump_dir); build/tests/conic_replay then solves exactly those problems with
   every backend and reports time and a common accuracy measure (see tests/conic_replay.cpp).

Requires a build with -DCARDILLO_WITH_SCS=ON -DCARDILLO_BUILD_TESTS=ON (and -DCARDILLO_WITH_MOSEK=ON
for MOSEK). Run from the repository root (scene assets are resolved relative to it):

    python3 tools/conic_backend_benchmark.py --out bench_out [--scenes chain jenga] [--steps 200]
"""

import argparse
import csv
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# executable -> base config (relative to examples/scenes)
SCENES = {
    "stacked_spheres": "stacked_spheres/scene.config",
    "chain": "chain/scene.config",
    "cardhouse": "cardhouse/scene.config",
    "jenga": "jenga/scene.config",
    "leaning_tower": "leaningTower/scene.config",
    "slinky": "slinky/scene_mosek.config",
    "domino": "domino/scene_mosek.config",
}

SOLVERS = {  # solver.type -> timer-name prefix in the timing breakdown
    "clarabel": "Clarabel",
    "qoco": "QOCO",
    "conicxx": "ConicXX",
    "mosek": "MOSEK",
    "scs": "SCS",
}

ROW = re.compile(r"^(?P<name>.+?)\s+(?P<calls>\d+)\s+(?P<excl>[0-9.eE+-]+)\s+(?P<pct>[0-9.eE+-]+)\s+(?P<avg>[0-9.eE+-]+)\s*$")


def read_dt(cfg: Path) -> float:
    for line in cfg.read_text().splitlines():
        m = re.match(r"\s*sim\.dt\s*=\s*([^#\s]+)", line)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"no sim.dt in {cfg}")


def write_config(base: Path, out: Path, solver: str, steps: int, tol: float, extra: dict) -> None:
    dt = read_dt(base)
    overrides = {
        "solver.type": solver,
        "sim.T": repr(steps * dt * (1 + 1e-9)),
        "output.interval_steps": "0",
        "debug.pj": "false",
        "pj.convergence_csv_dir": "",
        "pj.tol_abs": repr(tol),
        "pj.tol_rel": repr(tol),
        "mosek.tol_rel_gap": repr(tol),
        "mosek.tol_pfeas": repr(tol),
        "mosek.tol_dfeas": repr(tol),
        **extra,
    }
    text = base.read_text() + "\n# --- conic_backend_benchmark overrides ---\n"
    text += "".join(f"{k}={v}\n" for k, v in overrides.items())
    out.write_text(text)


def parse_breakdown(stdout: str, prefix: str):
    solver_s, total_s, solve_calls = 0.0, None, 0
    for line in stdout.splitlines():
        if line.startswith("Total (inclusive)"):
            total_s = float(line.split()[2])
            continue
        m = ROW.match(line)
        if m and m.group("name").startswith(prefix):
            solver_s += float(m.group("excl"))
            if m.group("name").startswith(prefix + " Solve"):
                solve_calls = int(m.group("calls"))
    return solver_s, total_s, solve_calls


def run_scene(exe: str, base: Path, args, out_dir: Path, solvers):
    rows = []
    dump_dir = out_dir / exe / "dumps"
    dump_dir.mkdir(parents=True, exist_ok=True)
    for solver in solvers:
        cfg = out_dir / exe / f"{solver}.config"
        extra = {}
        if solver == "scs":
            extra = {"scs.dump_dir": str(dump_dir), "scs.dump_every": str(args.dump_every), "scs.stats_csv": str(out_dir / exe / "scs_stats.csv")}
        write_config(base, cfg, solver, args.steps, args.tol, extra)
        t0 = time.time()
        p = subprocess.run([str(args.build / "examples" / exe), str(cfg)], cwd=ROOT, capture_output=True, text=True, timeout=args.timeout)
        wall = time.time() - t0
        (out_dir / exe / f"{solver}.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
        solver_s, total_s, calls = parse_breakdown(p.stdout, SOLVERS[solver])
        ok = p.returncode == 0
        rows.append({"scene": exe, "solver": solver, "ok": ok, "steps_solved": calls, "solver_s": solver_s, "total_s": total_s, "wall_s": wall,
                     "ms_per_step": 1e3 * solver_s / calls if calls else float("nan")})
        status = "ok" if ok else f"FAILED (exit {p.returncode}): " + (p.stderr.strip().splitlines() or ["?"])[-1]
        print(f"  {solver:9s} {status:10s} backend {solver_s:9.3f} s over {calls} steps", flush=True)
    return rows


def scs_stats(path: Path):
    if not path.exists():
        return {}
    r = list(csv.DictReader(path.open()))
    if not r:
        return {}
    g = lambda k: [float(x[k]) for x in r]
    return {
        "n": max(g("n")), "m": max(g("m")), "soc": max(g("n_soc")),
        "iters_med": statistics.median(g("iters")), "iters_max": max(g("iters")),
        "setup_ms_med": statistics.median(g("setup_ms")), "solve_ms_med": statistics.median(g("solve_ms")),
        "linsys_share": sum(g("lin_sys_ms")) / max(1e-12, sum(g("solve_ms"))),
        "scale_updates_mean": statistics.mean(g("scale_updates")),
        "aa_acc": sum(g("aa_accepted")), "aa_rej": sum(g("aa_rejected")),
        "inaccurate": sum(1 for x in r if x["status"] != "0"),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", type=Path, default=ROOT / "build")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--scenes", nargs="*", default=list(SCENES))
    ap.add_argument("--solvers", nargs="*", default=list(SOLVERS))
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--dump-every", type=int, default=20)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--ipm-tol", type=float, default=None, help="IPM tolerance in the replay (default: --tol)")
    ap.add_argument("--timeout", type=float, default=3600)
    args = ap.parse_args()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)

    e2e, replay_txt = [], {}
    for exe in args.scenes:
        base = ROOT / "examples" / "scenes" / SCENES[exe]
        (args.out / exe).mkdir(parents=True, exist_ok=True)
        print(f"== {exe} ({args.steps} steps)", flush=True)
        e2e += run_scene(exe, base, args, args.out, args.solvers)
        dumps = sorted((args.out / exe / "dumps").glob("step_*.bin"), key=lambda p: int(p.stem.split("_")[1]))
        if dumps and (args.build / "tests" / "conic_replay").exists():
            p = subprocess.run([str(args.build / "tests" / "conic_replay"), "--ipm-tol", repr(args.ipm_tol or args.tol), "--csv", str(args.out / exe / "replay.csv"), *map(str, dumps)],
                               capture_output=True, text=True, timeout=args.timeout)
            replay_txt[exe] = p.stdout + (p.stderr if p.returncode else "")
            print(replay_txt[exe], flush=True)

    with (args.out / "end_to_end.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(e2e[0]))
        w.writeheader()
        w.writerows(e2e)

    lines = [f"# Conic backend benchmark\n", f"steps per scene: {args.steps}, tolerance: {args.tol:g}, replay every {args.dump_every}th step\n",
             "## End to end (backend time per step, ms)\n", "| scene | " + " | ".join(args.solvers) + " |", "|---|" + "---|" * len(args.solvers)]
    for exe in args.scenes:
        cells = []
        for s in args.solvers:
            r = next((x for x in e2e if x["scene"] == exe and x["solver"] == s), None)
            cells.append("—" if r is None else (f"{r['ms_per_step']:.2f}" if r["ok"] else "failed"))
        lines.append(f"| {exe} | " + " | ".join(cells) + " |")
    lines.append("\n## SCS statistics (end-to-end run)\n")
    lines.append("| scene | n | m | SOC | iters med/max | setup ms | solve ms | linsys share | scale upd/step | AA acc/rej | inaccurate steps |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for exe in args.scenes:
        s = scs_stats(args.out / exe / "scs_stats.csv")
        if s:
            lines.append(f"| {exe} | {s['n']:.0f} | {s['m']:.0f} | {s['soc']:.0f} | {s['iters_med']:.0f}/{s['iters_max']:.0f} | {s['setup_ms_med']:.2f} | {s['solve_ms_med']:.2f} | "
                         f"{100 * s['linsys_share']:.0f}% | {s['scale_updates_mean']:.1f} | {s['aa_acc']:.0f}/{s['aa_rej']:.0f} | {s['inaccurate']} |")
    lines.append("\n## Same-problem replay\n")
    for exe, txt in replay_txt.items():
        lines += [f"### {exe}\n", "```", txt.rstrip(), "```\n"]
    (args.out / "summary.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {args.out / 'summary.md'}")


if __name__ == "__main__":
    sys.exit(main())
