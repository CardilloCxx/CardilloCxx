#!/usr/bin/env python3
import subprocess
import sys
import os
import csv
import time
import re
from pathlib import Path

# Configuration
RANKS = [1,2,4,6,8,10,12]
BINARY = str(Path(__file__).resolve().parents[1] / 'build' / 'bin' / 'collisions_demo')
BUILD_DIR = str(Path(__file__).resolve().parents[1] / 'build')
OUTPUT_DIR = str(Path(__file__).resolve().parent)
CSV_PATH = os.path.join(OUTPUT_DIR, 'benchmark_results.csv')
PNG_PATH = os.path.join(OUTPUT_DIR, 'benchmark_results.png')
# OpenMPI setting to avoid vader single copy issues seen previously
MPI_ENV = {
    'OMPI_MCA_btl_vader_single_copy_mechanism': 'none',
    'GMON_OUT_PREFIX': 'gmon.out'
}

# Patterns to parse gprof flat profile lines for functions of interest
# We'll aggregate by function name contains substring
SOLVER_PAT = re.compile(r'ProjectedJacobiSolver::iterateWithPreliminaryVelocity')
COMM_U_PAT = re.compile(r'Communication::exchangeBodyVelocitiesOwnerPush')
COMM_P_PAT = re.compile(r'Communication::exchangePercussionsOwnerPush')
REPLICATE_PAT = re.compile(r'Communication::replicateAllBodyVelocities')

# Helper to run a command and capture output

def run(cmd, cwd=None, env=None):
    e = os.environ.copy()
    if env:
        e.update(env)
    start = time.perf_counter()
    proc = subprocess.run(cmd, cwd=cwd, env=e, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    end = time.perf_counter()
    return proc.returncode, proc.stdout, proc.stderr, (end - start)

# Ensure build exists
if not os.path.isfile(BINARY):
    print(f"Binary not found at {BINARY}. Please build the project first.")
    sys.exit(1)

results = []

for ranks in RANKS:
    print(f"Running collisions_demo with {ranks} ranks...")
    # Clean previous gmon files
    for f in Path(BUILD_DIR).glob('gmon.out*'):
        try:
            f.unlink()
        except Exception:
            pass

    # Choose launcher and command
    if ranks == 1:
        cmd = [BINARY]
    else:
        cmd = ['mpirun', '-np', str(ranks), BINARY]

    # Measure total wall time for the run
    rc, out, err, wall = run(cmd, cwd=BUILD_DIR, env=MPI_ENV)
    if rc != 0:
        print(f"Run failed for ranks={ranks}: rc={rc}\nSTDERR:\n{err}")
        # still continue to next to collect what we can

    # gprof: if per-rank files exist (GMON_OUT_PREFIX), aggregate all, else use single
    print("Parsing gprof output...")
    gprof_cmd = ['gprof', BINARY]
    gmon_files = sorted([str(p) for p in Path(BUILD_DIR).glob('gmon.out*')])
    if gmon_files:
        gprof_cmd += gmon_files
    rc2, gout, gerr, _ = run(gprof_cmd, cwd=BUILD_DIR)
    if rc2 != 0:
        print(f"gprof failed for ranks={ranks}: rc={rc2}\nSTDERR:\n{gerr}")
        gout = ''

    # Aggregate times from flat profile
    total_solver = 0.0
    total_comm = 0.0

    for line in gout.splitlines():
        # Expected lines start with percent and seconds; we grab the self seconds column
        # Example: ' 61.23    454.52   454.52    60000 ... name'
        parts = line.strip().split()
        if len(parts) < 6:
            continue
        # Find function name at end of line
        name = line.strip()
        # Self seconds is the 3rd column
        try:
            self_sec = float(parts[2])
        except ValueError:
            continue
        if SOLVER_PAT.search(name):
            total_solver += self_sec
        elif COMM_U_PAT.search(name) or COMM_P_PAT.search(name) or REPLICATE_PAT.search(name):
            total_comm += self_sec

    results.append({
        'ranks': ranks,
        'wall_time_s': wall,
        'solver_time_s': total_solver,
        'communication_time_s': total_comm
    })

# Write CSV
with open(CSV_PATH, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['ranks','wall_time_s','solver_time_s','communication_time_s'])
    w.writeheader()
    for row in results:
        w.writerow(row)
print(f"Wrote CSV to {CSV_PATH}")

# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    ranks = [r['ranks'] for r in results]
    wall = [r['wall_time_s'] for r in results]
    solver = [r['solver_time_s'] for r in results]
    comm = [r['communication_time_s'] for r in results]

    plt.figure(figsize=(9,6))
    plt.plot(ranks, wall, '-o', label='Wall time (s)')
    plt.plot(ranks, solver, '-o', label='Solver time from gprof (s)')
    plt.plot(ranks, comm, '-o', label='Communication time from gprof (s)')
    plt.xlabel('Ranks')
    plt.ylabel('Time (seconds)')
    plt.title('collisions_demo runtime vs MPI ranks')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PNG_PATH, dpi=150)
    print(f"Wrote plot to {PNG_PATH}")
except Exception as e:
    print(f"Plotting failed: {e}\nYou can still use the CSV at {CSV_PATH}")
