#!/bin/bash
# Sweeps per step and solver seconds per step from the timing breakdown of run_scene.sh logs.
# usage: summarize.sh NAME...   (RUN_DIR as for run_scene.sh)
R=${RUN_DIR:-$(pwd)/runs}
for name in "$@"; do
  f=$R/$name.log
  steps=$(grep -E "^Integration " "$f" | awk '{print $2/2}')
  sw=$(grep -E "Sweep " "$f" | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) {print $i; exit}}')
  tot=$(grep "Total (inclusive)" "$f" | awk '{print $3}')
  solver=$(grep -E "^(Sparse Schur|Naive Schur|Condensed|Projected Jacobi|Build & Factor)" "$f" | awk '{for(i=1;i<=NF;i++) if ($i ~ /^[0-9]+\.[0-9]+$/) {s+=$i; break}} END {print s}')
  echo "$name: steps=$steps sweeps/step=$(python3 -c "print(round($sw/$steps,1))") solver_s/step=$(python3 -c "print(round($solver/$steps,4))") total_s=$tot"
done
