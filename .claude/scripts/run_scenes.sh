#!/bin/bash
# Short smoke runs of example scenes with a reduced simulation time and without VTK output.
#
# Usage: .claude/scripts/run_scenes.sh OUT_DIR target:scene_dir:config:T [...]
#   e.g. .claude/scripts/run_scenes.sh /tmp/smoke chain:chain:scene.config:0.1 \
#            torsional_oscillator:torsional_oscillator:scene.config:1.0
# Runs from the repository root (scenes load meshes relative to the cwd), with build/ in the root.
# Prints one line per scene: exit code, last solver iteration count and the first error message.
# Timeout per scene: SCENE_TIMEOUT seconds (default 900).

set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OUT=$1
shift
mkdir -p "$OUT"
for sc in "$@"; do
    IFS=: read -r target dir cfg T <<< "$sc"
    bin="$ROOT/build/examples/$target"
    [ -x "$bin" ] || { echo "== $target: no binary $bin"; continue; }
    copy="$OUT/$target.config"
    cp "$ROOT/examples/scenes/$dir/$cfg" "$copy"
    printf "\nsim.T=%s\noutput.interval_steps=0\n" "$T" >> "$copy"
    (cd "$ROOT" && timeout "${SCENE_TIMEOUT:-900}" "$bin" "$copy" > "$OUT/$target.log" 2>&1; echo "exit $?" >> "$OUT/$target.log")
    iters=$(tr '\r' '\n' < "$OUT/$target.log" | grep -o 'iters=[0-9]*' | tail -1)
    err=$(grep -m1 -o 'what():.*\|unrecognized config key.*\|Divergence detected.*' "$OUT/$target.log")
    echo "== $target ($cfg, T=$T): $(tail -1 "$OUT/$target.log") $iters $err"
done
