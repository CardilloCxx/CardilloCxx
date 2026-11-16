#!/bin/bash
# Short scene run with config overrides and a final-state fingerprint.
# usage: run_scene.sh NAME TARGET SCENE_DIR BASE_CONFIG T [key=value ...]
#   ROOT    repository root with build/ (default: current directory)
#   RUN_DIR where configs and logs go (default: ./runs)
#   TMO     timeout in seconds (default: 1800)
name=$1; target=$2; dir=$3; cfg=$4; T=$5; shift 5
ROOT=${ROOT:-$(pwd)}
R=${RUN_DIR:-$ROOT/runs}
mkdir -p "$R"
cp "$ROOT/examples/scenes/$dir/$cfg" "$R/$name.cfg"
printf "\nsim.T=%s\noutput.interval_steps=0\n" "$T" >> "$R/$name.cfg"
for kv in "$@"; do echo "$kv" >> "$R/$name.cfg"; done
cd "$ROOT" || exit 1
s=$(date +%s.%N)
CARDILLO_DUMP_STATE=1 timeout "${TMO:-1800}" "build/examples/$target" "$R/$name.cfg" > "$R/$name.log" 2>&1
ec=$?
e=$(date +%s.%N)
wall=$(python3 -c "print(f'{$e-$s:.1f}')")
state=$(grep -o "totalKE=[^ ]* posNormSum=[^ ]* maxConstraintViolation=[^ ]*" "$R/$name.log")
echo "$name exit=$ec wall=${wall}s $state"
