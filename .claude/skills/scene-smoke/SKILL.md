---
name: scene-smoke
description: Check a change against the example scenes with short runs (reduced sim.T, no VTK). Use after changes to the solvers, the assembler, constraints, collision or scene code, before handing a change or PR over for review.
---

# Scene smoke runs

1. Build first: `cmake --build build -j$(nproc)` (after changes to `config::Config` or other shared
   headers, all scene executables must be rebuilt) and run `cd build && ctest`.
2. Pick the scenes the change can affect and run them briefly:

   ```bash
   .claude/scripts/run_scenes.sh <out_dir> \
       torsional_oscillator:torsional_oscillator:scene.config:1.0 \
       pendulum:pendulum:scene.config:0.5 \
       chain:chain:scene.config:0.1 \
       cable:cable:scene.config:0.05 \
       wilberforce:wilberforce:scene_condensed.config:0.05
   ```

   Arguments are `target:scene_dir:config_file:sim_T`; targets are listed in
   `examples/CMakeLists.txt` (e.g. the `hangbride` directory builds `hangbridge`, `3DPrinter` builds
   `three_d_printer`). Use `<out_dir>` in the scratchpad/temp directory, not in the repository.
3. Keep runs short (≤ 0.1 s simulated for heavy scenes such as cable, net, 3DPrinter, conveyor_merge);
   long runs can take minutes to hours. Use `SCENE_TIMEOUT=<s>` to cap each run.
4. Interpret against `.claude/docs/known-issues.md`: a failure listed there is not a regression.
   For a suspected regression, run the same config on a build of `origin/main` before blaming the change.
5. For quantitative checks, compare the tracked CSV (`<output.folder>/<prefix>_tracked.csv`) of a run
   with VTK enabled against a reference run, e.g. the parcel of `conveyor_merge` after 0.6 s.
