# PVD Texture Coordinates Remapping Script

A Python utility script designed to remap texture coordinates for simulation results stored in `.pvd` and associated data files for ParaView visualization.

---

## Usage Example

```bash
python pvd_texture_coords.py --mode all chute_geo.pvd --aspect-ratio 2.44 --rotation 90 --scale 45.0 --frame-index 5300 --center 0.0 0.05 --flip-u --output-dir ./textured
```

---

## Notes & Performance Tips

- **Processing Time (`--mode all`)**: Running in `--mode all` copies and processes all `.pvd` / data files in the sequence, which can take some time depending on dataset size.
- **Previewing (`--mode single`)**: Use `--mode single` (e.g., targeting frame index `5300`) to quickly remap only a single frame. This allows you to inspect the results in ParaView and verify alignment.
