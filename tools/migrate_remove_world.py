#!/usr/bin/env python3
"""
Safe codemod to remove local 'sys' aliases and replace 'sys.' with 'engine.'.
Backs up files to <file>.bak and reports files needing manual attention (World::, toWorld).
Run from repo root: python3 tools/migrate_remove_world.py
"""

import re
import sys
from pathlib import Path
from shutil import copy2

# Inventory: adjust or leave empty to auto-find
inventory = [
    "CardilloMPI/examples/scenes/double_pendulum/DoublePendulumScene.hpp",
    "CardilloMPI/examples/scenes/hangbride/HangbrideScene.hpp",
    "CardilloMPI/examples/scenes/constraintTest/ConstraintTestScene.hpp",
    "CardilloMPI/examples/scenes/spaghetti/SpaghettiScene.hpp",
    "CardilloMPI/examples/scenes/dzhanibekov/dzhanibekov.hpp",
    "CardilloMPI/examples/scenes/euler_disk/EulerDiskScene.hpp",
    "CardilloMPI/examples/scenes/gears/GearsScene.hpp",
    "CardilloMPI/examples/scenes/rail/RailScene.hpp",
    "CardilloMPI/examples/scenes/rodAssembly/RodAssemblyScene.hpp",
    "CardilloMPI/examples/scenes/metronome/MetronomeScene.hpp",
    "CardilloMPI/examples/scenes/domino/DominoScene.hpp",
    "CardilloMPI/examples/scenes/unilateral/UnilateralScene.hpp",
    "CardilloMPI/examples/scenes/net/NetScene.hpp",
    "CardilloMPI/examples/scenes/parcel/ParcelScene.hpp",
    "CardilloMPI/examples/scenes/chain/ChainScene.hpp",
    "CardilloMPI/examples/scenes/cantilever/CantileverScene.hpp",
    "CardilloMPI/examples/scenes/wilberforce/WilberforcePendulum.hpp",
    "CardilloMPI/examples/scenes/slinky/SlinkyScene.hpp",
    "CardilloMPI/examples/scenes/cardhouse/CardhouseScene.hpp",
    "CardilloMPI/examples/scenes/sphere_packing/scene.hpp",
    "CardilloMPI/examples/scenes/stacked_spheres/StackedSpheresScene.hpp",
    "CardilloMPI/examples/scenes/rotating_ball/RotatingBallScene.hpp",
    "CardilloMPI/examples/scenes/strandbeest/StrandbeestScene.hpp",
    "CardilloMPI/examples/scenes/heightmap/HeightmapScene.hpp",
    "CardilloMPI/examples/scenes/fabric/FabricScene.hpp",
    "CardilloMPI/examples/scenes/discreteRod/DiscreteRodScene.hpp",
    "CardilloMPI/examples/scenes/painleve/painleveScene.hpp",
    "CardilloMPI/examples/scenes/jenga/scene.hpp",
    "CardilloMPI/examples/scenes/woodpecker/WoodpeckerScene.hpp",
    "CardilloMPI/examples/scenes/leaningTower/LeaningTowerScene.hpp",
    "CardilloMPI/examples/scenes/softbody/SoftbodyTestScene.hpp",
]

# If inventory is empty, auto-discover under CardilloMPI/examples
if not inventory:
    repo_root = Path('.').resolve()
    inventory = [str(p.relative_to(repo_root)) for p in repo_root.glob('CardilloMPI/examples/**') if p.suffix in ('.hpp', '.h', '.cpp')]


def safe_read(path: Path):
    return path.read_text(encoding="utf-8")


def safe_write(path: Path, text: str):
    path.write_text(text, encoding="utf-8")


def process_file(path: Path):
    txt = safe_read(path)
    original = txt

    # Remove common 'sys' alias declarations
    txt, n1 = re.subn(r'^\s*auto\s*&\s*sys\s*=\s*engine\.world\(\)\s*;\s*$\n?', '', txt, flags=re.MULTILINE)
    txt, n2 = re.subn(r'^\s*auto\s*sys\s*=\s*engine\.world\(\)\s*;\s*$\n?', '', txt, flags=re.MULTILINE)
    txt, n3 = re.subn(r'^\s*auto\s*&\s*sys\s*=\s*engine\.system\(\)\s*;\s*$\n?', '', txt, flags=re.MULTILINE)
    txt, n4 = re.subn(r'^\s*auto\s*sys\s*=\s*engine\.system\(\)\s*;\s*$\n?', '', txt, flags=re.MULTILINE)

    # Replace sys. -> engine.
    txt, n5 = re.subn(r'\bsys\.', 'engine.', txt)

    changed = txt != original
    if changed:
        bak = path.with_suffix(path.suffix + ".bak")
        copy2(path, bak)
        safe_write(path, txt)

    # Flags for manual review
    needs_manual = False
    reasons = []
    if re.search(r'\bWorld::', txt):
        needs_manual = True
        reasons.append("World::")
    if re.search(r'\btoWorld\s*\(', txt):
        needs_manual = True
        reasons.append("toWorld(")
    if re.search(r'\bcardillo::World\b', txt):
        needs_manual = True
        reasons.append("cardillo::World")
    if re.search(r'\bengine\.world\(', txt):
        needs_manual = True
        reasons.append("engine.world(")

    return {
        "path": str(path),
        "modified": changed,
        "deleted_aliases": n1 + n2 + n3 + n4,
        "replacements": n5,
        "needs_manual": needs_manual,
        "reasons": reasons
    }


def main():
    repo_root = Path('.').resolve()
    results = []
    for p in inventory:
        path = repo_root / p
        if not path.exists():
            print(f"SKIP (missing): {p}")
            continue
        print(f"Processing: {p}")
        r = process_file(path)
        results.append(r)

    print("\nSummary")
    manual = [r for r in results if r["needs_manual"]]
    modified = [r for r in results if r["modified"]]
    print(f"Files modified: {len(modified)}")
    for r in modified:
        print(f" - {r['path']}: aliases_removed={r['deleted_aliases']}, replacements={r['replacements']}")

    if manual:
        print("\nFiles needing manual review (World:: / toWorld / engine.world remnants):")
        for r in manual:
            print(f" - {r['path']}: {', '.join(r['reasons'])}")

    else:
        print("\nNo manual-review flags detected.")


if __name__ == '__main__':
    main()
