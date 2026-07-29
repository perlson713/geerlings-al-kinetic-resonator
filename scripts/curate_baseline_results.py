#!/usr/bin/env python3
"""Copy compact PEC baseline artifacts and normalize machine-local paths."""

from __future__ import annotations

import json
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "build" / "design_a"
NGSOLVE = SOURCE / "ngsolve"
OUTPUT = ROOT / "results" / "baseline_pec"


FILES = {
    SOURCE / "layout.svg": OUTPUT / "layout.svg",
    SOURCE / "layout-preview.png": OUTPUT / "layout-preview.png",
    NGSOLVE / "analysis_summary.json": OUTPUT / "analysis_summary.json",
    NGSOLVE / "convergence.csv": OUTPUT / "convergence.csv",
    NGSOLVE / "convergence.png": OUTPUT / "convergence.png",
    NGSOLVE / "finer" / "results.json": OUTPUT / "finer_results.json",
    NGSOLVE / "finer" / "modes.csv": OUTPUT / "finer_modes.csv",
    NGSOLVE / "finer" / "mesh_summary.json": OUTPUT / "finer_mesh_summary.json",
    NGSOLVE / "finer" / "mode_03_E_xy.png": OUTPUT / "mode_03_E_xy.png",
}


def _portable(value):
    if isinstance(value, dict):
        return {key: _portable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_portable(item) for item in value]
    if isinstance(value, str):
        normalized = value.replace("\\", "/")
        root = ROOT.as_posix()
        if normalized.casefold().startswith(root.casefold()):
            return "$PROJECT_ROOT" + normalized[len(root) :]
    return value


def main() -> int:
    missing = [str(path) for path in FILES if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing baseline artifact(s): " + ", ".join(missing))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for source, target in FILES.items():
        if source.suffix == ".json":
            data = json.loads(source.read_text(encoding="utf-8"))
            target.write_text(
                json.dumps(_portable(data), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
                newline="\n",
            )
        else:
            shutil.copy2(source, target)
    manifest = {
        "analysis_type": "full_wave_PEC_eigenmode_baseline",
        "source_root": "$PROJECT_ROOT/build/design_a",
        "published_files": sorted(path.name for path in FILES.values()),
        "raw_archive_release_tag": "em-results-v1",
    }
    (OUTPUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"curated {len(FILES)} files in {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
