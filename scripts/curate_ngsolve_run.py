#!/usr/bin/env python3
"""Copy compact NGSolve artifacts while replacing machine-local project paths."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]


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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    names = ("results.json", "mesh_summary.json", "modes.csv", "mode_03_E_xy.png")
    for name in names:
        origin = source / name
        target = output / name
        if name.endswith(".json"):
            data = json.loads(origin.read_text(encoding="utf-8"))
            target.write_text(
                json.dumps(_portable(data), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
                newline="\n",
            )
        else:
            shutil.copy2(origin, target)
    print(f"curated {len(names)} files from {source} to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
