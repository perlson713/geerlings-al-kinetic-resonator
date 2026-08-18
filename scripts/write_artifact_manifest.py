#!/usr/bin/env python3
"""Write SHA-256 hashes for tracked and pending repository files."""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "provenance" / "ARTIFACTS.sha256"


def _repository_files() -> list[Path]:
    raw = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT,
        stderr=subprocess.DEVNULL,
    )
    paths = [Path(item.decode("utf-8")) for item in raw.split(b"\0") if item]
    return sorted(path for path in paths if path.as_posix() != "provenance/ARTIFACTS.sha256")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    lines = [f"{_sha256(ROOT / path)}  {path.as_posix()}" for path in _repository_files()]
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(lines)} hashes to {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
