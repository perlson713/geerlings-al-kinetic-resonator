"""Parse Palace eig.csv output without third-party dependencies."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import re
from typing import Iterable

from .errors import ExternalToolError


@dataclass(frozen=True, slots=True)
class EigenmodeResult:
    mode: int
    frequency_ghz: float
    imaginary_frequency_ghz: float
    quality_factor: float
    backward_error: float
    absolute_error: float

    def as_dict(self) -> dict[str, int | float]:
        return asdict(self)


def _normalise_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.strip().lower())


def _find_column(headers: list[str], *needles: str) -> str:
    normalised = {_normalise_header(header): header for header in headers}
    for key, original in normalised.items():
        if all(needle in key for needle in needles):
            return original
    raise ExternalToolError(
        f"eig.csv is missing a column matching {needles!r}; got {headers!r}"
    )


def read_eigenmodes(path: str | Path) -> list[EigenmodeResult]:
    source = Path(path)
    try:
        handle = source.open(newline="", encoding="utf-8-sig")
    except OSError as exc:
        raise ExternalToolError(f"Unable to open Palace results {source}: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        headers = reader.fieldnames or []
        mode_col = _find_column(headers, "m") if len(headers) == 1 else headers[0]
        real_col = _find_column(headers, "ref", "ghz")
        imag_col = _find_column(headers, "imf", "ghz")
        q_col = _find_column(headers, "q")
        backward_col = _find_column(headers, "error", "bkwd")
        absolute_col = _find_column(headers, "error", "abs")
        modes: list[EigenmodeResult] = []
        try:
            for row in reader:
                modes.append(
                    EigenmodeResult(
                        mode=int(float(row[mode_col])),
                        frequency_ghz=float(row[real_col]),
                        imaginary_frequency_ghz=float(row[imag_col]),
                        quality_factor=float(row[q_col]),
                        backward_error=float(row[backward_col]),
                        absolute_error=float(row[absolute_col]),
                    )
                )
        except (TypeError, ValueError, KeyError) as exc:
            raise ExternalToolError(f"Invalid row in Palace eig.csv: {exc}") from exc
    if not modes:
        raise ExternalToolError(f"No eigenmodes were found in {source}")
    return sorted(modes, key=lambda item: item.frequency_ghz)


def filter_frequency_window(
    modes: Iterable[EigenmodeResult],
    minimum_ghz: float | None = None,
    maximum_ghz: float | None = None,
) -> list[EigenmodeResult]:
    return [
        mode
        for mode in modes
        if (minimum_ghz is None or mode.frequency_ghz >= minimum_ghz)
        and (maximum_ghz is None or mode.frequency_ghz <= maximum_ghz)
    ]


def format_modes(modes: Iterable[EigenmodeResult]) -> str:
    rows = list(modes)
    headings = ("mode", "Re(f) GHz", "Im(f) GHz", "Q", "abs. error")
    values = [
        (
            str(row.mode),
            f"{row.frequency_ghz:.9g}",
            f"{row.imaginary_frequency_ghz:.3g}",
            f"{row.quality_factor:.6g}",
            f"{row.absolute_error:.3g}",
        )
        for row in rows
    ]
    widths = [
        max(len(headings[index]), *(len(row[index]) for row in values))
        for index in range(len(headings))
    ]
    line = "  ".join(headings[i].ljust(widths[i]) for i in range(len(headings)))
    separator = "  ".join("-" * width for width in widths)
    body = [
        "  ".join(row[i].rjust(widths[i]) for i in range(len(headings)))
        for row in values
    ]
    return "\n".join([line, separator, *body])


def write_results_json(path: str | Path, modes: Iterable[EigenmodeResult]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps([mode.as_dict() for mode in modes], indent=2) + "\n",
        encoding="utf-8",
    )
    return target
