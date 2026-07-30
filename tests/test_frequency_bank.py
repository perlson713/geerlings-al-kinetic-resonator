from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import unittest

from geerlings_resonator.config import load_project
from geerlings_resonator.geometry import generate_layout


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "generate_frequency_bank", ROOT / "scripts" / "generate_frequency_bank.py"
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FrequencyBankTests(unittest.TestCase):
    def test_even_grid_is_centered_and_ten_mhz_apart(self) -> None:
        targets = MODULE.symmetric_targets(9.0, 10.0, 8)
        self.assertEqual(len(targets), 8)
        self.assertAlmostEqual(sum(targets) / len(targets), 9.0)
        self.assertEqual([round(value, 3) for value in targets], [
            8.965,
            8.975,
            8.985,
            8.995,
            9.005,
            9.015,
            9.025,
            9.035,
        ])
        for lower, upper in zip(targets, targets[1:]):
            self.assertAlmostEqual((upper - lower) * 1.0e3, 10.0)

    def test_invalid_grid_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.symmetric_targets(9.0, 10.0, 0)
        with self.assertRaises(ValueError):
            MODULE.symmetric_targets(9.0, 0.0, 8)

    def test_generated_bank_preserves_topology_and_spacing(self) -> None:
        result_root = ROOT / "results" / "frequency_bank_8x"
        with (result_root / "frequency_table.csv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        summary = json.loads((result_root / "results.json").read_text(encoding="utf-8"))
        self.assertEqual(len(rows), 8)
        self.assertEqual(summary["count"], 8)
        widths = [float(row["nominal_width_um"]) for row in rows]
        predicted = [float(row["predicted_al200_ghz"]) for row in rows]
        self.assertTrue(all(left > right for left, right in zip(widths, widths[1:])))
        for lower, upper in zip(predicted, predicted[1:]):
            self.assertAlmostEqual((upper - lower) * 1.0e3, 10.0, delta=0.05)
        for row in rows:
            self.assertLess(abs(float(row["target_error_al200_mhz"])), 0.03)
            pattern_id = row["pattern_id"]
            pattern_root = result_root / "patterns" / pattern_id
            manifest = json.loads((pattern_root / "manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["topology"]["fig1_idc_bottom"])
            self.assertTrue(manifest["topology"]["full_width_top_return"])
            self.assertTrue(manifest["topology"]["constant_radius_folds"])
            self.assertGreater((pattern_root / "layout.gds").stat().st_size, 1024)
            project = load_project(
                ROOT / "configs" / "design_fig1_square_9ghz_al200.toml",
                overrides=[f"resonator.nominal_width_um={row['nominal_width_um']}"],
            )
            topology = MODULE._topology_guard(generate_layout(project.resonator))
            self.assertTrue(topology["fig1_idc_bottom"])


if __name__ == "__main__":
    unittest.main()
