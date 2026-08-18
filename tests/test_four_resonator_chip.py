from __future__ import annotations

import json
from pathlib import Path
import unittest

from geerlings_resonator.chip_layout import (
    merge_current_cavity_field_calibration,
    resonator_external_qc,
)


ROOT = Path(__file__).resolve().parents[1]


class FourResonatorChipTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        pin_calibration = json.loads(
            (ROOT / "configs" / "four_resonator_chip_fem_calibration.json").read_text(
                encoding="utf-8"
            )
        )
        field_result = json.loads(
            (
                ROOT
                / "results"
                / "four_resonator_chip"
                / "current_cavity_field_fem.json"
            ).read_text(encoding="utf-8")
        )
        cls.calibration = merge_current_cavity_field_calibration(
            pin_calibration, field_result
        )
        cls.field_result = field_result
        cls.summary = json.loads(
            (
                ROOT
                / "results"
                / "four_resonator_chip"
                / "optimization_summary.json"
            ).read_text(encoding="utf-8")
        )
        cls.square_summary = json.loads(
            (
                ROOT
                / "results"
                / "four_resonator_chip"
                / "centered_square_optimization_summary.json"
            ).read_text(encoding="utf-8")
        )

    def test_coupled_mode_formula(self) -> None:
        value = resonator_external_qc(9.0, 10.0, 1.0e3, 0.025)
        self.assertAlmostEqual(value, 1.44e6)
        with self.assertRaises(ValueError):
            resonator_external_qc(9.0, 10.0, 1.0e3, 0.0)

    def test_current_field_merge_has_complete_2d_grids(self) -> None:
        grids = self.calibration["field_grids"]
        self.assertEqual(set(grids), {"cavity_A", "cavity_B", "cavity_C"})
        for cavity in grids.values():
            self.assertEqual(set(cavity), {"left", "right"})
            for grid in cavity.values():
                self.assertEqual(len(grid["local_x_mm"]), 7)
                self.assertEqual(len(grid["y_mm"]), 7)
                self.assertEqual(
                    [len(row) for row in grid["e_in_plane_rms_normalized"]],
                    [7] * 7,
                )

    def test_rectangle_solution_is_uniform_and_clear(self) -> None:
        placements = self.summary["placements"]
        self.assertEqual(len(placements), 24)
        self.assertTrue(all(row["rotation_deg"] == 0.0 for row in placements))
        self.assertTrue(
            self.summary["constraints"]["centers_form_axis_aligned_rectangle"]
        )
        for cavity in ("cavity_A", "cavity_B", "cavity_C"):
            for chip in ("left", "right"):
                selected = [
                    row
                    for row in placements
                    if row["cavity"] == cavity and row["chip"] == chip
                ]
                self.assertEqual(len(selected), 4)
                xs = {round(float(row["chip_local_x_mm"]), 9) for row in selected}
                ys = {round(float(row["y_mm"]), 9) for row in selected}
                points = {
                    (
                        round(float(row["chip_local_x_mm"]), 9),
                        round(float(row["y_mm"]), 9),
                    )
                    for row in selected
                }
                self.assertEqual(len(xs), 2)
                self.assertEqual(len(ys), 2)
                self.assertEqual(points, {(x, y) for x in xs for y in ys})
                self.assertEqual(
                    {row["frequency_band"] for row in selected},
                    {"below_9ghz", "above_9ghz"},
                )
        ratios = self.summary["global_resonator_qc_max_to_min"]
        self.assertEqual(set(ratios), {"1e+06"})
        self.assertLess(ratios["1e+06"], 1.03)
        ratio_by_cavity = {
            row["cavity"]: row["resonator_qc_max_to_min"]
            for row in self.summary["pin_calibration"]
            if row["target_resonator_external_qc"] == 1.0e6
        }
        self.assertLess(ratio_by_cavity["cavity_A"], 1.005)
        self.assertLess(ratio_by_cavity["cavity_B"], 1.005)
        self.assertLess(ratio_by_cavity["cavity_C"], 1.027)
        self.assertGreaterEqual(
            self.summary["minimum_resonator_center_separation_mm"], 0.60 - 1.0e-9
        )
        self.assertGreaterEqual(
            self.summary["minimum_ground_cutout_edge_gap_mm"], 0.174 - 1.0e-6
        )
        tolerance = self.summary["position_tolerance_worst_case_qc_max_to_min"]
        self.assertLess(tolerance["independent_xy_error_envelope"]["10"], 1.031)

    def test_pin_lengths_match_fem_interpolation(self) -> None:
        nominal = [
            row
            for row in self.summary["pin_calibration"]
            if row["target_resonator_external_qc"] == 1.0e6
        ]
        by_case = {row["cavity"]: row for row in nominal}
        self.assertAlmostEqual(
            by_case["cavity_A"]["pin_length_mm"], 8.46, delta=0.02
        )
        self.assertAlmostEqual(
            by_case["cavity_B"]["pin_length_mm"], 8.93, delta=0.02
        )
        self.assertIn("cavity_C", by_case)
        self.assertGreaterEqual(by_case["cavity_C"]["pin_length_mm"], 7.75)
        self.assertLessEqual(by_case["cavity_C"]["pin_length_mm"], 11.0)

    def test_cavity_c_uses_15p2mm_current_field_and_explicit_pin_proxy(self) -> None:
        geometry = self.field_result["cases"]["cavity_C"]["geometry"]
        self.assertEqual(geometry["cavity_height_mm"], 15.2)
        self.assertEqual(
            self.calibration["provenance"]["pin_anchor_aliases_applied"],
            {"cavity_C": "cavity_B"},
        )

    def test_centered_square_is_exact_and_centered(self) -> None:
        self.assertEqual(
            self.square_summary["schema_id"],
            "geerlings.four_resonator_chip.centered_square_qc.v1",
        )
        self.assertTrue(self.square_summary["constraints"]["centered_on_chip"])
        self.assertEqual(self.square_summary["constraints"]["square_side_mm"], 0.8)
        placements = self.square_summary["placements"]
        for cavity in ("cavity_A", "cavity_B", "cavity_C"):
            for chip in ("left", "right"):
                selected = [
                    row
                    for row in placements
                    if row["cavity"] == cavity and row["chip"] == chip
                ]
                self.assertEqual(
                    {float(row["chip_local_x_mm"]) for row in selected},
                    {-0.4, 0.4},
                )
                self.assertEqual(
                    {float(row["y_mm"]) for row in selected}, {-0.4, 0.4}
                )
        self.assertAlmostEqual(
            self.square_summary["minimum_resonator_center_separation_mm"], 0.8
        )
        self.assertLess(
            self.square_summary["global_resonator_qc_max_to_min"]["1e+06"],
            1.17,
        )

    def test_committed_gds_readback_is_complete(self) -> None:
        result = ROOT / "results" / "four_resonator_chip"
        verification = json.loads(
            (result / "gds_readback_verification.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(verification), 6)
        for name, row in verification.items():
            self.assertIn("below9_", name)
            self.assertIn("above9_", name)
            self.assertEqual(len(row["below_9ghz_patterns"]), 2)
            self.assertEqual(len(row["above_9ghz_patterns"]), 2)
            self.assertGreater((result / name).stat().st_size, 100_000)
            self.assertEqual(row["top_cells"], 1)
            self.assertEqual(row["bounding_box_um"], [5050.0, 5050.0])
            self.assertGreater(row["shape_counts"]["1/0"], 1000)
            self.assertEqual(row["shape_counts"]["100/0"], 1)
            self.assertEqual(row["shape_counts"]["101/0"], 4)

    def test_centered_square_gds_readback_is_complete(self) -> None:
        result = ROOT / "results" / "four_resonator_chip"
        verification = json.loads(
            (result / "centered_square_gds_readback_verification.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(len(verification), 6)
        for name, row in verification.items():
            self.assertIn("centered_square", name)
            self.assertTrue((result / name).is_file())
            self.assertEqual(row["top_cells"], 1)
            self.assertEqual(row["bounding_box_um"], [5050.0, 5050.0])
            self.assertEqual(row["shape_counts"]["101/0"], 4)

    def test_fullwave_scope_is_not_overstated(self) -> None:
        provenance = self.summary["fem_calibration_provenance"]
        self.assertTrue(provenance["current_depth_field_revalidated"])
        self.assertFalse(provenance["current_geometry_port_sweep_revalidated"])
        self.assertFalse(
            provenance["current_field_fem"]["validity"][
                "sma_bore_and_pin_included"
            ]
        )
        self.assertIn("coupled-mode", self.summary["validity"]["resonator_external_qc"])


if __name__ == "__main__":
    unittest.main()
