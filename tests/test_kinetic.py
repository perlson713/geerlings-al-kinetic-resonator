from __future__ import annotations

import unittest

from geerlings_resonator.config import load_project
from geerlings_resonator.geometry import generate_layout

from geerlings_resonator.kinetic import (
    AluminiumLondonModel,
    ReducedOrderResonator,
    evaluate_point,
    inclusive_thicknesses,
)


class KineticInductanceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.material = AluminiumLondonModel()
        self.resonator = ReducedOrderResonator(
            pec_frequency_ghz=10.416985438861417,
            modal_impedance_ohm=100.0,
            branch_squares=618.9624168,
            meander_energy_fraction=0.98,
        )

    def test_published_al_fit_sheet_inductance(self) -> None:
        self.assertAlmostEqual(
            self.material.sheet_inductance_h(100.0) * 1.0e15,
            116.694047,
            places=5,
        )
        self.assertAlmostEqual(
            self.material.sheet_inductance_h(200.0) * 1.0e15,
            70.421856,
            places=5,
        )

    def test_frequency_increases_with_thickness_but_stays_below_pec(self) -> None:
        thin = evaluate_point(100.0, material=self.material, resonator=self.resonator)
        thick = evaluate_point(200.0, material=self.material, resonator=self.resonator)
        self.assertLess(thin.frequency_ghz, thick.frequency_ghz)
        self.assertLess(thick.frequency_ghz, self.resonator.pec_frequency_ghz)
        self.assertAlmostEqual(thin.frequency_ghz, 10.174466773, places=7)
        self.assertAlmostEqual(thick.frequency_ghz, 10.268589644, places=7)

    def test_inclusive_thickness_grid(self) -> None:
        self.assertEqual(
            inclusive_thicknesses(100.0, 200.0, 10.0),
            [float(value) for value in range(100, 201, 10)],
        )

    def test_invalid_values_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            self.material.sheet_inductance_h(0)
        with self.assertRaises(ValueError):
            inclusive_thicknesses(200, 100, 10)

    def test_square_9ghz_layout_and_two_allowed_thicknesses(self) -> None:
        project = load_project("configs/design_square_9ghz_al200.toml")
        layout = generate_layout(project.resonator)
        cutout = next(poly for poly in layout.etch_polygons if poly.role == "ground_cutout")
        self.assertAlmostEqual(cutout.bounds.width, 413.0)
        self.assertAlmostEqual(cutout.bounds.height, 420.0)
        self.assertLess(abs(cutout.bounds.width / cutout.bounds.height - 1.0), 0.02)
        self.assertIsNone(layout.centerlines.feedline)
        self.assertEqual(layout.ports, ())

        resonator = ReducedOrderResonator(
            pec_frequency_ghz=9.098137314852767,
            modal_impedance_ohm=100.0,
            branch_squares=569.5881997512568,
            meander_energy_fraction=0.98,
        )
        thin = evaluate_point(150.0, material=self.material, resonator=resonator)
        thick = evaluate_point(200.0, material=self.material, resonator=resonator)
        self.assertAlmostEqual(thin.frequency_ghz, 8.9755633219, places=7)
        self.assertAlmostEqual(thick.frequency_ghz, 8.9935313119, places=7)
        self.assertLess(abs(thin.frequency_ghz - 9.0), 0.03)
        self.assertLess(abs(thick.frequency_ghz - 9.0), 0.01)


if __name__ == "__main__":
    unittest.main()
