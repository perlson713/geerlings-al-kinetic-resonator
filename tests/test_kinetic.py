from __future__ import annotations

import unittest

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


if __name__ == "__main__":
    unittest.main()
