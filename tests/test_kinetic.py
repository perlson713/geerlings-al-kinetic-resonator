from __future__ import annotations

import unittest
from math import hypot

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

        # Lock the Fig. 1-style hairpins: seven constant-radius semicircles per
        # side, mirrored about x=0. Parameter tuning may change straight-run
        # lengths or turn count, but must not distort these fold shapes.
        inductor = layout.centerlines.inductor
        self.assertAlmostEqual(inductor.width_um, 5.0)
        self.assertEqual(len(inductor.points), 184)
        left_points = inductor.points[:92]
        for row in range(7):
            radius = 12.5
            center_x = -25.0 if row % 2 == 0 else -181.5
            center_y = 222.5 + row * 25.0 + radius
            turn_start = 1 + row * 13
            for point in left_points[turn_start : turn_start + 13]:
                self.assertAlmostEqual(
                    hypot(point.x - center_x, point.y - center_y), radius, places=9
                )
        for left, mirrored_right in zip(left_points, reversed(inductor.points[92:])):
            self.assertAlmostEqual(left.x, -mirrored_right.x, places=9)
            self.assertAlmostEqual(left.y, mirrored_right.y, places=9)

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

    def test_approved_fig1_topology_idc_bottom_and_frequency(self) -> None:
        project = load_project("configs/design_fig1_square_9ghz_al200.toml")
        layout = generate_layout(project.resonator)
        cutout = next(poly for poly in layout.etch_polygons if poly.role == "ground_cutout")
        self.assertAlmostEqual(cutout.bounds.width, 424.0)
        self.assertAlmostEqual(cutout.bounds.height, 395.0)
        self.assertLess(abs(cutout.bounds.width / cutout.bounds.height - 1.0), 0.08)
        self.assertEqual(project.resonator.inductor_turns % 2, 0)
        self.assertIsNone(layout.centerlines.feedline)
        self.assertEqual(layout.ports, ())

        # An even turn count preserves the paper's full-width top return.
        inductor = layout.centerlines.inductor
        self.assertEqual(len(inductor.points), 158)
        left_top = inductor.points[78]
        right_top = inductor.points[79]
        self.assertAlmostEqual(left_top.x, -187.0)
        self.assertAlmostEqual(right_top.x, 187.0)
        self.assertAlmostEqual(left_top.y, 372.5)
        self.assertAlmostEqual(right_top.y, 372.5)

        # Fig. 1 IDC bottom: the lowest finger belongs to the left bus; the
        # first right-bus finger is one 15 um pitch above it. Both buses retain
        # their straight tails down to the common y=0 lower edge.
        fingers = layout.centerlines.capacitor_fingers
        left_bus, right_bus = layout.centerlines.capacitor_buses
        self.assertEqual(fingers[0].name, "capacitor_finger_00_left")
        self.assertEqual(fingers[1].name, "capacitor_finger_01_right")
        self.assertAlmostEqual(fingers[0].points[0].x, left_bus.points[0].x)
        self.assertAlmostEqual(fingers[0].points[0].y, 2.5)
        self.assertAlmostEqual(fingers[1].points[0].x, right_bus.points[0].x)
        self.assertAlmostEqual(fingers[1].points[0].y, 17.5)
        self.assertAlmostEqual(left_bus.points[0].y, 0.0)
        self.assertAlmostEqual(right_bus.points[0].y, 0.0)

        resonator = ReducedOrderResonator(
            pec_frequency_ghz=9.111777769251562,
            modal_impedance_ohm=100.0,
            branch_squares=569.3219026641857,
            meander_energy_fraction=0.98,
        )
        thin = evaluate_point(150.0, material=self.material, resonator=resonator)
        thick = evaluate_point(200.0, material=self.material, resonator=resonator)
        self.assertAlmostEqual(thin.frequency_ghz, 8.9888959880, places=7)
        self.assertAlmostEqual(thick.frequency_ghz, 9.0069087787, places=7)
        self.assertLess(abs(thin.frequency_ghz - 9.0), 0.012)
        self.assertLess(abs(thick.frequency_ghz - 9.0), 0.007)


if __name__ == "__main__":
    unittest.main()
