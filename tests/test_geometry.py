from __future__ import annotations

import unittest

from geerlings_resonator.geometry import generate_layout
from geerlings_resonator.parameters import DESIGN_A, DESIGN_B, DESIGN_C
from geerlings_resonator.svg import render_svg


class GeometryTests(unittest.TestCase):
    def test_design_a_fig1_topology_and_counts(self) -> None:
        layout = generate_layout(DESIGN_A)
        self.assertTrue(layout.topology.is_parallel_lc)
        self.assertEqual(len(layout.centerlines.capacitor_fingers), 14)
        horizontal = [
            (first, second)
            for first, second in zip(
                layout.centerlines.inductor.points,
                layout.centerlines.inductor.points[1:],
            )
            if first.y == second.y
        ]
        self.assertEqual(len(horizontal), 21)
        self.assertIsNone(layout.centerlines.feedline)
        self.assertEqual(layout.ports, ())

    def test_optional_feedline_and_ports(self) -> None:
        parameters = DESIGN_A.with_updates(
            include_feedline=True, include_ports=True
        )
        layout = generate_layout(parameters)
        self.assertEqual([port.layer for port in layout.ports], [201, 202])
        feedline = layout.metal_by_role()["feedline"][0].bounds
        coupling_ground = next(
            polygon.bounds
            for polygon in layout.metal_by_role()["ground"]
            if polygon.name == "ground_coupling_strip"
        )
        for port in layout.ports:
            self.assertEqual(port.polygon.bounds.y_min, coupling_ground.y_max)
            self.assertEqual(port.polygon.bounds.y_max, feedline.y_min)
        left_bridge = next(
            polygon.bounds
            for polygon in layout.metal_by_role()["ground"]
            if polygon.name == "ground_cpw_bridge_left"
        )
        right_bridge = next(
            polygon.bounds
            for polygon in layout.metal_by_role()["ground"]
            if polygon.name == "ground_cpw_bridge_right"
        )
        self.assertEqual(
            feedline.x_min - left_bridge.x_max,
            parameters.feedline_end_clearance_um,
        )
        self.assertEqual(
            right_bridge.x_min - feedline.x_max,
            parameters.feedline_end_clearance_um,
        )
        self.assertEqual(layout.get_port("left").polygon.bounds.x_min, feedline.x_min)
        self.assertEqual(layout.get_port("right").polygon.bounds.x_max, feedline.x_max)

    def test_paper_presets(self) -> None:
        self.assertEqual(
            (DESIGN_A.g_c_um, DESIGN_A.g_l_um, DESIGN_A.g_r_um, DESIGN_A.w_um, DESIGN_A.z0_ohm),
            (10, 20, 10, 5, 100),
        )
        self.assertEqual(
            (DESIGN_B.g_c_um, DESIGN_B.g_l_um, DESIGN_B.w_um, DESIGN_B.z0_ohm),
            (20, 5, 10, 200),
        )
        self.assertEqual(DESIGN_C.capacitor_trace_width_um, 40)
        self.assertEqual(DESIGN_C.inductor_trace_width_um, 10)
        self.assertEqual(DESIGN_C.z0_ohm, 300)

    def test_svg_is_deterministic_and_labels_assumptions(self) -> None:
        first = render_svg(DESIGN_A, show_etch_guides=True)
        second = render_svg(DESIGN_A, show_etch_guides=True)
        self.assertEqual(first, second)
        self.assertIn("14 IDC fingers", first)
        self.assertNotIn('data-role="feedline"', first)
        self.assertNotIn('data-layer="201"', first)

    def test_port_marker_layers_are_distinct(self) -> None:
        with self.assertRaises(ValueError):
            DESIGN_A.with_updates(
                include_feedline=True, include_ports=True, right_port_layer=201
            )
        with self.assertRaises(ValueError):
            DESIGN_A.with_updates(
                include_feedline=True, include_ports=True, left_port_layer=1
            )

    def test_ports_can_be_omitted_for_a_bare_layout(self) -> None:
        layout = generate_layout(
            DESIGN_A.with_updates(include_feedline=True, include_ports=False)
        )
        self.assertEqual(layout.ports, ())

    def test_feedline_can_be_removed_for_isolated_eigenmode(self) -> None:
        parameters = DESIGN_A.with_updates(
            include_feedline=False, include_ports=False
        )
        layout = generate_layout(parameters)
        self.assertIsNone(layout.centerlines.feedline)
        self.assertNotIn("feedline", layout.metal_by_role())
        self.assertEqual(
            [polygon.role for polygon in layout.etch_polygons], ["ground_cutout"]
        )
        self.assertEqual(layout.ports, ())
        with self.assertRaises(ValueError):
            DESIGN_A.with_updates(include_feedline=False, include_ports=True)


if __name__ == "__main__":
    unittest.main()
