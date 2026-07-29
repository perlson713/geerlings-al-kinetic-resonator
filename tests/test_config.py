from __future__ import annotations

from pathlib import Path
import unittest

from geerlings_resonator.config import (
    apply_overrides,
    load_project,
    project_from_mapping,
)


ROOT = Path(__file__).resolve().parents[1]


class ProjectConfigTests(unittest.TestCase):
    def test_design_a_paper_values_and_counts(self) -> None:
        config = load_project(ROOT / "configs" / "design_a.toml")
        self.assertEqual(config.resonator.g_c_um, 10.0)
        self.assertEqual(config.resonator.g_l_um, 20.0)
        self.assertEqual(config.resonator.g_r_um, 10.0)
        self.assertEqual(config.resonator.capacitor_finger_count, 14)
        self.assertEqual(config.resonator.inductor_turns, 10)
        self.assertFalse(config.resonator.include_feedline)
        self.assertFalse(config.resonator.include_ports)
        # The Design-A working variant uses the user's nominal 175 nm Al film.
        self.assertEqual(config.stackup.nb_thickness_um, 0.175)

    def test_dotted_overrides_are_typed(self) -> None:
        config = load_project(
            ROOT / "configs" / "design_a.toml",
            overrides=["resonator.gC=40", "eigenmode.modes=10", "mesh.boundary=PEC"],
        )
        self.assertEqual(config.resonator.g_c_um, 40)
        self.assertEqual(config.eigenmode.modes, 10)
        self.assertEqual(config.mesh.boundary, "PEC")

    def test_override_requires_section_and_key(self) -> None:
        with self.assertRaises(ValueError):
            apply_overrides({}, ["bad=1"])

    def test_gds_and_stackup_metal_layers_must_match(self) -> None:
        with self.assertRaises(ValueError):
            project_from_mapping(
                {
                    "layers": {"metal_layer": 7},
                    "stackup": {"metal_layer": 8},
                }
            )


if __name__ == "__main__":
    unittest.main()
