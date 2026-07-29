from __future__ import annotations

import unittest
import xml.etree.ElementTree as ET

from geerlings_resonator.em import (
    EigenmodeParameters,
    MeshParameters,
    StackupParameters,
    ideal_lumped_targets,
    stackup_xml,
)
from geerlings_resonator.errors import ConfigurationError


class StackupTests(unittest.TestCase):
    def test_stackup_is_valid_gds2palace_xml(self) -> None:
        params = StackupParameters()
        root = ET.fromstring(stackup_xml(params))
        self.assertEqual(root.tag, "Stackup")
        layer = root.find("./ELayers/Layers/Layer")
        self.assertIsNotNone(layer)
        assert layer is not None
        self.assertEqual(layer.attrib["Name"], "Nb")
        self.assertEqual(layer.attrib["Layer"], "1")
        substrate = root.find("./ELayers/Layers/Substrate")
        self.assertEqual(substrate.attrib["Offset"], "430")  # type: ignore[union-attr]

    def test_invalid_stackup_is_rejected(self) -> None:
        with self.assertRaises(ConfigurationError):
            StackupParameters(nb_thickness_um=0).validate()

    def test_mesh_and_eigenmode_validation(self) -> None:
        MeshParameters().validate()
        EigenmodeParameters().validate()
        with self.assertRaises(ConfigurationError):
            MeshParameters(cells_per_wavelength=9).validate()
        with self.assertRaises(ConfigurationError):
            MeshParameters(boundary="PMC").validate()
        with self.assertRaises(ConfigurationError):
            EigenmodeParameters(target_frequency_ghz=12, search_max_frequency_ghz=8).validate()

    def test_design_a_six_ghz_lumped_targets(self) -> None:
        target = ideal_lumped_targets(6.0, 100.0)
        self.assertAlmostEqual(target.inductance_nh, 2.652582, places=5)
        self.assertAlmostEqual(target.capacitance_ff, 265.2582, places=3)


if __name__ == "__main__":
    unittest.main()
