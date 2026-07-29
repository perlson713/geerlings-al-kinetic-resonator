from __future__ import annotations

import ast
import unittest

from geerlings_resonator.em import EigenmodeParameters, MeshParameters
from geerlings_resonator.setupem import setupem_argv, setupem_model_source


class SetupEmTests(unittest.TestCase):
    def test_command_line(self) -> None:
        self.assertEqual(
            setupem_argv("layout.gds", "stackup.xml"),
            ["setupEM", "-gdsfile", "layout.gds", "-xmlfile", "stackup.xml"],
        )

    def test_generated_model_is_python_and_has_ports(self) -> None:
        source = setupem_model_source(
            gds_filename="layout.gds",
            stackup_filename="stackup.xml",
            mesh=MeshParameters(),
            eigenmode=EigenmodeParameters(),
            left_port_layer=301,
            right_port_layer=302,
            metal_datatype=5,
            port_datatype=6,
        )
        ast.parse(source)
        self.assertIn("source_layernum=301", source)
        self.assertIn("source_layernum=302", source)
        self.assertIn("purposelist=[5, 6]", source)
        self.assertIn('target_layername="Nb"', source)
        self.assertIn("margin = 150.0", source)
        port_lines = [
            line for line in source.splitlines() if "simulation_port(" in line
        ]
        self.assertEqual(len(port_lines), 2)
        for line in port_lines:
            arguments = line.split("simulation_port(", 1)[1].rstrip(")")
            ast.parse(f"f({arguments})", mode="eval")


if __name__ == "__main__":
    unittest.main()
