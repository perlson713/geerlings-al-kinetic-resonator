from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

from geerlings_resonator.em import EigenmodeParameters, MeshParameters
from geerlings_resonator.mesh import build_mesh_and_eigenmode_config


class _Ports:
    def __init__(self) -> None:
        self.portlayers: list[int] = []

    def add_port(self, port: SimpleNamespace) -> None:
        self.portlayers.append(port.source_layernum)


class _Metals:
    def getlayernumbers(self) -> list[int]:
        return [1]


class MeshBridgeTests(unittest.TestCase):
    def test_headless_bridge_and_generated_paths(self) -> None:
        captured: dict[str, object] = {}
        fake = ModuleType("gds2palace")
        fake.stackup_reader = SimpleNamespace(
            read_substrate=lambda _: (object(), object(), _Metals())
        )

        def read_gds(*args, **kwargs):  # type: ignore[no-untyped-def]
            captured["purposelist"] = kwargs["purposelist"]
            captured["layers"] = args[1]
            return object()

        fake.gds_reader = SimpleNamespace(read_gds=read_gds)

        def simulation_port(**kwargs):  # type: ignore[no-untyped-def]
            return SimpleNamespace(**kwargs)

        def create_palace(_, settings):  # type: ignore[no-untyped-def]
            captured["settings"] = settings
            target = Path(settings["sim_path"])
            mesh_name = f"{settings['model_basename']}.msh"
            (target / mesh_name).write_text("mesh", encoding="utf-8")
            config = {
                "Problem": {"Type": "Driven", "Output": "output/driven"},
                "Model": {"Mesh": mesh_name, "L0": 1.0e-6},
                "Domains": {"Materials": [{"Attributes": [1], "Permittivity": 1.0}]},
                "Boundaries": {
                    "Conductivity": [
                        {"Attributes": [11], "Conductivity": 1.0e10, "Thickness": 0.2}
                    ],
                    "Absorbing": {"Attributes": [20], "Order": 2},
                    "LumpedPort": [
                        {"Index": 1, "R": 50.0, "Attributes": [31], "Excitation": False}
                    ],
                },
                "Solver": {
                    "Order": 2,
                    "Driven": {"Samples": []},
                    "Linear": {"Type": "Default", "Tol": 1.0e-6, "MaxIts": 400},
                },
            }
            config_path = target / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            return str(config_path), "output/driven"

        fake.simulation_setup = SimpleNamespace(
            all_simulation_ports=_Ports,
            simulation_port=simulation_port,
            create_palace=create_palace,
        )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gds = root / "layout.gds"
            stackup = root / "stackup.xml"
            gds.write_bytes(b"gds")
            stackup.write_text("<Stackup/>", encoding="utf-8")
            with patch.dict(sys.modules, {"gds2palace": fake}):
                mesh_path, eigen_path = build_mesh_and_eigenmode_config(
                    gds,
                    stackup,
                    root / "palace",
                    MeshParameters(),
                    EigenmodeParameters(keep_passive_ports=True),
                    left_port_layer=301,
                    right_port_layer=302,
                    metal_datatype=5,
                    port_datatype=6,
                )
            eigen = json.loads(eigen_path.read_text(encoding="utf-8"))

        self.assertEqual(captured["purposelist"], [5, 6])
        self.assertEqual(captured["layers"], [1, 301, 302])
        settings = captured["settings"]
        assert isinstance(settings, dict)
        self.assertTrue(settings["nogui"])
        self.assertTrue(mesh_path.name.endswith(".msh"))
        self.assertEqual(eigen["Problem"]["Type"], "Eigenmode")
        self.assertEqual(eigen["Boundaries"]["Absorbing"]["Order"], 2)
        self.assertIn("LumpedPort", eigen["Boundaries"])


if __name__ == "__main__":
    unittest.main()
