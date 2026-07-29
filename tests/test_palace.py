from __future__ import annotations

import unittest

from geerlings_resonator.em import EigenmodeParameters
from geerlings_resonator.errors import ConfigurationError
from geerlings_resonator.palace import (
    driven_to_eigenmode,
    palace_argv,
    validate_eigenmode_config,
)


def driven_fixture() -> dict:
    return {
        "Problem": {"Type": "Driven", "Verbose": 3, "Output": "output/driven"},
        "Model": {"Mesh": "fig1.msh", "L0": 1e-6},
        "Domains": {
            "Materials": [
                {"Attributes": [1], "Permittivity": 1.0},
                {"Attributes": [2], "Permittivity": 9.4, "LossTan": 1e-6},
            ]
        },
        "Boundaries": {
            "Conductivity": [
                {"Attributes": [11], "Conductivity": 1e10, "Thickness": 0.2},
                {"Attributes": [12], "Conductivity": 1e10, "Thickness": 0.2},
            ],
            "Absorbing": {"Attributes": [20], "Order": 2},
            "LumpedPort": [],
        },
        "Solver": {
            "Order": 2,
            "Device": "CPU",
            "Driven": {"Samples": []},
            "Linear": {"Type": "Default", "Tol": 1e-6, "MaxIts": 400},
        },
    }


class PalaceConfigTests(unittest.TestCase):
    def test_driven_to_pec_eigenmode(self) -> None:
        converted = driven_to_eigenmode(
            driven_fixture(), EigenmodeParameters(target_frequency_ghz=4.0)
        )
        self.assertEqual(converted["Problem"]["Type"], "Eigenmode")
        self.assertNotIn("Driven", converted["Solver"])
        self.assertEqual(converted["Solver"]["Eigenmode"]["Target"], 4.0)
        self.assertEqual(converted["Solver"]["Eigenmode"]["TargetUpper"], 12.0)
        self.assertEqual(converted["Boundaries"]["Absorbing"]["Order"], 2)
        self.assertNotIn("Conductivity", converted["Boundaries"])
        self.assertEqual(converted["Boundaries"]["PEC"]["Attributes"], [11, 12])
        validate_eigenmode_config(converted)

    def test_existing_pec_is_preserved(self) -> None:
        source = driven_fixture()
        source["Boundaries"]["PEC"] = {"Attributes": [30]}
        converted = driven_to_eigenmode(source, EigenmodeParameters())
        self.assertEqual(converted["Boundaries"]["PEC"]["Attributes"], [30, 11, 12])

    def test_passive_ports_can_be_preserved(self) -> None:
        source = driven_fixture()
        source["Boundaries"]["LumpedPort"] = [
            {"Index": 1, "R": 50.0, "Attributes": [41], "Excitation": False}
        ]
        converted = driven_to_eigenmode(
            source, EigenmodeParameters(keep_passive_ports=True)
        )
        self.assertEqual(converted["Boundaries"]["LumpedPort"][0]["R"], 50.0)

    def test_missing_metal_is_rejected(self) -> None:
        source = driven_fixture()
        source["Boundaries"].pop("Conductivity")
        with self.assertRaises(ConfigurationError):
            driven_to_eigenmode(source, EigenmodeParameters())

    def test_palace_command(self) -> None:
        self.assertEqual(
            palace_argv("config.json", mpi_processes=6),
            ["palace", "-np", "6", "config.json"],
        )
        self.assertEqual(
            palace_argv("config.json", dry_run=True),
            ["palace", "--dry-run", "config.json"],
        )


if __name__ == "__main__":
    unittest.main()
