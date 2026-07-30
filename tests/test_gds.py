from __future__ import annotations

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from geerlings_resonator.config import load_project
from geerlings_resonator.gds import write_gds
from geerlings_resonator.geometry import generate_layout


@unittest.skipUnless(importlib.util.find_spec("gdsfactory"), "gdsfactory not installed")
class GDSExportTests(unittest.TestCase):
    def test_current_gdsfactory_can_write_numeric_layer_layout(self) -> None:
        layout = generate_layout(
            load_project("configs/design_fig1_square_9ghz_al200.toml").resonator
        )
        with TemporaryDirectory() as directory:
            output = Path(directory) / "layout.gds"
            write_gds(layout, output, cell_name="fig1_gds_test", include_ports=False)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 1024)


if __name__ == "__main__":
    unittest.main()
