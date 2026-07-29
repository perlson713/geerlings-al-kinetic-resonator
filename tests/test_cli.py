from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from geerlings_resonator.cli import build_parser


class CliTests(unittest.TestCase):
    def test_parser_accepts_geometry_override(self) -> None:
        args = build_parser().parse_args(
            ["build", "--set", "resonator.gC=20", "--step"]
        )
        self.assertEqual(args.command, "build")
        self.assertEqual(args.set_values, ["resonator.gC=20"])
        self.assertTrue(args.step)

    def test_sweep_parser(self) -> None:
        args = build_parser().parse_args(
            ["sweep", "resonator.gC", "3", "5", "10", "--mesh"]
        )
        self.assertEqual(args.values, ["3", "5", "10"])
        self.assertTrue(args.mesh)


if __name__ == "__main__":
    unittest.main()
