from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from geerlings_resonator.results import (
    filter_frequency_window,
    format_modes,
    read_eigenmodes,
)


SAMPLE = """        m,                Re{f} (GHz),                Im{f} (GHz),                          Q,              Error (Bkwd.),               Error (Abs.)
 1.00e+00,        +4.099115457610e+00,        +1.104722583251e-04,        +1.855269151390e+04,        +3.275563953950e-14,        +1.108615314934e-08
 2.00e+00,        +5.603265962190e+00,        +3.541371834327e-04,        +7.911151716785e+03,        +3.490195627919e-14,        +1.181261238013e-08
"""


class ResultTests(unittest.TestCase):
    def test_read_filter_and_format(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "eig.csv"
            path.write_text(SAMPLE, encoding="utf-8")
            modes = read_eigenmodes(path)
        self.assertEqual([mode.mode for mode in modes], [1, 2])
        self.assertAlmostEqual(modes[0].frequency_ghz, 4.099115457610)
        selected = filter_frequency_window(modes, 5.0, 6.0)
        self.assertEqual([mode.mode for mode in selected], [2])
        table = format_modes(selected)
        self.assertIn("Re(f) GHz", table)
        self.assertIn("5.60326596", table)


if __name__ == "__main__":
    unittest.main()
