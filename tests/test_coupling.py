from __future__ import annotations

import unittest

from geerlings_resonator.config import load_project
from geerlings_resonator.coupling import (
    pair_coupling_estimate,
    resonator_dipole,
)


class DirectCouplingTests(unittest.TestCase):
    def test_far_field_pair_coupling_scales_as_inverse_cube(self) -> None:
        config = load_project("configs/design_fig1_square_9ghz_al200.toml")
        first = resonator_dipole(
            "P01", config.resonator, frequency_ghz=8.965, impedance_ohm=100.0
        )
        second = resonator_dipole(
            "P02",
            config.resonator.with_updates(nominal_width_um=425.5),
            frequency_ghz=8.975,
            impedance_ohm=100.0,
        )
        near = pair_coupling_estimate(
            first,
            second,
            separation_m=2.0e-3,
            effective_relative_permittivity=5.2,
        )
        far = pair_coupling_estimate(
            first,
            second,
            separation_m=4.0e-3,
            effective_relative_permittivity=5.2,
        )
        self.assertAlmostEqual(near.total_coupling_hz / far.total_coupling_hz, 8.0)
        self.assertGreater(far.electric_coupling_hz, far.magnetic_coupling_hz)


if __name__ == "__main__":
    unittest.main()
