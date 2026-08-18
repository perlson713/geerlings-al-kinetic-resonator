"""Reduced-order direct coupling estimates for separated lumped resonators.

The estimator retains far-field electric- and magnetic-dipole terms and adds
their magnitudes.  It is intended for spacing selection, not as a replacement
for a simultaneous multi-resonator full-wave solve.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

from .geometry import CenterlineGeometry, generate_centerlines
from .parameters import ResonatorParameters


MU0_H_PER_M = 4.0e-7 * math.pi
EPS0_F_PER_M = 8.8541878128e-12


@dataclass(frozen=True, slots=True)
class ResonatorDipole:
    """LC and far-field geometry parameters for one resonator pattern."""

    pattern_id: str
    frequency_hz: float
    impedance_ohm: float
    inductance_h: float
    capacitance_f: float
    electric_centroid_separation_m: float
    magnetic_loop_area_m2: float
    footprint_width_m: float
    footprint_height_m: float

    def to_dict(self) -> dict[str, float | str]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PairCouplingEstimate:
    """Conservative no-cancellation estimate for a resonator pair."""

    first_pattern: str
    second_pattern: str
    separation_m: float
    electric_coupling_hz: float
    magnetic_coupling_hz: float
    total_coupling_hz: float

    def to_dict(self) -> dict[str, float | str]:
        return asdict(self)


def _polyline_length_centroid(line) -> tuple[float, float, float]:
    length = moment_x = moment_y = 0.0
    for first, second in zip(line.points, line.points[1:]):
        segment = math.hypot(second.x - first.x, second.y - first.y)
        length += segment
        moment_x += segment * 0.5 * (first.x + second.x)
        moment_y += segment * 0.5 * (first.y + second.y)
    if length <= 0.0:
        raise ValueError(f"polyline {line.name!r} has no length")
    return length, moment_x, moment_y


def _electrode_centroid(lines: Iterable) -> tuple[float, float]:
    values = [_polyline_length_centroid(line) for line in lines]
    total = sum(value[0] for value in values)
    return (
        sum(value[1] for value in values) / total,
        sum(value[2] for value in values) / total,
    )


def _dipole_geometry(centerlines: CenterlineGeometry) -> tuple[float, float]:
    left = _electrode_centroid(
        (centerlines.capacitor_buses[0], *centerlines.capacitor_fingers[0::2])
    )
    right = _electrode_centroid(
        (centerlines.capacitor_buses[1], *centerlines.capacitor_fingers[1::2])
    )
    electric_separation_um = math.hypot(right[0] - left[0], right[1] - left[1])

    points = centerlines.inductor.points
    signed_area_um2 = 0.5 * sum(
        first.x * second.y - second.x * first.y
        for first, second in zip(points, points[1:] + points[:1])
    )
    return electric_separation_um * 1.0e-6, abs(signed_area_um2) * 1.0e-12


def resonator_dipole(
    pattern_id: str,
    parameters: ResonatorParameters,
    *,
    frequency_ghz: float,
    impedance_ohm: float,
) -> ResonatorDipole:
    """Build an LC dipole descriptor from the generated resonator geometry."""

    if frequency_ghz <= 0.0 or impedance_ohm <= 0.0:
        raise ValueError("frequency and impedance must be positive")
    centerlines = generate_centerlines(parameters)
    electric_separation_m, magnetic_area_m2 = _dipole_geometry(centerlines)
    frequency_hz = frequency_ghz * 1.0e9
    omega = 2.0 * math.pi * frequency_hz
    return ResonatorDipole(
        pattern_id=str(pattern_id),
        frequency_hz=frequency_hz,
        impedance_ohm=impedance_ohm,
        inductance_h=impedance_ohm / omega,
        capacitance_f=1.0 / (omega * impedance_ohm),
        electric_centroid_separation_m=electric_separation_m,
        magnetic_loop_area_m2=magnetic_area_m2,
        footprint_width_m=parameters.nominal_width_um * 1.0e-6,
        footprint_height_m=(
            parameters.capacitor_finger_count * parameters.capacitor_trace_width_um
            + (parameters.capacitor_finger_count - 1) * parameters.g_c_um
            + parameters.g_l_um
            + parameters.inductor_trace_width_um
            + parameters.inductor_turns
            * (parameters.inductor_trace_width_um + parameters.g_l_um)
            + 2.0 * parameters.g_r_um
        )
        * 1.0e-6,
    )


def pair_coupling_estimate(
    first: ResonatorDipole,
    second: ResonatorDipole,
    *,
    separation_m: float,
    effective_relative_permittivity: float,
    electric_orientation_factor: float = 2.0,
    magnetic_orientation_factor: float = 1.0,
) -> PairCouplingEstimate:
    """Estimate |J_e|+|J_m| in hertz for two separated resonators.

    The electric factor of two is the maximum point-dipole orientation factor.
    Adding magnitudes deliberately does not claim electric/magnetic cancellation.
    """

    if separation_m <= 0.0 or effective_relative_permittivity <= 0.0:
        raise ValueError("separation and effective permittivity must be positive")
    if electric_orientation_factor < 0.0 or magnetic_orientation_factor < 0.0:
        raise ValueError("orientation factors must be non-negative")

    inverse_r3 = separation_m**-3
    mutual_capacitance = (
        first.capacitance_f
        * second.capacitance_f
        * first.electric_centroid_separation_m
        * second.electric_centroid_separation_m
        * electric_orientation_factor
        * inverse_r3
        / (4.0 * math.pi * EPS0_F_PER_M * effective_relative_permittivity)
    )
    mutual_inductance = (
        MU0_H_PER_M
        * first.magnetic_loop_area_m2
        * second.magnetic_loop_area_m2
        * magnetic_orientation_factor
        * inverse_r3
        / (4.0 * math.pi)
    )
    mean_frequency = math.sqrt(first.frequency_hz * second.frequency_hz)
    electric_hz = 0.5 * mean_frequency * abs(mutual_capacitance) / math.sqrt(
        first.capacitance_f * second.capacitance_f
    )
    magnetic_hz = 0.5 * mean_frequency * abs(mutual_inductance) / math.sqrt(
        first.inductance_h * second.inductance_h
    )
    return PairCouplingEstimate(
        first_pattern=first.pattern_id,
        second_pattern=second.pattern_id,
        separation_m=separation_m,
        electric_coupling_hz=electric_hz,
        magnetic_coupling_hz=magnetic_hz,
        total_coupling_hz=electric_hz + magnetic_hz,
    )
