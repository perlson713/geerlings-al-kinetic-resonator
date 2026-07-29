"""Reduced-order superconducting aluminium kinetic-inductance correction.

The functions in this module combine a converged PEC eigenfrequency with a
low-temperature thin-film sheet inductance.  This is intentionally not
presented as a full-wave surface-impedance eigenanalysis: the PEC current
profile is reduced to an effective number of squares and a modal impedance.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math


MU0_H_PER_M = 4.0e-7 * math.pi


def _positive(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return number


@dataclass(frozen=True, slots=True)
class AluminiumLondonModel:
    """Low-temperature Al-film fit from Lopez-Nunez et al. (2025)."""

    london_penetration_nm: float = 15.7
    clean_coherence_length_nm: float = 1600.0
    scattering_factor: float = 1.26

    def validate(self) -> None:
        _positive("london_penetration_nm", self.london_penetration_nm)
        _positive("clean_coherence_length_nm", self.clean_coherence_length_nm)
        _positive("scattering_factor", self.scattering_factor)

    def penetration_depth_nm(self, thickness_nm: float) -> float:
        """Return lambda(d)=a*lambda_L*sqrt(xi_0/d), with mean free path l=d."""

        self.validate()
        thickness = _positive("thickness_nm", thickness_nm)
        return (
            self.scattering_factor
            * self.london_penetration_nm
            * math.sqrt(self.clean_coherence_length_nm / thickness)
        )

    def sheet_inductance_h(self, thickness_nm: float) -> float:
        """Return L_square=mu0*lambda*coth(d/lambda), in henries/square."""

        thickness = _positive("thickness_nm", thickness_nm)
        penetration = self.penetration_depth_nm(thickness)
        effective_penetration_nm = penetration / math.tanh(thickness / penetration)
        return MU0_H_PER_M * effective_penetration_nm * 1.0e-9


@dataclass(frozen=True, slots=True)
class ReducedOrderResonator:
    """PEC modal baseline and frozen-current kinetic-inductance reduction."""

    pec_frequency_ghz: float
    modal_impedance_ohm: float
    branch_squares: float
    meander_energy_fraction: float = 0.98

    def validate(self) -> None:
        _positive("pec_frequency_ghz", self.pec_frequency_ghz)
        _positive("modal_impedance_ohm", self.modal_impedance_ohm)
        _positive("branch_squares", self.branch_squares)
        fraction = float(self.meander_energy_fraction)
        if not math.isfinite(fraction) or not 0.0 < fraction <= 1.0:
            raise ValueError("meander_energy_fraction must be in (0, 1]")

    @property
    def effective_squares(self) -> float:
        self.validate()
        return self.branch_squares / self.meander_energy_fraction

    @property
    def geometric_inductance_h(self) -> float:
        self.validate()
        omega = 2.0 * math.pi * self.pec_frequency_ghz * 1.0e9
        return self.modal_impedance_ohm / omega

    @property
    def capacitance_f(self) -> float:
        self.validate()
        omega = 2.0 * math.pi * self.pec_frequency_ghz * 1.0e9
        return 1.0 / (omega * self.modal_impedance_ohm)


@dataclass(frozen=True, slots=True)
class KineticSweepPoint:
    thickness_nm: float
    penetration_depth_nm: float
    sheet_inductance_fh_per_sq: float
    effective_squares: float
    kinetic_inductance_ph: float
    geometric_inductance_ph: float
    kinetic_fraction: float
    frequency_ghz: float
    frequency_shift_mhz: float
    corrected_modal_impedance_ohm: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


def evaluate_point(
    thickness_nm: float,
    *,
    material: AluminiumLondonModel,
    resonator: ReducedOrderResonator,
) -> KineticSweepPoint:
    """Evaluate one thickness using f=f_PEC/sqrt(1+Lk/Lg)."""

    thickness = _positive("thickness_nm", thickness_nm)
    resonator.validate()
    penetration = material.penetration_depth_nm(thickness)
    sheet_h = material.sheet_inductance_h(thickness)
    kinetic_h = sheet_h * resonator.effective_squares
    geometric_h = resonator.geometric_inductance_h
    ratio = kinetic_h / geometric_h
    corrected_ghz = resonator.pec_frequency_ghz / math.sqrt(1.0 + ratio)
    return KineticSweepPoint(
        thickness_nm=thickness,
        penetration_depth_nm=penetration,
        sheet_inductance_fh_per_sq=sheet_h * 1.0e15,
        effective_squares=resonator.effective_squares,
        kinetic_inductance_ph=kinetic_h * 1.0e12,
        geometric_inductance_ph=geometric_h * 1.0e12,
        kinetic_fraction=kinetic_h / (geometric_h + kinetic_h),
        frequency_ghz=corrected_ghz,
        frequency_shift_mhz=(corrected_ghz - resonator.pec_frequency_ghz) * 1.0e3,
        corrected_modal_impedance_ohm=(
            resonator.modal_impedance_ohm * math.sqrt(1.0 + ratio)
        ),
    )


def inclusive_thicknesses(start_nm: float, stop_nm: float, step_nm: float) -> list[float]:
    start = _positive("start_nm", start_nm)
    stop = _positive("stop_nm", stop_nm)
    step = _positive("step_nm", step_nm)
    if stop < start:
        raise ValueError("stop_nm must be >= start_nm")
    count = int(math.floor((stop - start) / step + 1.0e-12))
    values = [start + index * step for index in range(count + 1)]
    if not math.isclose(values[-1], stop, rel_tol=0.0, abs_tol=1.0e-9):
        values.append(stop)
    return values


__all__ = [
    "AluminiumLondonModel",
    "KineticSweepPoint",
    "MU0_H_PER_M",
    "ReducedOrderResonator",
    "evaluate_point",
    "inclusive_thicknesses",
]
