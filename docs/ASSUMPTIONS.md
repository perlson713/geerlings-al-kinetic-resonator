# Reproduction boundary

This implementation is a topology-accurate, parameterized reconstruction of
Fig. 1, not a recovered Yale mask file.

The supplied `configs/design_a.toml` is the requested isolated-resonator
variant: the upper CPW feedline and its ports are disabled. CPW parameters are
retained only for an optional Driven/S21 comparison model.

## Article-derived values

- Design A: `gC=10 µm`, `gL=20 µm`, `gR=10 µm`, `w=5 µm`, `Z0=100 Ω`
- Design B: `gC=20 µm`, `gL=5 µm`, `gR=10 µm`, `w=10 µm`, `Z0=200 Ω`
- Design C: `gC=80 µm`, `gL=10 µm`, `gR=10 µm`, `wL=10 µm`,
  `wC=40 µm`, `Z0=300 Ω`
- 200 nm etched Nb on c-plane sapphire
- tested resonators at 5-8 GHz
- Fig. 1 horizontal resonator scale: 300 µm

## Figure-derived counts

- 14 IDC fingers: 7 on each electrode
- 21 horizontal inductor spans: 10 per side plus one top bridge

## Initial assumptions that must be calibrated

- exact finger and meander horizontal lengths
- substrate thickness, permittivity tensor, and loss tangent
- optional CPW center width/gap and package geometry
- optional finite CPW open-end clearance and ground-reconnection bridge
- optional ground shield width and feedline-parallel coupling length
- enclosure height and absorbing-boundary distance
- Nb surface impedance/London penetration depth below `Tc`

The TOML file exposes every item above. A fabrication-targeted result should
replace the defaults with process measurements and perform mesh/box
convergence plus a Driven-versus-Eigenmode cross-check.
