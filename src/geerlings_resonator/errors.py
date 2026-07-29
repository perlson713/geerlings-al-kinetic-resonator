"""Project-specific exceptions."""

from __future__ import annotations


class GeerlingsError(RuntimeError):
    """Base exception for actionable workflow failures."""


class OptionalDependencyError(GeerlingsError):
    """Raised when an explicitly requested optional stage is unavailable."""


class ConfigurationError(GeerlingsError, ValueError):
    """Raised for invalid geometry or simulation settings."""


class ExternalToolError(GeerlingsError):
    """Raised when Gmsh, Palace, or another external program fails."""
