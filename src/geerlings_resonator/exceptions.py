"""Package-specific exceptions."""


class MissingOptionalDependencyError(ImportError):
    """Raised when an explicitly requested export backend is unavailable."""


__all__ = ["MissingOptionalDependencyError"]
