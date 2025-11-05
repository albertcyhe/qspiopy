"""Domain-specific exceptions for the offline runtime."""

from __future__ import annotations


class OfflineError(RuntimeError):
    """Base class for offline runtime errors."""


class ConfigError(OfflineError):
    """Raised when solver or scenario configuration is invalid."""


class SnapshotError(OfflineError):
    """Raised when frozen snapshot assets are missing or malformed."""


class NumericsError(OfflineError):
    """Raised when the numerical solver fails to converge."""


class AlignmentFail(OfflineError):
    """Raised when semantic alignment checks exceed tolerances."""


__all__ = [
    "OfflineError",
    "ConfigError",
    "SnapshotError",
    "NumericsError",
    "AlignmentFail",
]
