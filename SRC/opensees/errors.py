class XaraError(Exception):
    """Base for all OpenSees errors."""

class ModelError(XaraError):
    """The model is in an invalid state (e.g., missing nodes, duplicate tags)."""

class AnalysisError(XaraError):
    """Analysis failed to converge or was misconfigured."""

class ArgumentError(XaraError):
    """Wrong arguments passed to a command."""
