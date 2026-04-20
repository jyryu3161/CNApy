"""Exception hierarchy for the ecmodel package.

Design Ref: §6.2 — EcModelError base with specific subclasses.
"""

from __future__ import annotations


class EcModelError(Exception):
    """Base exception for ecmodel-related errors."""


class EcYamlError(EcModelError):
    """GECKO YAML parsing or serialization error."""


class GprParseError(EcModelError):
    """Cannot parse a gene_reaction_rule string."""

    def __init__(self, gpr: str, reason: str = ""):
        self.gpr = gpr
        self.reason = reason
        msg = f"Failed to parse GPR: {gpr!r}"
        if reason:
            msg += f" — {reason}"
        super().__init__(msg)


class SchemaVersionError(EcModelError):
    """Unknown or unsupported ECModelData schema version."""
