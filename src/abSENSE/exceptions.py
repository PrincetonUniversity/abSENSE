"""Custom exceptions for the package."""
from __future__ import annotations


class AbsenseException(Exception):
    """Base exception class."""


class MissingSpeciesException(AbsenseException):
    """A species is missing from the analysis."""


class MissingGeneException(AbsenseException):
    """A gene is missing from the analysis."""
