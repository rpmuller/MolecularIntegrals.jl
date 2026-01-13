"""
Python wrapper for MolecularIntegrals.jl

Provides access to electron repulsion integrals computed by the Julia package,
with a pyquante2-compatible interface for easy comparison.
"""

from .molecular_integrals import (
    coulomb_repulsion,
    compare_with_pyquante2,
    PGBF,
    CGBF,
    Atom,
    Basis,
    MolecularIntegrals,
)

__all__ = [
    "coulomb_repulsion",
    "compare_with_pyquante2",
    "PGBF",
    "CGBF",
    "Atom",
    "Basis",
    "MolecularIntegrals",
]

__version__ = "0.1.0"
