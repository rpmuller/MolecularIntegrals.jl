"""
Python wrapper for MolecularIntegrals.jl electron repulsion integrals.

This module provides a Python interface to the Julia MolecularIntegrals.jl package,
with a function signature compatible with pyquante2 for easy comparison.

Usage:
    from molecular_integrals import coulomb_repulsion, MolecularIntegrals

    # Low-level primitive ERI (pyquante2-compatible signature)
    eri = coulomb_repulsion(
        (xa, ya, za), norma, (la, ma, na), alphaa,
        (xb, yb, zb), normb, (lb, mb, nb), alphab,
        (xc, yc, zc), normc, (lc, mc, nc), alphac,
        (xd, yd, zd), normd, (ld, md, nd), alphad
    )

    # High-level interface
    mi = MolecularIntegrals()
    eri = mi.coulomb(a, b, c, d)  # CGBF objects
    all_eris = mi.all_twoe_ints(basis)

Requirements:
    - Julia (https://julialang.org/)
    - juliacall: pip install juliacall
    - MolecularIntegrals.jl (this package)
"""

import os
import sys
from pathlib import Path
from typing import Tuple, List, Optional, Union
import numpy as np

# Find the package root directory
_PACKAGE_ROOT = Path(__file__).parent.parent.absolute()


def _init_julia():
    """Initialize Julia and load MolecularIntegrals.jl."""
    try:
        from juliacall import Main as jl
    except ImportError:
        raise ImportError(
            "juliacall is required. Install with: pip install juliacall"
        )

    # Add package to Julia's load path and import
    jl.seval(f'push!(LOAD_PATH, "{_PACKAGE_ROOT}")')
    jl.seval('using MolecularIntegrals')

    return jl


# Lazy initialization
_jl = None


def _get_julia():
    """Get or initialize the Julia runtime."""
    global _jl
    if _jl is None:
        _jl = _init_julia()
    return _jl


def coulomb_repulsion(
    origin_a: Tuple[float, float, float], norm_a: float, powers_a: Tuple[int, int, int], alpha_a: float,
    origin_b: Tuple[float, float, float], norm_b: float, powers_b: Tuple[int, int, int], alpha_b: float,
    origin_c: Tuple[float, float, float], norm_c: float, powers_c: Tuple[int, int, int], alpha_c: float,
    origin_d: Tuple[float, float, float], norm_d: float, powers_d: Tuple[int, int, int], alpha_d: float,
) -> float:
    """
    Compute the electron repulsion integral between four primitive Gaussian basis functions.

    This function has the same signature as pyquante2's coulomb_repulsion function for
    easy comparison.

    Parameters
    ----------
    origin_a : tuple of float
        (x, y, z) coordinates of basis function a
    norm_a : float
        Normalization constant for basis function a
    powers_a : tuple of int
        (l, m, n) angular momentum powers for basis function a
    alpha_a : float
        Gaussian exponent for basis function a
    (same for b, c, d)

    Returns
    -------
    float
        The electron repulsion integral <ab|cd>

    Example
    -------
    >>> # s-type functions on hydrogen atoms
    >>> eri = coulomb_repulsion(
    ...     (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    ...     (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    ...     (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    ...     (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    ... )
    """
    jl = _get_julia()

    xa, ya, za = origin_a
    la, ma, na = powers_a
    xb, yb, zb = origin_b
    lb, mb, nb = powers_b
    xc, yc, zc = origin_c
    lc, mc, nc = powers_c
    xd, yd, zd = origin_d
    ld, md, nd = powers_d

    # Call the Julia coulomb function with normalization factors included
    # The Julia function signature is:
    # coulomb(aexpn,ax,ay,az,aI,aJ,aK, bexpn,bx,by,bz,bI,bJ,bK, ...)
    raw_integral = jl.MolecularIntegrals.coulomb(
        alpha_a, xa, ya, za, la, ma, na,
        alpha_b, xb, yb, zb, lb, mb, nb,
        alpha_c, xc, yc, zc, lc, mc, nc,
        alpha_d, xd, yd, zd, ld, md, nd,
    )

    # Apply normalization factors (pyquante2 passes norms separately)
    return norm_a * norm_b * norm_c * norm_d * raw_integral


class PGBF:
    """
    Primitive Gaussian Basis Function wrapper.

    Wraps a Julia PGBF object for use in Python.

    Parameters
    ----------
    expn : float
        Gaussian exponent
    origin : tuple of float
        (x, y, z) coordinates
    powers : tuple of int
        (I, J, K) angular momentum powers
    norm : float, optional
        Normalization constant (default: computed automatically)
    """

    def __init__(
        self,
        expn: float,
        origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
        powers: Tuple[int, int, int] = (0, 0, 0),
        norm: Optional[float] = None,
    ):
        self.jl = _get_julia()
        x, y, z = origin
        I, J, K = powers

        # Use the Julia pgbf helper function which normalizes automatically
        self._julia_obj = self.jl.MolecularIntegrals.pgbf(expn, x, y, z, I, J, K)

        if norm is not None:
            self._julia_obj.norm = norm

    @property
    def expn(self) -> float:
        return float(self._julia_obj.expn)

    @property
    def origin(self) -> Tuple[float, float, float]:
        xyz = self._julia_obj.xyz
        return (float(xyz[1]), float(xyz[2]), float(xyz[3]))

    @property
    def powers(self) -> Tuple[int, int, int]:
        return (int(self._julia_obj.I), int(self._julia_obj.J), int(self._julia_obj.K))

    @property
    def norm(self) -> float:
        return float(self._julia_obj.norm)

    def __repr__(self):
        return f"PGBF(expn={self.expn}, origin={self.origin}, powers={self.powers}, norm={self.norm})"


class CGBF:
    """
    Contracted Gaussian Basis Function wrapper.

    Wraps a Julia CGBF object for use in Python.

    Parameters
    ----------
    origin : tuple of float
        (x, y, z) coordinates
    powers : tuple of int
        (I, J, K) angular momentum powers
    exponents : list of float, optional
        Gaussian exponents for primitives
    coefficients : list of float, optional
        Contraction coefficients
    """

    def __init__(
        self,
        origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
        powers: Tuple[int, int, int] = (0, 0, 0),
        exponents: Optional[List[float]] = None,
        coefficients: Optional[List[float]] = None,
    ):
        self.jl = _get_julia()
        x, y, z = origin
        I, J, K = powers

        # Create empty CGBF
        self._julia_obj = self.jl.MolecularIntegrals.cgbf(x, y, z, I, J, K)

        # Add primitives if provided
        if exponents is not None and coefficients is not None:
            for expn, coef in zip(exponents, coefficients):
                self.jl.MolecularIntegrals.addbf_b(self._julia_obj, expn, coef)

    @classmethod
    def _from_julia(cls, julia_obj):
        """Create a CGBF wrapper from an existing Julia CGBF object."""
        instance = cls.__new__(cls)
        instance.jl = _get_julia()
        instance._julia_obj = julia_obj
        return instance

    @property
    def origin(self) -> Tuple[float, float, float]:
        xyz = self._julia_obj.xyz
        return (float(xyz[1]), float(xyz[2]), float(xyz[3]))

    @property
    def powers(self) -> Tuple[int, int, int]:
        return (int(self._julia_obj.I), int(self._julia_obj.J), int(self._julia_obj.K))

    @property
    def norm(self) -> float:
        return float(self._julia_obj.norm)

    @property
    def exponents(self) -> List[float]:
        return [float(p.expn) for p in self._julia_obj.pgbfs]

    @property
    def coefficients(self) -> List[float]:
        return [float(c) for c in self._julia_obj.coefs]

    def __repr__(self):
        return f"CGBF(origin={self.origin}, powers={self.powers}, nprims={len(self.exponents)})"


class Atom:
    """
    Atom wrapper.

    Parameters
    ----------
    atno : int
        Atomic number
    xyz : tuple of float
        (x, y, z) coordinates in Bohr
    """

    def __init__(self, atno: int, xyz: Tuple[float, float, float]):
        self.jl = _get_julia()
        x, y, z = xyz
        # Create Julia MVector for coordinates
        self._julia_obj = self.jl.MolecularIntegrals.Atom(
            atno,
            self.jl.StaticArrays.MVector(x, y, z)
        )

    @classmethod
    def _from_julia(cls, julia_obj):
        """Create an Atom wrapper from an existing Julia Atom object."""
        instance = cls.__new__(cls)
        instance.jl = _get_julia()
        instance._julia_obj = julia_obj
        return instance

    @property
    def atno(self) -> int:
        return int(self._julia_obj.atno)

    @property
    def xyz(self) -> Tuple[float, float, float]:
        xyz = self._julia_obj.xyz
        return (float(xyz[1]), float(xyz[2]), float(xyz[3]))

    def __repr__(self):
        return f"Atom(atno={self.atno}, xyz={self.xyz})"


class Basis:
    """
    Basis set wrapper.

    Parameters
    ----------
    atoms : list of Atom
        List of atoms
    name : str
        Basis set name (e.g., 'sto3g', '6-31g')
    """

    def __init__(self, atoms: List[Atom], name: str = "sto3g"):
        self.jl = _get_julia()

        # Convert Python atoms to Julia vector
        julia_atoms = self.jl.seval("Atom[]")
        for atom in atoms:
            self.jl.push_b(julia_atoms, atom._julia_obj)

        # Build the basis
        self._julia_obj = self.jl.MolecularIntegrals.build_basis(julia_atoms, name)

    @classmethod
    def _from_julia(cls, julia_obj):
        """Create a Basis wrapper from an existing Julia Basis object."""
        instance = cls.__new__(cls)
        instance.jl = _get_julia()
        instance._julia_obj = julia_obj
        return instance

    def __len__(self) -> int:
        return int(self.jl.length(self._julia_obj))

    def __getitem__(self, i: int) -> CGBF:
        # Julia is 1-indexed
        return CGBF._from_julia(self._julia_obj.cgbfs[i + 1])

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __repr__(self):
        return f"Basis(nbf={len(self)})"


class MolecularIntegrals:
    """
    High-level interface to MolecularIntegrals.jl.

    This class provides convenient access to electron repulsion integrals
    and other molecular integrals.

    Example
    -------
    >>> from molecular_integrals import MolecularIntegrals, Atom, Basis
    >>>
    >>> # Create H2 molecule
    >>> atoms = [Atom(1, (0.0, 0.0, 0.0)), Atom(1, (1.4, 0.0, 0.0))]
    >>> basis = Basis(atoms, "sto3g")
    >>>
    >>> mi = MolecularIntegrals()
    >>>
    >>> # Compute single ERI
    >>> eri = mi.coulomb(basis[0], basis[0], basis[1], basis[1])
    >>>
    >>> # Compute all ERIs
    >>> all_eris = mi.all_twoe_ints(basis)
    """

    def __init__(self):
        self.jl = _get_julia()

    def coulomb_pgbf(self, a: PGBF, b: PGBF, c: PGBF, d: PGBF) -> float:
        """
        Compute ERI between four primitive Gaussian basis functions.

        Parameters
        ----------
        a, b, c, d : PGBF
            Primitive Gaussian basis functions

        Returns
        -------
        float
            The electron repulsion integral <ab|cd>
        """
        return float(self.jl.MolecularIntegrals.coulomb(
            a._julia_obj, b._julia_obj, c._julia_obj, d._julia_obj
        ))

    def coulomb(self, a: CGBF, b: CGBF, c: CGBF, d: CGBF) -> float:
        """
        Compute ERI between four contracted Gaussian basis functions.

        Parameters
        ----------
        a, b, c, d : CGBF
            Contracted Gaussian basis functions

        Returns
        -------
        float
            The electron repulsion integral <ab|cd>
        """
        return float(self.jl.MolecularIntegrals.coulomb(
            a._julia_obj, b._julia_obj, c._julia_obj, d._julia_obj
        ))

    def all_twoe_ints(self, basis: Basis, method: str = "default") -> np.ndarray:
        """
        Compute all two-electron integrals for a basis set.

        Parameters
        ----------
        basis : Basis
            The basis set
        method : str
            Method to use: 'default', 'hgp', or 'rys'

        Returns
        -------
        np.ndarray
            1D array of unique ERIs using 8-fold symmetry
        """
        if method == "default":
            eri_func = self.jl.MolecularIntegrals.coulomb
        elif method == "hgp":
            # Use HGP method if available
            eri_func = getattr(self.jl.MolecularIntegrals, 'coulomb',
                              self.jl.MolecularIntegrals.coulomb)
        elif method == "rys":
            eri_func = self.jl.MolecularIntegrals.coulomb_rys
        else:
            raise ValueError(f"Unknown method: {method}")

        julia_ints = self.jl.MolecularIntegrals.all_twoe_ints(
            basis._julia_obj.cgbfs, eri_func
        )

        return np.array([float(x) for x in julia_ints])

    def all_twoe_ints_4d(self, basis: Basis, method: str = "default") -> np.ndarray:
        """
        Compute all two-electron integrals as a 4D array.

        This is less memory-efficient but more convenient for comparison
        with pyquante2.

        Parameters
        ----------
        basis : Basis
            The basis set
        method : str
            Method to use: 'default', 'hgp', or 'rys'

        Returns
        -------
        np.ndarray
            4D array of shape (nbf, nbf, nbf, nbf)
        """
        n = len(basis)
        ints_1d = self.all_twoe_ints(basis, method)
        ints_4d = np.zeros((n, n, n, n))

        # Expand 1D packed array to 4D using symmetry
        # The Julia iindex function uses 1-based indexing, so we implement it here
        def triangle(i):
            return i * (i + 1) // 2

        def iindex(i, j, k, l):
            ij = triangle(max(i, j)) + min(i, j)
            kl = triangle(max(k, l)) + min(k, l)
            return triangle(max(ij, kl)) + min(ij, kl)

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        idx = iindex(i, j, k, l)
                        ints_4d[i, j, k, l] = ints_1d[idx]

        return ints_4d

    def overlap(self, a: CGBF, b: CGBF) -> float:
        """Compute overlap integral between two CGBFs."""
        return float(self.jl.MolecularIntegrals.overlap(a._julia_obj, b._julia_obj))

    def kinetic(self, a: CGBF, b: CGBF) -> float:
        """Compute kinetic energy integral between two CGBFs."""
        return float(self.jl.MolecularIntegrals.kinetic(a._julia_obj, b._julia_obj))


# Convenience functions for quick comparisons

def compare_with_pyquante2(
    origin_a, norm_a, powers_a, alpha_a,
    origin_b, norm_b, powers_b, alpha_b,
    origin_c, norm_c, powers_c, alpha_c,
    origin_d, norm_d, powers_d, alpha_d,
    verbose: bool = True
) -> dict:
    """
    Compare ERI values between MolecularIntegrals.jl and pyquante2.

    Returns a dictionary with both values and the difference.
    """
    julia_eri = coulomb_repulsion(
        origin_a, norm_a, powers_a, alpha_a,
        origin_b, norm_b, powers_b, alpha_b,
        origin_c, norm_c, powers_c, alpha_c,
        origin_d, norm_d, powers_d, alpha_d,
    )

    try:
        # Try to import pyquante2
        try:
            from pyquante2.ints.hgp import coulomb_repulsion as pq_coulomb
        except ImportError:
            from pyquante2.ints.rys import coulomb_repulsion as pq_coulomb

        pyquante_eri = pq_coulomb(
            origin_a, norm_a, powers_a, alpha_a,
            origin_b, norm_b, powers_b, alpha_b,
            origin_c, norm_c, powers_c, alpha_c,
            origin_d, norm_d, powers_d, alpha_d,
        )

        diff = abs(julia_eri - pyquante_eri)

        if verbose:
            print(f"MolecularIntegrals.jl: {julia_eri:.12f}")
            print(f"pyquante2:            {pyquante_eri:.12f}")
            print(f"Difference:           {diff:.2e}")

        return {
            "julia": julia_eri,
            "pyquante2": pyquante_eri,
            "difference": diff,
        }

    except ImportError:
        if verbose:
            print(f"MolecularIntegrals.jl: {julia_eri:.12f}")
            print("pyquante2 not installed - cannot compare")

        return {
            "julia": julia_eri,
            "pyquante2": None,
            "difference": None,
        }


if __name__ == "__main__":
    # Simple test
    print("Testing MolecularIntegrals.jl Python wrapper...")

    # Test primitive ERI (s-type functions)
    print("\n1. Testing primitive ERI (s-type):")
    eri = coulomb_repulsion(
        (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
        (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
        (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
        (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    )
    print(f"   (ss|ss) = {eri:.10f}")

    print("\n2. Testing PGBF wrapper:")
    p = PGBF(1.0, (0.0, 0.0, 0.0), (0, 0, 0))
    print(f"   {p}")

    print("\n3. Testing Atom and Basis:")
    atoms = [Atom(1, (0.0, 0.0, 0.0)), Atom(1, (1.4, 0.0, 0.0))]
    basis = Basis(atoms, "sto3g")
    print(f"   H2 with STO-3G: {basis}")

    print("\n4. Testing contracted ERI:")
    mi = MolecularIntegrals()
    eri = mi.coulomb(basis[0], basis[0], basis[1], basis[1])
    print(f"   (11|22) = {eri:.10f}")

    print("\n5. Testing all_twoe_ints:")
    all_eris = mi.all_twoe_ints(basis)
    print(f"   Number of unique ERIs: {len(all_eris)}")
    print(f"   ERIs: {all_eris}")

    print("\n6. Compare with pyquante2:")
    compare_with_pyquante2(
        (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
        (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
        (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
        (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    )

    print("\nAll tests completed!")
