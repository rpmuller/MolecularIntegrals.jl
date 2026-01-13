#!/usr/bin/env python
"""
Example script comparing electron repulsion integrals between
MolecularIntegrals.jl and pyquante2.

This script demonstrates how to use the Python wrapper to compute ERIs
and compare them with pyquante2 for validation.

Usage:
    python compare_with_pyquante2.py

Requirements:
    - juliacall: pip install juliacall
    - numpy: pip install numpy
    - pyquante2 (optional): pip install pyquante2
"""

import sys
import time
from pathlib import Path

# Add parent directory to path for local import
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from molecular_integrals import (
    coulomb_repulsion,
    compare_with_pyquante2,
    Atom,
    Basis,
    MolecularIntegrals,
    PGBF,
)


def test_primitive_integrals():
    """Test primitive Gaussian ERIs."""
    print("=" * 60)
    print("Testing Primitive Gaussian ERIs")
    print("=" * 60)

    test_cases = [
        # (description, a_params, b_params, c_params, d_params)
        (
            "(ss|ss) at same center, exponent=1.0",
            ((0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
        ),
        (
            "(ss|ss) separated by 1.4 bohr",
            ((0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
        ),
        (
            "(sp|sp) mixed angular momentum",
            ((0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((0.0, 0.0, 0.0), 1.0, (1, 0, 0), 1.0),  # px
            ((1.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((1.0, 0.0, 0.0), 1.0, (1, 0, 0), 1.0),  # px
        ),
        (
            "(pp|pp) p-orbitals",
            ((0.0, 0.0, 0.0), 1.0, (1, 0, 0), 0.5),
            ((0.0, 0.0, 0.0), 1.0, (0, 1, 0), 0.5),
            ((1.0, 1.0, 0.0), 1.0, (1, 0, 0), 0.5),
            ((1.0, 1.0, 0.0), 1.0, (0, 1, 0), 0.5),
        ),
        (
            "(dd|ss) d-orbital with s",
            ((0.0, 0.0, 0.0), 1.0, (2, 0, 0), 0.8),  # dxx
            ((0.0, 0.0, 0.0), 1.0, (1, 1, 0), 0.8),  # dxy
            ((2.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
            ((2.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0),
        ),
    ]

    results = []
    for desc, a, b, c, d in test_cases:
        print(f"\n{desc}:")
        result = compare_with_pyquante2(*a, *b, *c, *d)
        results.append(result)

    return results


def test_h2_molecule():
    """Test ERIs for H2 molecule with STO-3G basis."""
    print("\n" + "=" * 60)
    print("Testing H2 Molecule with STO-3G Basis")
    print("=" * 60)

    # Create H2 molecule (bond length ~1.4 bohr)
    atoms = [
        Atom(1, (0.0, 0.0, 0.0)),
        Atom(1, (1.4, 0.0, 0.0)),
    ]
    basis = Basis(atoms, "sto3g")

    print(f"\nBasis set: {basis}")
    print(f"Number of basis functions: {len(basis)}")

    mi = MolecularIntegrals()

    # Compute all ERIs
    print("\nComputing all two-electron integrals...")
    start = time.time()
    all_eris = mi.all_twoe_ints(basis)
    julia_time = time.time() - start

    print(f"Time: {julia_time:.4f} s")
    print(f"Number of unique integrals: {len(all_eris)}")
    print(f"Integrals: {all_eris}")

    # Get 4D array for easier inspection
    eris_4d = mi.all_twoe_ints_4d(basis)
    print(f"\n4D tensor shape: {eris_4d.shape}")

    # Show some specific integrals
    print("\nSpecific integrals:")
    print(f"  (00|00) = {eris_4d[0,0,0,0]:.10f}")
    print(f"  (00|11) = {eris_4d[0,0,1,1]:.10f}")
    print(f"  (01|01) = {eris_4d[0,1,0,1]:.10f}")
    print(f"  (11|11) = {eris_4d[1,1,1,1]:.10f}")

    return all_eris, eris_4d


def test_water_molecule():
    """Test ERIs for H2O molecule with STO-3G basis."""
    print("\n" + "=" * 60)
    print("Testing H2O Molecule with STO-3G Basis")
    print("=" * 60)

    # Water molecule geometry (approximate, in bohr)
    # O at origin, H atoms at ~1.8 bohr with 104.5 degree angle
    atoms = [
        Atom(8, (0.0, 0.0, 0.0)),  # O
        Atom(1, (1.43, 1.11, 0.0)),  # H
        Atom(1, (-1.43, 1.11, 0.0)),  # H
    ]
    basis = Basis(atoms, "sto3g")

    print(f"\nBasis set: {basis}")
    print(f"Number of basis functions: {len(basis)}")

    mi = MolecularIntegrals()

    # Compute all ERIs
    print("\nComputing all two-electron integrals...")
    start = time.time()
    all_eris = mi.all_twoe_ints(basis)
    julia_time = time.time() - start

    print(f"Time: {julia_time:.4f} s")
    print(f"Number of unique integrals: {len(all_eris)}")

    # Show statistics
    print(f"\nIntegral statistics:")
    print(f"  Max: {np.max(all_eris):.10f}")
    print(f"  Min: {np.min(all_eris):.10f}")
    print(f"  Mean: {np.mean(all_eris):.10f}")
    print(f"  Std: {np.std(all_eris):.10f}")

    return all_eris


def compare_all_h2_eris_with_pyquante2():
    """Compare all H2 ERIs with pyquante2."""
    print("\n" + "=" * 60)
    print("Full Comparison: H2 ERIs vs pyquante2")
    print("=" * 60)

    try:
        from pyquante2.geo.molecule import molecule
        from pyquante2.basisset import basisset
        from pyquante2.ints.integrals import twoe_integrals

        # Create H2 with pyquante2
        h2_pq = molecule([
            (1, 0.0, 0.0, 0.0),
            (1, 1.4/0.52918, 0.0, 0.0),  # Convert bohr to angstrom
        ], units='Angstrom')
        bfs_pq = basisset(h2_pq, 'sto-3g')
        eris_pq = twoe_integrals(bfs_pq)

        # Create H2 with MolecularIntegrals.jl
        atoms = [Atom(1, (0.0, 0.0, 0.0)), Atom(1, (1.4, 0.0, 0.0))]
        basis = Basis(atoms, "sto3g")
        mi = MolecularIntegrals()
        eris_jl = mi.all_twoe_ints_4d(basis)

        # Compare
        diff = np.abs(eris_pq._2e_ints - eris_jl)
        max_diff = np.max(diff)
        mean_diff = np.mean(diff)

        print(f"\nComparison results:")
        print(f"  Max difference: {max_diff:.2e}")
        print(f"  Mean difference: {mean_diff:.2e}")

        if max_diff < 1e-10:
            print("\n  PASS: Integrals agree to machine precision!")
        elif max_diff < 1e-6:
            print("\n  PASS: Integrals agree to 6 decimal places")
        else:
            print("\n  WARNING: Significant differences found")

        return {"max_diff": max_diff, "mean_diff": mean_diff}

    except ImportError:
        print("\npyquante2 not installed. Install with: pip install pyquante2")
        print("Skipping full comparison.")
        return None


def benchmark_eri_computation():
    """Benchmark ERI computation for various system sizes."""
    print("\n" + "=" * 60)
    print("Benchmarking ERI Computation")
    print("=" * 60)

    mi = MolecularIntegrals()

    test_systems = [
        ("H2 (2 bf)", [Atom(1, (0.0, 0.0, 0.0)), Atom(1, (1.4, 0.0, 0.0))], "sto3g"),
        ("H2O (7 bf)", [
            Atom(8, (0.0, 0.0, 0.0)),
            Atom(1, (1.43, 1.11, 0.0)),
            Atom(1, (-1.43, 1.11, 0.0)),
        ], "sto3g"),
        ("CH4 (9 bf)", [
            Atom(6, (0.0, 0.0, 0.0)),
            Atom(1, (1.19, 1.19, 1.19)),
            Atom(1, (-1.19, -1.19, 1.19)),
            Atom(1, (-1.19, 1.19, -1.19)),
            Atom(1, (1.19, -1.19, -1.19)),
        ], "sto3g"),
    ]

    print(f"\n{'System':<15} {'NBF':>6} {'N_ERI':>10} {'Time (s)':>12}")
    print("-" * 45)

    for name, atoms, basis_name in test_systems:
        basis = Basis(atoms, basis_name)
        n = len(basis)
        n_eri = n * (n + 1) // 2
        n_eri = n_eri * (n_eri + 1) // 2

        start = time.time()
        _ = mi.all_twoe_ints(basis)
        elapsed = time.time() - start

        print(f"{name:<15} {n:>6} {n_eri:>10} {elapsed:>12.6f}")


def main():
    """Run all tests and comparisons."""
    print("\n" + "#" * 60)
    print("# MolecularIntegrals.jl Python Wrapper - Comparison Tests")
    print("#" * 60)

    # Run tests
    test_primitive_integrals()
    test_h2_molecule()
    test_water_molecule()
    compare_all_h2_eris_with_pyquante2()
    benchmark_eri_computation()

    print("\n" + "=" * 60)
    print("All tests completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
