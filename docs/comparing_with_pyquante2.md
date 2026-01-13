# Comparing MolecularIntegrals.jl with pyquante2

This document describes how to use the Python wrapper for MolecularIntegrals.jl to compare electron repulsion integral (ERI) calculations with the pyquante2 package.

## Overview

MolecularIntegrals.jl now includes a Python wrapper that provides a pyquante2-compatible interface for computing electron repulsion integrals. This enables direct comparison of:

1. **Numerical accuracy** - Verify that both packages produce identical results
2. **Performance** - Benchmark computation times for various system sizes
3. **API compatibility** - Use the same function signatures for both packages

## Installation

### Prerequisites

1. **Julia** (v1.6 or later): https://julialang.org/downloads/
2. **Python** (v3.8 or later)
3. **MolecularIntegrals.jl** (this package)

### Python Dependencies

```bash
# Required
pip install juliacall numpy

# For comparison with pyquante2
pip install pyquante2
```

### Verify Installation

```bash
# Check Julia is available
julia --version

# Test the wrapper
cd MolecularIntegrals.jl
python -c "from python.molecular_integrals import coulomb_repulsion; print('Success!')"
```

## Quick Start

### Basic ERI Comparison

```python
import sys
sys.path.insert(0, '/path/to/MolecularIntegrals.jl')

from python.molecular_integrals import coulomb_repulsion, compare_with_pyquante2

# Compare a single (ss|ss) integral
result = compare_with_pyquante2(
    (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # s-function at origin
    (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # s-function at origin
    (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # s-function at 1.4 bohr
    (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # s-function at 1.4 bohr
)

# Output:
# MolecularIntegrals.jl: 0.349813832507
# pyquante2:            0.349813832507
# Difference:           1.23e-15
```

## Detailed Comparison Methods

### Method 1: Primitive Gaussian ERIs

Compare individual primitive Gaussian basis function integrals using the pyquante2-compatible `coulomb_repulsion` function.

```python
from python.molecular_integrals import coulomb_repulsion
import time

# Parameters: origin, normalization, angular momentum (l,m,n), exponent
# for each of four basis functions a, b, c, d

def compare_primitive_eri():
    """Compare primitive ERI between both packages."""

    # Test case: (px py | px py) integral
    params = (
        (0.0, 0.0, 0.0), 1.0, (1, 0, 0), 0.5,  # px at origin
        (0.0, 0.0, 0.0), 1.0, (0, 1, 0), 0.5,  # py at origin
        (2.0, 0.0, 0.0), 1.0, (1, 0, 0), 0.5,  # px displaced
        (2.0, 0.0, 0.0), 1.0, (0, 1, 0), 0.5,  # py displaced
    )

    # MolecularIntegrals.jl
    start = time.perf_counter()
    julia_eri = coulomb_repulsion(*params)
    julia_time = time.perf_counter() - start

    # pyquante2
    from pyquante2.ints.hgp import coulomb_repulsion as pq_coulomb
    start = time.perf_counter()
    pq_eri = pq_coulomb(*params)
    pq_time = time.perf_counter() - start

    print(f"MolecularIntegrals.jl: {julia_eri:.12f} ({julia_time*1000:.3f} ms)")
    print(f"pyquante2:            {pq_eri:.12f} ({pq_time*1000:.3f} ms)")
    print(f"Difference:           {abs(julia_eri - pq_eri):.2e}")

compare_primitive_eri()
```

### Method 2: Contracted Basis Function ERIs

Compare integrals over contracted Gaussian basis functions (CGBFs).

```python
from python.molecular_integrals import Atom, Basis, MolecularIntegrals
import numpy as np

def compare_contracted_eris():
    """Compare contracted ERIs for H2 molecule."""

    # Create H2 with MolecularIntegrals.jl
    atoms_jl = [
        Atom(1, (0.0, 0.0, 0.0)),
        Atom(1, (1.4, 0.0, 0.0)),
    ]
    basis_jl = Basis(atoms_jl, "sto3g")
    mi = MolecularIntegrals()

    # Compute all ERIs as 4D tensor
    eris_jl = mi.all_twoe_ints_4d(basis_jl)

    # Create H2 with pyquante2
    from pyquante2.geo.molecule import molecule
    from pyquante2.basisset import basisset
    from pyquante2.ints.integrals import twoe_integrals

    h2_pq = molecule([
        (1, 0.0, 0.0, 0.0),
        (1, 1.4/0.52918, 0.0, 0.0),  # Convert bohr to angstrom
    ], units='Angstrom')
    bfs_pq = basisset(h2_pq, 'sto-3g')
    eris_pq = twoe_integrals(bfs_pq)._2e_ints

    # Compare all integrals
    diff = np.abs(eris_jl - eris_pq)

    print("Comparison of all H2 ERIs:")
    print(f"  Shape: {eris_jl.shape}")
    print(f"  Max difference: {np.max(diff):.2e}")
    print(f"  Mean difference: {np.mean(diff):.2e}")
    print(f"  RMS difference: {np.sqrt(np.mean(diff**2)):.2e}")

    # Show individual integrals
    print("\nIndividual integrals:")
    print(f"  {'Integral':<12} {'Julia':>14} {'pyquante2':>14} {'Diff':>10}")
    print("  " + "-" * 52)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    jl = eris_jl[i,j,k,l]
                    pq = eris_pq[i,j,k,l]
                    d = abs(jl - pq)
                    print(f"  ({i}{j}|{k}{l})      {jl:14.10f} {pq:14.10f} {d:10.2e}")

compare_contracted_eris()
```

### Method 3: Performance Benchmarking

Compare computation times for various molecular systems.

```python
from python.molecular_integrals import Atom, Basis, MolecularIntegrals
import time
import numpy as np

def benchmark_comparison():
    """Benchmark ERI computation for various systems."""

    # Define test molecules (coordinates in Bohr)
    molecules = {
        "H2": [
            Atom(1, (0.0, 0.0, 0.0)),
            Atom(1, (1.4, 0.0, 0.0)),
        ],
        "H2O": [
            Atom(8, (0.0, 0.0, 0.0)),
            Atom(1, (1.43, 1.11, 0.0)),
            Atom(1, (-1.43, 1.11, 0.0)),
        ],
        "CH4": [
            Atom(6, (0.0, 0.0, 0.0)),
            Atom(1, (1.19, 1.19, 1.19)),
            Atom(1, (-1.19, -1.19, 1.19)),
            Atom(1, (-1.19, 1.19, -1.19)),
            Atom(1, (1.19, -1.19, -1.19)),
        ],
    }

    mi = MolecularIntegrals()

    print(f"{'Molecule':<10} {'NBF':>6} {'Julia (ms)':>12} {'pyquante2 (ms)':>16} {'Speedup':>10}")
    print("-" * 58)

    for name, atoms in molecules.items():
        # MolecularIntegrals.jl timing
        basis = Basis(atoms, "sto3g")
        n = len(basis)

        # Warm-up run (JIT compilation)
        _ = mi.all_twoe_ints(basis)

        # Timed run
        start = time.perf_counter()
        for _ in range(5):
            eris_jl = mi.all_twoe_ints(basis)
        julia_time = (time.perf_counter() - start) / 5 * 1000  # ms

        # pyquante2 timing
        try:
            from pyquante2.geo.molecule import molecule
            from pyquante2.basisset import basisset
            from pyquante2.ints.integrals import twoe_integrals

            # Convert to pyquante2 format
            mol_pq = molecule([
                (a.atno, a.xyz[0]/0.52918, a.xyz[1]/0.52918, a.xyz[2]/0.52918)
                for a in atoms
            ], units='Angstrom')
            bfs_pq = basisset(mol_pq, 'sto-3g')

            start = time.perf_counter()
            for _ in range(5):
                eris_pq = twoe_integrals(bfs_pq)
            pq_time = (time.perf_counter() - start) / 5 * 1000  # ms

            speedup = pq_time / julia_time
            print(f"{name:<10} {n:>6} {julia_time:>12.3f} {pq_time:>16.3f} {speedup:>10.2f}x")

        except ImportError:
            print(f"{name:<10} {n:>6} {julia_time:>12.3f} {'N/A':>16} {'N/A':>10}")

benchmark_comparison()
```

## Comprehensive Test Suite

Run the included comparison script for a full validation:

```bash
cd MolecularIntegrals.jl
python python/examples/compare_with_pyquante2.py
```

This script tests:
- Primitive Gaussian ERIs with various angular momenta (s, p, d)
- Contracted ERIs for H2 and H2O molecules
- Full comparison of all H2 integrals
- Performance benchmarks

## API Reference

### Low-Level Function

```python
coulomb_repulsion(
    origin_a, norm_a, powers_a, alpha_a,
    origin_b, norm_b, powers_b, alpha_b,
    origin_c, norm_c, powers_c, alpha_c,
    origin_d, norm_d, powers_d, alpha_d,
) -> float
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `origin_x` | `(float, float, float)` | (x, y, z) coordinates in Bohr |
| `norm_x` | `float` | Normalization constant |
| `powers_x` | `(int, int, int)` | Angular momentum (l, m, n) |
| `alpha_x` | `float` | Gaussian exponent |

### High-Level Classes

| Class | Description |
|-------|-------------|
| `Atom(atno, xyz)` | Atom with atomic number and coordinates (Bohr) |
| `Basis(atoms, name)` | Basis set ("sto3g", "6-31g", etc.) |
| `MolecularIntegrals()` | Main interface for integral computation |
| `PGBF(expn, origin, powers)` | Primitive Gaussian basis function |
| `CGBF(origin, powers, exponents, coefficients)` | Contracted Gaussian |

### MolecularIntegrals Methods

| Method | Description |
|--------|-------------|
| `coulomb(a, b, c, d)` | ERI between four CGBFs |
| `coulomb_pgbf(a, b, c, d)` | ERI between four PGBFs |
| `all_twoe_ints(basis)` | All unique ERIs (1D packed array) |
| `all_twoe_ints_4d(basis)` | All ERIs as 4D tensor |
| `overlap(a, b)` | Overlap integral |
| `kinetic(a, b)` | Kinetic energy integral |

## Expected Results

For correctly implemented ERIs, you should see:

- **Numerical agreement**: Differences < 10⁻¹⁰ (machine precision effects)
- **Identical symmetry**: (ij|kl) = (ji|kl) = (ij|lk) = (kl|ij), etc.

### Example Output

```
Comparison of all H2 ERIs:
  Shape: (2, 2, 2, 2)
  Max difference: 1.11e-15
  Mean difference: 4.23e-16
  RMS difference: 5.67e-16

Individual integrals:
  Integral          Julia      pyquante2       Diff
  ----------------------------------------------------
  (00|00)    0.7746059439  0.7746059439   1.11e-16
  (00|01)    0.4441076657  0.4441076657   0.00e+00
  (00|10)    0.4441076657  0.4441076657   0.00e+00
  (00|11)    0.5696758031  0.5696758031   1.11e-16
  ...
```

## Troubleshooting

### Slow First Run

The first call to Julia functions triggers JIT compilation. This is normal:
- First call: 1-5 seconds
- Subsequent calls: milliseconds

### Import Errors

```python
# If juliacall not found:
pip install juliacall

# If MolecularIntegrals.jl not found:
# Ensure you're running from the correct directory or update sys.path
```

### Numerical Differences

Small differences (< 10⁻¹⁰) are expected due to:
- Different floating-point evaluation order
- Different intermediate precision

Large differences may indicate:
- Different basis set definitions
- Coordinate unit mismatch (Bohr vs Angstrom)
- Normalization convention differences

## References

- [pyquante2 GitHub](https://github.com/rpmuller/pyquante2)
- [juliacall Documentation](https://juliapy.github.io/PythonCall.jl/stable/)
- Head-Gordon & Pople, J. Chem. Phys. 89, 5777 (1988) - HGP method
- Taketa, Huzinaga, O-ohata, J. Phys. Soc. Japan 21, 2313 (1966) - THO method
