# Python Wrapper for MolecularIntegrals.jl

This directory contains a Python wrapper for the electron repulsion integrals (ERIs)
computed by MolecularIntegrals.jl. The wrapper provides a pyquante2-compatible interface
for easy comparison between the two packages.

## Installation

### Prerequisites

1. **Julia**: Install from https://julialang.org/downloads/
2. **MolecularIntegrals.jl**: This package (the parent directory)

### Python Dependencies

```bash
pip install juliacall numpy
```

For comparison with pyquante2:
```bash
pip install pyquante2
```

### Usage

You can either:

1. **Use directly** by adding the `python` directory to your Python path
2. **Install as a package** with `pip install -e .` from this directory

## Quick Start

### Low-level primitive ERI (pyquante2-compatible)

```python
from molecular_integrals import coulomb_repulsion

# Compute (ss|ss) integral
eri = coulomb_repulsion(
    (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # function a: origin, norm, powers, exponent
    (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # function b
    (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # function c
    (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,  # function d
)
print(f"ERI = {eri}")
```

### High-level interface with basis sets

```python
from molecular_integrals import Atom, Basis, MolecularIntegrals

# Create H2 molecule
atoms = [
    Atom(1, (0.0, 0.0, 0.0)),    # H at origin
    Atom(1, (1.4, 0.0, 0.0)),    # H at 1.4 bohr
]

# Build STO-3G basis
basis = Basis(atoms, "sto3g")

# Create integral engine
mi = MolecularIntegrals()

# Compute single ERI
eri = mi.coulomb(basis[0], basis[0], basis[1], basis[1])

# Compute all ERIs
all_eris = mi.all_twoe_ints(basis)

# Get 4D tensor for easier indexing
eris_4d = mi.all_twoe_ints_4d(basis)
print(f"(00|11) = {eris_4d[0,0,1,1]}")
```

### Comparison with pyquante2

```python
from molecular_integrals import compare_with_pyquante2

# Compare a single integral
result = compare_with_pyquante2(
    (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    (0.0, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
    (1.4, 0.0, 0.0), 1.0, (0, 0, 0), 1.0,
)
# Output:
# MolecularIntegrals.jl: 0.349813832507
# pyquante2:            0.349813832507
# Difference:           1.23e-15
```

## API Reference

### `coulomb_repulsion`

Compute ERI between four primitive Gaussians (pyquante2-compatible signature).

```python
coulomb_repulsion(
    origin_a, norm_a, powers_a, alpha_a,
    origin_b, norm_b, powers_b, alpha_b,
    origin_c, norm_c, powers_c, alpha_c,
    origin_d, norm_d, powers_d, alpha_d,
) -> float
```

### `PGBF`

Primitive Gaussian Basis Function wrapper.

```python
p = PGBF(expn=1.0, origin=(0., 0., 0.), powers=(0, 0, 0))
```

### `CGBF`

Contracted Gaussian Basis Function wrapper.

```python
c = CGBF(
    origin=(0., 0., 0.),
    powers=(0, 0, 0),
    exponents=[3.42525, 0.62391, 0.16886],
    coefficients=[0.15433, 0.53533, 0.44463],
)
```

### `Atom`

Atom representation.

```python
h = Atom(atno=1, xyz=(0.0, 0.0, 0.0))  # Hydrogen at origin (coordinates in Bohr)
```

### `Basis`

Basis set for a molecule.

```python
basis = Basis(atoms, "sto3g")  # Supported: "sto3g", "6-31g", etc.
```

### `MolecularIntegrals`

High-level interface to the integral engine.

```python
mi = MolecularIntegrals()

# ERIs
eri = mi.coulomb(cgbf_a, cgbf_b, cgbf_c, cgbf_d)
eri = mi.coulomb_pgbf(pgbf_a, pgbf_b, pgbf_c, pgbf_d)

# All ERIs
all_eris = mi.all_twoe_ints(basis)
all_eris_4d = mi.all_twoe_ints_4d(basis)

# One-electron integrals
s = mi.overlap(cgbf_a, cgbf_b)
t = mi.kinetic(cgbf_a, cgbf_b)
```

## Examples

See the `examples/` directory for complete examples:

- `compare_with_pyquante2.py`: Comprehensive comparison between MolecularIntegrals.jl and pyquante2

## Running Tests

```bash
python examples/compare_with_pyquante2.py
```

Or run the module directly:

```bash
python -m molecular_integrals
```

## Performance Notes

The first call to any function will be slow due to Julia JIT compilation.
Subsequent calls will be fast. For best performance:

1. Use `all_twoe_ints()` for batch computation (uses Julia's threading)
2. Avoid calling `coulomb_repulsion()` in tight Python loops
3. Consider using the Julia package directly for performance-critical applications

## Troubleshooting

### Julia not found

Ensure Julia is in your PATH:
```bash
julia --version
```

### MolecularIntegrals.jl not found

The wrapper automatically adds the parent directory to Julia's LOAD_PATH.
If you moved the Python files, update `_PACKAGE_ROOT` in `molecular_integrals.py`.

### Slow first call

This is expected due to Julia's JIT compilation. Subsequent calls will be fast.
